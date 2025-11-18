mod output;
mod stats;
mod utils;

use std::path::Path;
use std::sync::Arc;

use anyhow::Result;
use ncbi_vdb_sys::SraReader;
use output::{build_segment_writer, BoxedSegmentWriter};
use parking_lot::Mutex;

use crate::cli::{DumpOutput, FilterOptions, InputOptions, OutputFormat};
use crate::output::{build_path_name, OutputFileType};
use crate::prefetch::identify_url;
use crate::RECORD_CAPACITY;

use crate::utils::get_num_records;
use stats::ProcessStatistics;
use utils::write_segment_to_buffer_set;

fn launch_threads(
    path: &str,
    num_threads: u64,
    records_per_thread: u64,
    remainder: u64,
    writer: Arc<Mutex<BoxedSegmentWriter>>,
    filter_opts: FilterOptions,
    format: OutputFormat,
    accession_prefix: Option<String>,
    include_sid: bool,
) -> Result<ProcessStatistics> {
    // Segments included in the output
    let segment_set = if filter_opts.include.is_empty() {
        None
    } else {
        // checking a small vector should be faster than a HashSet
        let set: Vec<usize> = filter_opts.include.clone();
        Some(set)
    };

    let mut handles = Vec::new();
    for i in 0..num_threads {
        let segment_set = segment_set.clone();

        let start = (i * records_per_thread) + 1;
        let stop = if i == num_threads - 1 {
            start + records_per_thread + remainder - 1
        } else {
            start + records_per_thread - 1
        };
        let path = path.to_string();
        let shared_writer = writer.clone();
        let accession_prefix = accession_prefix.clone();

        let handle = std::thread::spawn(move || -> Result<ProcessStatistics> {
            let reader = SraReader::new(&path)?;

            // Initialize local buffers and counters
            let mut stats = ProcessStatistics::default();
            let mut local_buffers = shared_writer.lock().generate_local_buffers();
            let mut counts = vec![0; local_buffers.len()];

            // Iterate over record spots and write to buffers
            for (idx, record) in reader.into_range_iter(start as i64, stop)?.enumerate() {
                let record = record?;

                // Iterate over segments in the record
                for segment in record.into_iter() {
                    // Skip segment if outside of set
                    if let Some(ref set) = segment_set {
                        if !set.contains(&segment.sid()) {
                            continue;
                        }
                    }

                    // Skip technical segments if required
                    if filter_opts.skip_technical && segment.is_technical() {
                        // Increment filter statistics
                        stats.inc_filter_type(segment.sid());
                        continue;
                    }

                    // Skip reads if they are under the minimum read length
                    if segment.len() < filter_opts.min_read_len {
                        // Increment filter statistics
                        stats.inc_filter_size(segment.sid());
                        continue;
                    }

                    // Write the segment to the record set
                    write_segment_to_buffer_set(
                        &mut local_buffers,
                        &segment,
                        format,
                        accession_prefix.as_deref(),
                        include_sid,
                    )?;

                    if counts.len() == 1 {
                        counts[0] += 1;
                    } else {
                        counts[segment.sid()] += 1;
                    }

                    // Increment read statistics
                    stats.inc_reads(segment.sid());
                }

                // Handle buffer writes at specific intervals
                if idx > 0 && (idx % RECORD_CAPACITY == 0) {
                    shared_writer
                        .lock()
                        .write_all_buffers(&mut local_buffers, &mut counts)?;
                }

                // Increment record statistics
                stats.inc_spots();
            }

            // write remaining buffers
            shared_writer
                .lock()
                .write_all_buffers(&mut local_buffers, &mut counts)?;

            // Return thread-specific statistics
            Ok(stats)
        });

        // Collect all thread handles
        handles.push(handle);
    }

    // Collect all statistics
    let mut stats = ProcessStatistics::default();
    for handle in handles {
        let thread_stats = handle.join().expect("Thread panicked")?;
        stats = stats + thread_stats;
    }

    Ok(stats)
}

pub fn dump(
    input: &InputOptions,
    num_threads: u64,
    output_opts: &DumpOutput,
    filter_opts: FilterOptions,
) -> Result<()> {
    let accession = if !Path::new(&input.accession).exists() {
        eprintln!(
            "Identifying SRA data URL for Accession: {}",
            &input.accession
        );
        let runtime = tokio::runtime::Runtime::new()?;
        let url = runtime.block_on(identify_url(&input.accession, &input.options))?;
        eprintln!("Streaming SRA records from URL: {}", url);
        url
    } else {
        input.accession.to_string()
    };

    let num_records = get_num_records(&accession)?;

    // Adjust the number of records to process if a limit is provided
    let num_records = if let Some(limit) = filter_opts.limit {
        if limit > num_records {
            eprintln!("Warning: Provided spot limit ({}) is greater than the actual number of spots ({}). Will process the full archive.",
                limit, num_records);
        }
        num_records.min(limit)
    } else {
        num_records
    };

    // Calculate records per thread and remainder (final thread)
    let records_per_thread = num_records / num_threads;
    let remainder = num_records % num_threads;

    let writer = build_segment_writer(
        Some(&output_opts.outdir),
        &output_opts.prefix,
        output_opts.compression,
        output_opts.format,
        num_threads as usize,
        &filter_opts,
        output_opts.named_pipes,
        output_opts.split,
    )
    .map(|x| Arc::new(Mutex::new(x)))?;

    let included_segs = filter_opts.include.clone();
    let accession_name = Path::new(&input.accession)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(&input.accession)
        .to_string();
    let accession_prefix = Some(accession_name);

    // Launch worker threads
    let stats = launch_threads(
        &accession,
        num_threads,
        records_per_thread,
        remainder,
        writer,
        filter_opts,
        output_opts.format,
        accession_prefix,
        output_opts.include_sid(),
    )?;

    // Remove empty files
    if output_opts.split {
        let wrap = |x| {
            if output_opts.named_pipes {
                OutputFileType::NamedPipe(x)
            } else {
                OutputFileType::RegularFile(x)
            }
        };
        stats.reads_per_segment.iter().enumerate().try_for_each(
            |(seg_id, &count)| -> Result<()> {
                // if included_segs was non-empty, so we are applying a filter
                // and this segment was not included, then we didn't create a real
                // file for it.
                if !included_segs.is_empty() && !included_segs.contains(&seg_id) {
                    return Ok(());
                }
                if count == 0 || output_opts.named_pipes {
                    let path = build_path_name(
                        wrap(&output_opts.outdir),
                        &output_opts.prefix,
                        output_opts.compression,
                        output_opts.format,
                        seg_id,
                    );
                    if output_opts.keep_empty {
                        eprintln!("Warning => empty path: {}", path);
                    } else {
                        eprintln!("Removing empty path: {}", path);
                        std::fs::remove_file(path)?;
                    }
                }
                Ok(())
            },
        )?;
    }

    // Print all statistics
    stats.pprint(&mut std::io::stderr())?;

    Ok(())
}
