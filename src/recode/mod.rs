use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::sync::Arc;

use anyhow::{bail, Result};
use binseq::{
    bq::{BinseqHeaderBuilder, BinseqWriterBuilder},
    vbq::{VBinseqHeaderBuilder, VBinseqWriterBuilder},
    BitSize, Policy,
};
use ncbi_vdb_sys::SraReader;
use parking_lot::Mutex;

use crate::cli::{BinseqFlavor, RecodeArgs};
use crate::describe::describe_inner;
use crate::prefetch::identify_url;
use crate::utils::get_num_records;

const THREAD_UPDATE_INTERVAL: usize = 1024;

pub fn recode(args: &RecodeArgs) -> Result<()> {
    args.validate()?;
    let accession = if !Path::new(&args.input.accession).exists() {
        eprintln!(
            "Identifying SRA data URL for Accession: {}",
            &args.input.accession
        );
        let runtime = tokio::runtime::Runtime::new()?;
        let url = runtime.block_on(identify_url(&args.input.accession, &args.input.options))?;
        eprintln!("Streaming SRA records from URL: {}", url);
        url
    } else {
        args.input.accession.to_string()
    };

    match args.output.flavor {
        BinseqFlavor::Binseq => recode_to_binseq(
            &accession,
            &args.output.name(),
            args.primary_sid(),
            args.extended_sid(),
            args.output.bitsize(),
            args.runtime.threads(),
        ),
        BinseqFlavor::VBinseq => recode_to_vbinseq(
            &accession,
            &args.output.name(),
            args.primary_sid(),
            args.extended_sid(),
            args.output.bitsize(),
            args.output.block_size,
            args.runtime.threads(),
        ),
    }
}

fn recode_to_binseq(
    accession: &str,
    output_path: &str,
    primary_sid: usize,
    extended_sid: Option<usize>,
    bitsize: BitSize,
    num_threads: u64,
) -> Result<()> {
    let stats = describe_inner(accession, 0, 100)?;
    let sid_lengths = stats.segment_lengths();

    let slen = if sid_lengths[primary_sid].fract() == 0.0 {
        sid_lengths[primary_sid] as u32
    } else {
        bail!("Segment ID {primary_sid} shows variance in length. Cannot encode to BINSEQ (try VBINSEQ instead)")
    };

    let xlen = if let Some(extended_sid) = extended_sid {
        if sid_lengths[extended_sid].fract() == 0.0 {
            sid_lengths[extended_sid] as u32
        } else {
            bail!("Segment ID {extended_sid} shows variance in length. Cannot encode to BINSEQ (try VBINSEQ instead)")
        }
    } else {
        0
    };

    let output = File::create(output_path).map(BufWriter::new)?;
    let header = BinseqHeaderBuilder::new()
        .slen(slen)
        .xlen(xlen)
        .bitsize(bitsize)
        .flags(false)
        .build()?;
    let policy = Policy::RandomDraw;
    let g_writer = BinseqWriterBuilder::default()
        .header(header)
        .policy(policy)
        .build(output)?;
    let g_writer = Arc::new(Mutex::new(g_writer));

    let num_records = get_num_records(accession)?;
    let records_per_thread = num_records / num_threads;
    let remainder = num_records % num_threads;

    let mut handles = Vec::new();
    for tid in 0..num_threads {
        let start = (tid * records_per_thread) + 1; // 1-indexed
        let stop = if tid == num_threads - 1 {
            start + records_per_thread + remainder - 1
        } else {
            start + records_per_thread - 1
        };
        let t_accession = accession.to_string();
        let mut t_writer = BinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .policy(policy)
            .build(Vec::new())?;
        let g_writer = g_writer.clone();

        let handle = std::thread::spawn(move || -> Result<()> {
            let reader = SraReader::new(&t_accession)?;

            for (iter_index, record) in reader.into_range_iter(start as i64, stop)?.enumerate() {
                let record = record?;
                if xlen > 0 {
                    let primary_seg = record.get_segment(primary_sid).unwrap();
                    let extended_seg = record.get_segment(extended_sid.unwrap()).unwrap();
                    t_writer.write_paired_record(None, primary_seg.seq(), extended_seg.seq())?;
                } else {
                    let primary_seg = record.get_segment(primary_sid).unwrap();
                    t_writer.write_record(None, primary_seg.seq())?;
                }

                // Process records at a constant interval
                if iter_index % THREAD_UPDATE_INTERVAL == 0 {
                    {
                        let mut global = g_writer.lock();
                        global.ingest(&mut t_writer)?;
                        global.flush()?;
                    }
                }
            }

            // Process the last batch of records
            {
                let mut global = g_writer.lock();
                global.ingest(&mut t_writer)?;
                global.flush()?;
            }

            Ok(())
        });

        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap().unwrap();
    }

    Ok(())
}

fn recode_to_vbinseq(
    accession: &str,
    output_path: &str,
    primary_sid: usize,
    extended_sid: Option<usize>,
    bitsize: BitSize,
    block_size: usize,
    num_threads: u64,
) -> Result<()> {
    let output = File::create(output_path).map(BufWriter::new)?;
    let header = VBinseqHeaderBuilder::new()
        .block(block_size as u64)
        .flags(false)
        .qual(true)
        .headers(false)
        .bitsize(bitsize)
        .compressed(true)
        .paired(extended_sid.is_some())
        .build();

    let policy = Policy::RandomDraw;
    let g_writer = VBinseqWriterBuilder::default()
        .header(header)
        .policy(policy)
        .build(output)?;
    let g_writer = Arc::new(Mutex::new(g_writer));

    let num_records = get_num_records(accession)?;
    let records_per_thread = num_records / num_threads;
    let remainder = num_records % num_threads;

    let mut handles = Vec::new();
    for tid in 0..num_threads {
        let start = (tid * records_per_thread) + 1; // 1-indexed
        let stop = if tid == num_threads - 1 {
            start + records_per_thread + remainder - 1
        } else {
            start + records_per_thread - 1
        };
        let t_accession = accession.to_string();
        let mut t_writer = VBinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .policy(policy)
            .build(Vec::new())?;
        let g_writer = g_writer.clone();

        let handle = std::thread::spawn(move || -> Result<()> {
            let reader = SraReader::new(&t_accession)?;

            for (iter_index, record) in reader.into_range_iter(start as i64, stop)?.enumerate() {
                let record = record?;
                if let Some(extended_sid) = extended_sid {
                    let primary_seg = record.get_segment(primary_sid).unwrap();
                    let extended_seg = record.get_segment(extended_sid).unwrap();
                    t_writer.write_paired_record(
                        None,
                        None,
                        primary_seg.seq(),
                        Some(primary_seg.qual()),
                        None,
                        extended_seg.seq(),
                        Some(extended_seg.qual()),
                    )?;
                } else {
                    let primary_seg = record.get_segment(primary_sid).unwrap();
                    t_writer.write_record(
                        None,
                        None,
                        primary_seg.seq(),
                        Some(primary_seg.qual()),
                    )?;
                }

                // Process records at a constant interval
                if iter_index % THREAD_UPDATE_INTERVAL == 0 {
                    {
                        let mut global = g_writer.lock();
                        global.ingest(&mut t_writer)?;
                    }
                }
            }

            // Process the last batch of records
            {
                let mut global = g_writer.lock();
                global.ingest(&mut t_writer)?;
            }

            Ok(())
        });

        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap().unwrap();
    }

    g_writer.lock().finish()?;

    Ok(())
}
