use std::io::Write;

use anyhow::{bail, Result};
use ncbi_vdb_sys::Segment;

use crate::cli::OutputFormat;

pub fn write_segment_to_buffer_set(
    buffers: &mut [Vec<u8>],
    segment: &Segment<'_>,
    format: OutputFormat,
    accession_prefix: Option<&str>,
    include_sid: bool,
) -> Result<()> {
    if buffers.len() == 1 {
        // Interleaved output - single output handle
        let buffer = &mut buffers[0];
        match format {
            OutputFormat::Fasta => write_fasta(buffer, segment, accession_prefix, include_sid)?,
            OutputFormat::Fastq => write_fastq(buffer, segment, accession_prefix, include_sid)?,
        }
        Ok(())
    } else {
        if segment.sid() >= buffers.len() {
            bail!(
                "Provided Segment ID: {} is above the expected 4-segment counts",
                segment.sid()
            );
        }
        let seg_id = segment.sid();
        let buffer = &mut buffers[seg_id];
        match format {
            OutputFormat::Fasta => write_fasta(buffer, segment, accession_prefix, include_sid)?,
            OutputFormat::Fastq => write_fastq(buffer, segment, accession_prefix, include_sid)?,
        }
        Ok(())
    }
}

pub fn write_fastq<W: Write>(
    wtr: &mut W,
    segment: &Segment<'_>,
    accession_prefix: Option<&str>,
    include_sid: bool,
) -> Result<()> {
    if let Some(prefix) = accession_prefix {
        if include_sid {
            writeln!(wtr, "@{}.{}.{}", prefix, segment.rid(), segment.sid())
        } else {
            writeln!(wtr, "@{}.{}", prefix, segment.rid())
        }
    } else {
        writeln!(wtr, "@{}.{}", segment.rid(), segment.sid())
    }?;
    wtr.write_all(segment.seq())?;
    writeln!(wtr, "\n+")?;
    wtr.write_all(segment.qual())?;
    writeln!(wtr)?;
    Ok(())
}

pub fn write_fasta<W: Write>(
    wtr: &mut W,
    segment: &Segment<'_>,
    accession_prefix: Option<&str>,
    include_sid: bool,
) -> Result<()> {
    if let Some(prefix) = accession_prefix {
        if include_sid {
            writeln!(wtr, ">{}.{}.{}", prefix, segment.rid(), segment.sid())
        } else {
            writeln!(wtr, ">{}.{}", prefix, segment.rid())
        }
    } else {
        writeln!(wtr, ">{}.{}", segment.rid(), segment.sid())
    }?;
    wtr.write_all(segment.seq())?;
    writeln!(wtr)?;
    Ok(())
}
