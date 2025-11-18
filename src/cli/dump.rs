use crate::output::Compression;
use clap::Parser;

use super::{FilterOptions, InputOptions, RuntimeOptions};

#[derive(Parser, Debug)]
pub struct DumpArgs {
    #[clap(flatten)]
    pub input: InputOptions,

    #[clap(flatten)]
    pub filter: FilterOptions,

    #[clap(flatten)]
    pub output: DumpOutput,

    #[clap(flatten)]
    pub runtime: RuntimeOptions,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "OUTPUT OPTIONS")]
pub struct DumpOutput {
    /// Output directory
    ///
    /// Only used when splitting read segments to separate files
    #[clap(short = 'o', long, default_value = "output")]
    pub outdir: String,

    /// Treat "output" as a filename stem and write output to a named pipe
    /// named as <outdir>.<prefix><segment>.<ext>
    ///
    /// Only used when splitting read segments to separate files
    #[clap(short = 'n', long, requires = "split")]
    pub named_pipes: bool,

    /// Output Format
    #[clap(short = 'f', long, default_value = "q")]
    pub format: OutputFormat,

    /// Split read segments to separate files
    ///
    /// Default will output interleaved reads to stdout
    #[clap(short = 's', long)]
    pub split: bool,

    /// Prefix for segment files
    ///
    /// Output will follow the pattern: <outdir>/<prefix><segment>.<ext>
    #[clap(short = 'p', long, default_value = "seg_")]
    pub prefix: String,

    /// Include the segment ID in split files
    #[clap(long, requires = "split")]
    pub include_sid: bool,

    /// Compress output files
    ///
    /// [uncompressed, gzip, bgzip, zstd]
    #[clap(short = 'c', long, default_value = "u")]
    pub compression: Compression,

    /// Keep empty files
    ///
    /// By default empty files will be deleted
    #[clap(short = 'E', long)]
    pub keep_empty: bool,
}
impl DumpOutput {
    pub fn include_sid(&self) -> bool {
        if self.split {
            self.include_sid
        } else {
            true
        }
    }
}

#[derive(clap::ValueEnum, Clone, Copy, Debug)]
pub enum OutputFormat {
    #[clap(name = "q", help = "FASTQ")]
    Fastq,
    #[clap(name = "a", help = "FASTA")]
    Fasta,
}
impl OutputFormat {
    pub fn ext(&self) -> &str {
        match self {
            Self::Fasta => "fa",
            Self::Fastq => "fq",
        }
    }
}
