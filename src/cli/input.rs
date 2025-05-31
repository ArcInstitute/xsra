use std::fmt::Display;

use clap::{Parser, ValueEnum};

#[derive(Debug, Parser)]
#[clap(next_help_heading = "INPUT OPTIONS")]
pub struct InputOptions {
    /// SRA accession or path to discrete SRA file or Directory
    #[clap(name = "SRA accession", required = true)]
    pub accession: String,

    /// Accession options
    #[clap(flatten)]
    pub options: AccessionOptions,
}

#[derive(Debug, Parser)]
#[clap(next_help_heading = "INPUT OPTIONS")]
pub struct MultiInputOptions {
    /// SRA accession or path to discrete SRA file or Directory
    #[clap(name = "SRA accession(s)", required = true, num_args = 1..)]
    pub accessions: Vec<String>,

    /// Accession options
    #[clap(flatten)]
    pub options: AccessionOptions,
}
impl MultiInputOptions {
    pub fn accession_set(&self) -> &[String] {
        &self.accessions
    }
}

#[derive(Debug, Clone, Parser)]
#[clap(next_help_heading = "ACCESSION OPTIONS")]
pub struct AccessionOptions {
    /// Only download an SRA with complete quality scores
    ///
    /// Default: lite, but will fallback to full quality if unavailable
    #[clap(short = 'Q', long)]
    pub full_quality: bool,

    /// URL provider
    #[clap(short = 'P', long, default_value = "https")]
    pub provider: Provider,

    /// GCP project ID
    #[clap(short = 'G', long, required_if_eq("provider", "gcp"))]
    pub gcp_project_id: Option<String>,

    /// Maximum number of retries on request limiting before bailing out
    #[clap(long, default_value = "5")]
    pub retry_limit: usize,

    /// Delay in milliseconds between retries
    #[clap(long, default_value = "500")]
    pub retry_delay: usize,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Provider {
    Https,
    Gcp,
    Aws,
}
impl Display for Provider {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Https => f.write_str("https"),
            Self::Gcp => f.write_str("gcp"),
            Self::Aws => f.write_str("aws"),
        }
    }
}
impl Provider {
    pub fn url_prefix(&self) -> &str {
        match self {
            Self::Https => "https://",
            Self::Gcp => "gs://",
            Self::Aws => "s3://",
        }
    }
}
