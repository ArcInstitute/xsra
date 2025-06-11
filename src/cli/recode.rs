use super::{InputOptions, RuntimeOptions};
use anyhow::{bail, Result};
use clap::Parser;

#[derive(Parser, Debug)]
pub struct RecodeArgs {
    #[clap(flatten)]
    pub input: InputOptions,

    #[clap(flatten)]
    pub selection: SelectionOptions,

    #[clap(flatten)]
    pub runtime: RuntimeOptions,

    #[clap(flatten)]
    pub output: RecodeOutput,
}
impl RecodeArgs {
    pub fn validate(&self) -> Result<()> {
        match &self.selection.include.len() {
            0 => bail!(
                "Recoding requires including at least one spot segment (see 'xsra recode --help' for usage)"
            ),
            1 | 2 => Ok(()),
            _ => bail!("Recoding can only include one or two spot segments"),
        }
    }

    pub fn paired(&self) -> bool {
        self.selection.include.len() == 2
    }

    pub fn primary_sid(&self) -> usize {
        self.selection.include[0]
    }

    pub fn extended_sid(&self) -> Option<usize> {
        if self.paired() {
            Some(self.selection.include[1])
        } else {
            None
        }
    }
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "SELECTION OPTIONS")]
pub struct SelectionOptions {
    /// Only process up to N spots
    ///
    /// Note: This is not the number of individual segments, but rather the number of spots (or records) to process.
    #[clap(short = 'l', long)]
    pub limit: Option<usize>,

    /// Include specific segments (zero-indexed) as CSV
    ///
    /// I.e. to include the first and third segments, use "-I 0,2".
    ///
    /// The first entry is the primary spot segment, and the second entry is the extended spot segment.
    #[clap(short = 'I', long, num_args = 1..=2, value_delimiter = ',')]
    pub include: Vec<usize>,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "OUTPUT OPTIONS")]
pub struct RecodeOutput {
    /// BINSEQ output name (default: "output.{bq,vbq}")
    #[clap(short, long)]
    pub name: Option<String>,

    /// BINSEQ output flavor
    #[clap(short, long)]
    pub flavor: BinseqFlavor,

    /// VBQ virtual block size (in bytes)
    ///
    /// Only used by vbq
    #[clap(short = 'B', long, value_parser = parse_memory_size, default_value = "128K")]
    pub block_size: usize,
}
impl RecodeOutput {
    pub fn name(&self) -> String {
        if let Some(name) = &self.name {
            name.clone()
        } else {
            let ext = self.flavor.extension();
            format!("output.{}", ext)
        }
    }
}

#[derive(clap::ValueEnum, Clone, Copy, Debug)]
pub enum BinseqFlavor {
    #[clap(name = "b", help = "BINSEQ")]
    Binseq,
    #[clap(name = "v", help = "VBINSEQ")]
    VBinseq,
}
impl BinseqFlavor {
    pub fn extension(&self) -> &str {
        match self {
            BinseqFlavor::Binseq => "bq",
            BinseqFlavor::VBinseq => "vbq",
        }
    }
}

fn parse_memory_size(input: &str) -> Result<usize, String> {
    let input = input.trim().to_uppercase();
    let last_char = input.chars().last().unwrap_or('0');

    let (number_str, multiplier) = match last_char {
        'K' => (&input[..input.len() - 1], 1024),
        'M' => (&input[..input.len() - 1], 1024 * 1024),
        'G' => (&input[..input.len() - 1], 1024 * 1024 * 1024),
        _ if last_char.is_ascii_digit() => (input.as_str(), 1),
        _ => return Err(format!("Invalid memory size format: {input}")),
    };

    match number_str.parse::<usize>() {
        Ok(number) => Ok(number * multiplier),
        Err(_) => Err(format!("Failed to parse number: {number_str}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // parse_memory_size tests
    #[test]
    fn parse_memory_size_k_suffix() {
        assert_eq!(parse_memory_size("1K"), Ok(1024));
    }

    #[test]
    fn parse_memory_size_m_suffix() {
        assert_eq!(parse_memory_size("1M"), Ok(1024 * 1024));
    }

    #[test]
    fn parse_memory_size_g_suffix() {
        assert_eq!(parse_memory_size("1G"), Ok(1024 * 1024 * 1024));
    }

    #[test]
    fn parse_memory_size_ascii_digit() {
        assert_eq!(parse_memory_size("1024"), Ok(1024));
    }

    #[test]
    fn parse_memory_size_invalid_format() {
        assert!(parse_memory_size("512X").is_err());
    }

    #[test]
    fn parse_memory_size_invalid_number() {
        assert!(parse_memory_size("abcK").is_err());
    }
}
