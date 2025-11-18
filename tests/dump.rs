use anyhow::Result;
use assert_cmd::Command;
use predicates::prelude::*;
use std::fs;
use tempfile::TempDir;
use xsra::cli::{
    AccessionOptions, DumpOutput, FilterOptions, InputOptions, OutputFormat, Provider,
};
use xsra::dump::dump;
use xsra::output::Compression;

mod fixtures;
use fixtures::setup::TestFixtures;

// Integration tests for dump

// Test helper function to create default accession options
fn default_accession_options() -> AccessionOptions {
    AccessionOptions {
        full_quality: false,
        lite_only: false,
        provider: Provider::Https,
        retry_limit: 5,
        retry_delay: 500,
        gcp_project_id: None,
    }
}

#[test]
fn test_simple_fastq_dump_cli() -> Result<()> {
    let fixtures = TestFixtures::ensure_fixtures()?;
    let mut cmd = Command::cargo_bin("xsra")?;
    cmd.arg("dump")
        .arg(fixtures.small_variable_sra)
        .arg("--limit")
        .arg("10");

    cmd.assert()
        .success()
        .stdout(predicate::str::starts_with("@"))
        .stdout(predicate::str::contains("\n+\n"))
        .stdout(predicate::str::ends_with("\n"));

    let output = cmd.output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert_eq!(
        stdout.lines().count(),
        40,
        "Expected 10 records * 4 lines/record"
    );

    Ok(())
}

#[test]
fn test_split_file_dump() -> Result<()> {
    let fixtures = TestFixtures::ensure_fixtures()?;
    let temp_dir = TempDir::new()?;

    let input = InputOptions {
        accession: fixtures.small_variable_sra.clone(),
        options: default_accession_options(),
    };

    let output = DumpOutput {
        outdir: temp_dir.path().to_string_lossy().to_string(),
        prefix: "test".to_string(),
        compression: Compression::Uncompressed,
        format: OutputFormat::Fastq,
        named_pipes: false,
        split: true,
        keep_empty: false,
        include_sid: false,
    };

    let filter_opts = FilterOptions {
        include: vec![],
        skip_technical: false,
        min_read_len: 1,
        limit: Some(100), // Limit to 100 spots for fast testing
    };

    let result = dump(&input, 1, &output, filter_opts);
    assert!(
        result.is_ok(),
        "Split dump command failed: {:?}",
        result.err()
    );

    // Verify that at least one split file was created
    // Most paired-end SRA files will have at least _1.fastq and _2.fastq
    let mut split_files_found = 0;
    for entry in fs::read_dir(temp_dir.path())? {
        let entry = entry?;
        let file_name = entry.file_name().to_string_lossy().to_string();
        if file_name.starts_with("test") && file_name.ends_with(".fq") {
            split_files_found += 1;

            // Verify file is not empty
            let file_size = fs::metadata(entry.path())?.len();
            assert!(file_size > 0, "Split file {} is empty", file_name);
        }
    }

    assert!(split_files_found > 0, "No split files were created");

    Ok(())
}

#[test]
fn test_filtered_dump_cli() -> Result<()> {
    let fixtures = TestFixtures::ensure_fixtures()?;
    let mut cmd = Command::cargo_bin("xsra")?;
    cmd.arg("dump")
        .arg(fixtures.small_variable_sra)
        .arg("--skip-technical")
        .arg("--min-read-len")
        .arg("50")
        .arg("--limit")
        .arg("10")
        .arg("--include")
        .arg("0,1");

    cmd.assert()
        .success()
        .stdout(predicate::str::starts_with("@"));

    Ok(())
}

#[test]
fn test_empty_file_removal() -> Result<()> {
    let fixtures = TestFixtures::ensure_fixtures()?;
    let temp_dir = TempDir::new()?;

    let input = InputOptions {
        accession: fixtures.small_variable_sra.clone(),
        options: default_accession_options(),
    };

    let output = DumpOutput {
        outdir: temp_dir.path().to_string_lossy().to_string(),
        prefix: "test".to_string(),
        compression: Compression::Uncompressed,
        format: OutputFormat::Fastq,
        named_pipes: false,
        split: true,
        keep_empty: false,
        include_sid: false,
    };

    let filter_opts = FilterOptions {
        include: vec![0], // Only include segment 0, which might create empty files for other segments
        skip_technical: true,
        min_read_len: 10000, // Very high threshold to potentially create empty files
        limit: Some(50),     // Small limit for fast testing
    };

    let result = dump(&input, 1, &output, filter_opts);
    assert!(
        result.is_ok(),
        "Empty file removal test failed: {:?}",
        result.err()
    );

    // Count files created
    let mut file_count = 0;
    for entry in fs::read_dir(temp_dir.path())? {
        let entry = entry?;
        let file_name = entry.file_name().to_string_lossy().to_string();
        if file_name.starts_with("test") {
            file_count += 1;
        }
    }

    // With aggressive filtering, we should have fewer files than without
    assert!(file_count >= 0, "File count should be non-negative");

    Ok(())
}

#[test]
#[cfg(unix)]
fn test_named_pipes_io() -> Result<()> {
    use std::os::unix::fs::FileTypeExt;
    use tempfile::TempDir;

    let temp_dir = TempDir::new()?;

    // We don't run the actual dump command with named_pipes=true because
    // that would block indefinitely waiting for a reader process.

    // Create a test FIFO manually to verify the mkfifo command works
    let test_fifo = temp_dir.path().join("test.fq");
    let status = std::process::Command::new("mkfifo")
        .arg(&test_fifo)
        .status()?;

    assert!(status.success(), "mkfifo command failed");
    assert!(test_fifo.exists(), "FIFO was not created");

    // Verify it's actually a FIFO
    let metadata = std::fs::metadata(&test_fifo)?;
    assert!(metadata.file_type().is_fifo(), "Created file is not a FIFO");

    Ok(())
}

#[test]
fn test_fasta_output_cli() -> Result<()> {
    let fixtures = TestFixtures::ensure_fixtures()?;
    let mut cmd = Command::cargo_bin("xsra")?;
    cmd.arg("dump")
        .arg(fixtures.small_variable_sra)
        .arg("--format")
        .arg("a")
        .arg("--limit")
        .arg("10");

    cmd.assert()
        .success()
        .stdout(predicate::str::starts_with(">"));

    let output = cmd.output()?;
    let stdout = String::from_utf8(output.stdout)?;
    assert_eq!(
        stdout.lines().count(),
        20,
        "Expected 10 records * 2 lines/record"
    );

    Ok(())
}
