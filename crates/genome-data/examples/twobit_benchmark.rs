//! Replay the pinned hg38 interval benchmark and independent UCSC API checks.
//! Run with `cargo run --release --locked -p genome-data --example twobit_benchmark -- hg38.2bit`.
use anyhow::{ensure, Context, Result};
use genome_data::{
    sha256,
    twobit::{Strand, TwoBitReference},
};
use serde::Deserialize;
use serde_json::json;
use sha2::{Digest, Sha256};
use std::{env, path::PathBuf, time::Instant};

#[derive(Deserialize)]
struct ReferenceReceipt {
    bytes: u64,
    sha256: String,
}

#[derive(Deserialize)]
struct ValidationReceipt {
    ucsc_api_checks: Vec<SequenceCheck>,
    sequence_only_performance: QueryBench,
}

#[derive(Deserialize)]
struct SequenceCheck {
    chromosome: String,
    start0: u64,
    end0: u64,
    sequence_sha256_uppercase: String,
    sequence_sha256_softmasked: String,
}

#[derive(Deserialize)]
struct QueryBench {
    query_count: u64,
    bases_per_query: u64,
    output_sha256: String,
}

fn main() -> Result<()> {
    let path = env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .context("pass the pinned hg38.2bit file path after --")?;
    let identity: ReferenceReceipt = serde_json::from_str(include_str!(
        "../../../research/context-reference/hg38-reference.json"
    ))?;
    let validation: ValidationReceipt = serde_json::from_str(include_str!(
        "../../../research/context-reference/validation.json"
    ))?;
    let opened = Instant::now();
    let mut reference = TwoBitReference::open(&path)?;
    let open_seconds = opened.elapsed().as_secs_f64();
    ensure!(
        reference.file_bytes() == identity.bytes,
        "reference size differs from the pinned hg38 file"
    );
    let verifying = Instant::now();
    let reference_sha256 = reference.sha256()?;
    ensure!(
        reference_sha256 == identity.sha256,
        "reference checksum differs from the pinned hg38 file"
    );
    let reference_sha256_seconds = verifying.elapsed().as_secs_f64();

    for check in &validation.ucsc_api_checks {
        for (masked, expected) in [
            (false, &check.sequence_sha256_uppercase),
            (true, &check.sequence_sha256_softmasked),
        ] {
            let sequence = reference.fetch(
                &check.chromosome,
                check.start0,
                check.end0,
                Strand::Forward,
                masked,
            )?;
            ensure!(
                sha256(sequence.as_bytes()) == *expected,
                "UCSC API check failed for {}:{}-{} with masked={masked}",
                check.chromosome,
                check.start0,
                check.end0
            );
        }
    }

    let bench = validation.sequence_only_performance;
    let chr1_length = reference
        .contigs()
        .find(|contig| contig.name == "chr1")
        .context("chr1 missing from reference index")?
        .length;
    let start_floor = 10_000;
    let span = chr1_length
        .checked_sub(start_floor + bench.bases_per_query)
        .filter(|span| *span > 0)
        .context("chr1 is too short for the retained queries")?;
    let started = Instant::now();
    let mut digest = Sha256::new();
    for index in 0..bench.query_count {
        let start = start_floor + (index * 104_729) % span;
        let end = start + bench.bases_per_query;
        let sequence = reference.fetch("chr1", start, end, Strand::Forward, false)?;
        digest.update(sequence.as_bytes());
    }
    let query_seconds = started.elapsed().as_secs_f64();
    let output_sha256 = format!("{:x}", digest.finalize());
    ensure!(
        output_sha256 == bench.output_sha256,
        "query digest disagrees with the retained Python oracle"
    );
    let receipt = json!({
        "reference_sha256": reference_sha256,
        "contigs": reference.contigs().count(),
        "ucsc_windows_matched": validation.ucsc_api_checks.len(),
        "open_seconds": open_seconds,
        "reference_sha256_seconds": reference_sha256_seconds,
        "query_count": bench.query_count,
        "decoded_bases": bench.query_count * bench.bases_per_query,
        "query_seconds_including_output_sha256": query_seconds,
        "output_sha256": output_sha256,
    });
    println!("{}", serde_json::to_string_pretty(&receipt)?);
    Ok(())
}
