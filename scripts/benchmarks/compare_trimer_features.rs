//! Release comparison for two trimer feature implementations on one frozen panel.
use anyhow::{ensure, Result};
use genome_context_model::{feature_names, sequence_features};
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::{fs, hint::black_box, time::Instant};

#[allow(dead_code)]
#[path = "baseline/mod.rs"]
mod baseline;

const PANEL_SHA256: &str = "fce6b2a1aaf5aa4d530a2e298dccb840bc7c594a05b6831fdb0c6c0f9e2a9be3";

#[derive(Deserialize)]
struct Panel {
    records: Vec<Record>,
}

#[derive(Deserialize)]
struct Record {
    sequence: String,
}

fn time<T>(records: &[Record], kernel: impl Fn(&str) -> Result<T>) -> Result<f64> {
    let start = Instant::now();
    for record in records {
        black_box(kernel(black_box(record.sequence.as_str()))?);
    }
    Ok(start.elapsed().as_secs_f64())
}

fn main() -> Result<()> {
    let bytes = fs::read("datasets/context-panel/panel.json")?;
    ensure!(
        format!("{:x}", Sha256::digest(&bytes)) == PANEL_SHA256,
        "context panel SHA-256 mismatch"
    );
    let Panel { records } = serde_json::from_slice(&bytes)?;
    assert_eq!(records.len(), 55_276);
    assert_eq!(feature_names(), baseline::feature_names());
    for (index, record) in records.iter().enumerate() {
        assert_eq!(record.sequence.len(), 200);
        let old = baseline::sequence_features(&record.sequence)?;
        let new = sequence_features(&record.sequence)?;
        assert_eq!(old.length, new.length);
        assert_eq!(old.canonical_trimer_counts, new.canonical_trimer_counts);
        ensure!(
            old.vector()
                .iter()
                .map(|value| value.to_bits())
                .eq(new.vector().iter().map(|value| value.to_bits())),
            "feature bits differ at panel row {index}"
        );
    }

    let mut old_seconds = Vec::new();
    let mut new_seconds = Vec::new();
    for pass in 0..5 {
        if pass % 2 == 0 {
            old_seconds.push(time(&records, baseline::sequence_features)?);
            new_seconds.push(time(&records, sequence_features)?);
        } else {
            new_seconds.push(time(&records, sequence_features)?);
            old_seconds.push(time(&records, baseline::sequence_features)?);
        }
    }
    old_seconds.sort_by(f64::total_cmp);
    new_seconds.sort_by(f64::total_cmp);
    let receipt = serde_json::json!({
        "panel_sha256": PANEL_SHA256,
        "equal_panel_records": records.len(),
        "baseline_seconds": old_seconds,
        "packed_seconds": new_seconds,
        "median_speedup": old_seconds[2] / new_seconds[2],
    });
    println!("{}", serde_json::to_string_pretty(&receipt)?);
    Ok(())
}
