//! Reproduce the frozen study without requiring the integration CLI.
use anyhow::{ensure, Context, Result};
use genome_regulatory::{motif::MotifSet, run, InputRecord, Protocol};
use genome_regulatory_data::RegulatoryCohort;
use std::{fs, io::Write, path::Path};

fn main() -> Result<()> {
    let output = std::env::args()
        .nth(1)
        .context("usage: frozen_study NEW_OUTPUT_JSON")?;
    ensure!(!Path::new(&output).exists(), "output already exists");
    let protocol = Protocol::open(Path::new("experiments/regulatory-prediction-v1.json"))?;
    let motifs = MotifSet::open(Path::new("datasets/regulatory-motifs"))?;
    let cohort = RegulatoryCohort::open(Path::new("datasets/regulatory-cohort"))?;
    let records: Vec<_> = cohort
        .records()
        .iter()
        .map(|record| InputRecord {
            id: record.id.clone(),
            assay_id: record.assay_id.clone(),
            position_1based: record.target_position_1based,
            reference: record.reference.clone(),
            alternate: record.alternate.clone(),
            target_sequence: record.target_sequence.clone(),
            score: record.score,
            score_raw: record.score_raw.clone(),
            unique_tags: record.unique_tags,
            dna_sequences: record.dna_sequences,
            rna_sequences: record.rna_sequences,
            cell_context: record.cell_context.clone(),
        })
        .collect();
    eprintln!("Frozen protocol {}; {} supported source SNVs before quality selection; preparing fixed models and controls.",protocol.sha256(),records.len());
    let result = run(&records, &protocol, &motifs)?;
    let mut file = fs::OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&output)?;
    serde_json::to_writer_pretty(&mut file, &result)?;
    file.write_all(b"\n")?;
    file.sync_all()?;
    for analysis in &result.analyses {
        eprintln!("{}:", analysis.family);
        for model in &analysis.models {
            eprintln!(
                "  {} macro-position MSE {:.9}, MAE {:.9}",
                model.model, model.macro_assay_position_mse, model.macro_assay_position_mae
            );
        }
    }
    eprintln!("Saved complete deterministic result to {output}");
    Ok(())
}
