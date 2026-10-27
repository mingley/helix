//! Frozen regulatory study orchestration and immutable compressed artifacts.
use anyhow::{ensure, Context, Result};
use flate2::{read::GzDecoder, Compression, GzBuilder};
use genome_knowledge::{Ledger, Revision};
use genome_logic::{Derivation, Rule};
use genome_regulatory::{motif::MotifSet, InputRecord, Protocol, StudyResult};
use genome_regulatory_data::RegulatoryCohort;
use genome_types::Probability;
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::{Read, Write},
    path::{Path, PathBuf},
    process::Command,
};

#[derive(clap::Args)]
pub struct StudyArgs {
    #[arg(long, default_value = "experiments/regulatory-prediction-v1.json")]
    protocol: PathBuf,
    #[arg(long, default_value = "datasets/regulatory-cohort")]
    cohort: PathBuf,
    #[arg(long, default_value = "datasets/regulatory-motifs")]
    motifs: PathBuf,
    #[arg(long, default_value = "research/regulatory-forensics.json")]
    forensics: PathBuf,
    #[arg(long, default_value = "lean")]
    lean: PathBuf,
    #[arg(long)]
    out: PathBuf,
}

const FILES: [&str; 12] = [
    "cohort.json",
    "context.json",
    "findings.md",
    "forensics.json",
    "inputs.json",
    "ledger.json",
    "motifs.json",
    "proof.lean",
    "protocol.json",
    "resampling.json",
    "study.json.gz",
    "summary.json",
];

fn source_hashes() -> Result<BTreeMap<String, String>> {
    fn walk(path: &Path, out: &mut BTreeMap<String, String>) -> Result<()> {
        for entry in fs::read_dir(path)? {
            let entry = entry?;
            let p = entry.path();
            if entry.file_type()?.is_dir() {
                walk(&p, out)?;
            } else if p.extension().is_some_and(|s| s == "rs")
                || p.file_name().is_some_and(|s| s == "Cargo.toml")
            {
                out.insert(
                    p.to_string_lossy().into_owned(),
                    genome_data::sha256(&fs::read(p)?),
                );
            }
        }
        Ok(())
    }
    let mut out = BTreeMap::new();
    walk(Path::new("crates"), &mut out)?;
    for name in [
        "Cargo.lock",
        "Cargo.toml",
        "rust-toolchain.toml",
        "proofs/RegulatoryStudy.lean",
    ] {
        out.insert(name.into(), genome_data::sha256(&fs::read(name)?));
    }
    Ok(out)
}

fn summarize(result: &StudyResult) -> Value {
    let analyses:Vec<_>=result.analyses.iter().map(|a|json!({
        "family":a.family,"interpretation":a.interpretation,"fold_count":a.folds.len(),
        "fold_diagnostics":a.folds,
        "models":a.models.iter().map(|m|json!({
            "model":m.model,"feature_count":m.feature_names.len(),
            "prediction_count":m.predictions.len(),
            "macro_assay_position_mse":m.macro_assay_position_mse,
            "macro_assay_position_mae":m.macro_assay_position_mae,
            "pooled_row_metrics":m.pooled_row_metrics,"per_assay_metrics":m.per_assay_metrics
        })).collect::<Vec<_>>()
    })).collect();
    json!({"schema_version":1,"protocol_sha256":result.protocol_sha256,
        "input_records_sha256":result.input_records_sha256,"motif_manifest_sha256":result.motif_manifest_sha256,
        "inclusion":result.inclusion,"feature_schema":result.feature_schema,"metric_arithmetic":result.metric_arithmetic,"analyses":analyses,
        "limitations":result.limitations,
        "replication_policy":"Independent numerical comparison is a separate required verification bundle. This artifact alone does not claim reproduction or biological novelty.",
        "secondary_rank_limit":"Exact-tie Spearman can be sensitive to floating reduction order. Preserve attempted comparisons and distinguish replicated continuous losses from unresolved rank endpoints."})
}

pub fn run(args: &StudyArgs) -> Result<()> {
    ensure!(
        !args.out.exists(),
        "output already exists; use a new study directory"
    );
    let protocol = Protocol::open(&args.protocol)?;
    let cohort = RegulatoryCohort::open(&args.cohort)?;
    let motifs = MotifSet::open(&args.motifs)?;
    ensure!(
        cohort.manifest.frozen_protocol_sha256 == protocol.sha256(),
        "cohort gate differs from study protocol"
    );
    let forensic_bytes = fs::read(&args.forensics)?;
    let forensics: Value = serde_json::from_slice(&forensic_bytes)?;
    let inputs = json!({
        "protocol_sha256":protocol.sha256(),"cohort_manifest_sha256":cohort.manifest_sha256(),
        "motif_manifest_sha256":motifs.manifest_sha256(),
        "forensics_sha256":genome_data::sha256(&forensic_bytes),
        "implementation_files":source_hashes()?,
        "software":{
            "version":env!("CARGO_PKG_VERSION"),
            "commit":Command::new("git").args(["rev-parse","HEAD"]).output().ok().filter(|o|o.status.success()).map(|o|String::from_utf8_lossy(&o.stdout).trim().to_owned()),
            "working_tree_dirty":Command::new("git").args(["status","--porcelain"]).output().ok().map(|o|!o.status.success()||!o.stdout.is_empty())
        }
    });
    let records: Vec<_> = cohort
        .records()
        .iter()
        .map(|r| InputRecord {
            id: r.id.clone(),
            assay_id: r.assay_id.clone(),
            position_1based: r.target_position_1based,
            reference: r.reference.clone(),
            alternate: r.alternate.clone(),
            target_sequence: r.target_sequence.clone(),
            score: r.score,
            score_raw: r.score_raw.clone(),
            unique_tags: r.unique_tags,
            dna_sequences: r.dna_sequences,
            rna_sequences: r.rna_sequences,
            cell_context: r.cell_context.clone(),
        })
        .collect();
    let mut ledger = Ledger::default();
    ledger.append("freeze_protocol_and_inputs", &inputs)?;
    ledger.append("validate_metadata_selected_cohort", cohort.report())?;
    ledger.append("validate_prior_motif_assets", motifs.report())?;
    eprintln!("Executing six fixed models, four evaluations, 99 rotation controls and 2000 bootstrap draws.");
    let result = genome_regulatory::run(&records, &protocol, &motifs)?;
    let context = super::study_context::analyze(&cohort, protocol.as_value())?;
    let declared: BTreeSet<_> = result.inclusion["per_assay"]
        .as_array()
        .context("missing inclusion inventory")?
        .iter()
        .flat_map(|a| a["included_ids"].as_array().into_iter().flatten())
        .filter_map(Value::as_str)
        .collect();
    let context_ids: BTreeSet<_> = cohort
        .records()
        .iter()
        .filter(|r| super::study_context::quality_passes(r))
        .map(|r| r.id.as_str())
        .collect();
    ensure!(
        declared == context_ids,
        "model and context analysis eligibility differ"
    );
    let summary = summarize(&result);
    let resampling = json!({"bootstrap":result.resampling,"spatial_null":result.spatial_null});
    ledger.append("execute_frozen_evaluation", &summary)?;
    ledger.append("execute_fixed_resampling",json!({"sha256":genome_data::sha256(&super::pretty(&resampling)?),
        "spatial_rotations":99,"bootstrap_draws":2000,"interpretation":"conditional robustness diagnostics, not calibrated biological significance"}))?;
    ledger.append("compare_source_contexts", &context)?;
    ledger.append("retain_source_forensics", &forensics)?;
    let derivation=Derivation::derive(
        vec!["The committed source observations and metadata are the intended study inputs".into()],
        vec![Rule{
            premise:"The committed source observations and metadata are the intended study inputs".into(),
            conclusion:"The fixed study describes conditional prediction error and source-label consistency, not a proven biological mechanism".into(),
            assumption:"The implementations, source estimates and experimental-context descriptions satisfy their documented contracts; independent numerical and empirical checks remain distinct".into()
        }],
        "The fixed study describes conditional prediction error and source-label consistency, not a proven biological mechanism".into())?;
    let n = result.inclusion["included_rows"]
        .as_u64()
        .context("missing included rows")?;
    let ratio = Probability::new(u32::try_from(n)?, u32::try_from(records.len())?)?;
    let proof_source = genome_proof::regulatory_source(
        &derivation,
        &genome_logic::probability_obligations(&[ratio]),
    );
    let proof = genome_proof::check_regulatory(&args.lean, &proof_source)?;
    ledger.append("check_scoped_method_contracts", &proof)?;
    ledger.append("record_inference_limits",json!({
        "new_biological_discovery":false,
        "source_studies":1,
        "rank_endpoint_reproduction":"must be established separately; retain all failures and numerical corrections",
        "novelty_review":"requires the literature/claim ledger; a favorable robustness tail alone is insufficient"
    }))?;
    Ledger::verify(ledger.revisions())?;
    let findings = narrative(&summary, &context)?;
    if let Some(parent) = args.out.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    fs::create_dir(&args.out)?;
    let mut commitments = BTreeMap::new();
    let file = fs::OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(args.out.join("study.json.gz"))?;
    let mut gzip = GzBuilder::new()
        .mtime(0)
        .write(file, Compression::default());
    serde_json::to_writer(&mut gzip, &result)?;
    gzip.write_all(b"\n")?;
    gzip.finish()?.sync_all()?;
    commitments.insert(
        "study.json.gz".into(),
        genome_data::sha256(&fs::read(args.out.join("study.json.gz"))?),
    );
    for (name, bytes) in [
        ("cohort.json", super::pretty(cohort.report())?),
        ("context.json", super::pretty(&context)?),
        ("findings.md", findings.into_bytes()),
        ("forensics.json", forensic_bytes),
        ("inputs.json", super::pretty(&inputs)?),
        ("ledger.json", super::pretty(&ledger.revisions())?),
        ("motifs.json", super::pretty(&motifs.report())?),
        ("proof.lean", proof_source.into_bytes()),
        ("protocol.json", fs::read(&args.protocol)?),
        ("resampling.json", super::pretty(&resampling)?),
        ("summary.json", super::pretty(&summary)?),
    ] {
        super::write_new(&args.out.join(name), &bytes)?;
        commitments.insert(name.into(), genome_data::sha256(&bytes));
    }
    super::write_new(
        &args.out.join("artifacts.json"),
        &super::pretty(&super::Artifacts {
            schema_version: 3,
            files: commitments,
        })?,
    )?;
    for family in &result.analyses {
        println!("{}", family.family);
        for model in &family.models {
            println!(
                "  {}: macro position MSE {:.9}",
                model.model, model.macro_assay_position_mse
            );
        }
    }
    println!(
        "Saved complete compressed study and source-context diagnostics to {}",
        args.out.display()
    );
    println!("Independent Python comparison and novelty review are required before reproduction or novelty claims.");
    Ok(())
}

fn narrative(summary: &Value, context: &Value) -> Result<String> {
    let mut text=String::from("# Frozen regulatory study\n\nThis is a retrospective computational reanalysis of published reporter-effect estimates. It does not establish new biological mechanisms, independent laboratory replication, clinical risk, or whole-genome understanding.\n\n## Fixed evaluation results\n\n| Evaluation | Model | Macro position MSE | Macro position MAE |\n| --- | --- | ---: | ---: |\n");
    for analysis in summary["analyses"].as_array().context("missing analyses")? {
        for model in analysis["models"].as_array().context("missing models")? {
            text.push_str(&format!(
                "| {} | {} | {:.9} | {:.9} |\n",
                analysis["family"].as_str().unwrap(),
                model["model"].as_str().unwrap(),
                model["macro_assay_position_mse"].as_f64().unwrap(),
                model["macro_assay_position_mae"].as_f64().unwrap()
            ));
        }
    }
    text.push_str("\nPrimary loss first averages alleles within position, positions within assay, then assays equally. HBB was already explored. Row-random results are a deliberately easier/leaky diagnostic, not primary validation. Input quality keeps every passing allele, including incomplete positions. No score p-value or count is a biological predictor.\n\n## Same-sequence context diagnostics\n\n| Relationship | Matched variants | Coefficient MSE | Conditional context-blind label MSE bound |\n| --- | ---: | ---: | ---: |\n");
    for pair in context["pairs"].as_array().context("missing pairs")? {
        text.push_str(&format!(
            "| {} | {} | {:.9} | {:.9} |\n",
            pair["relationship"].as_str().unwrap(),
            pair["matched_variants"],
            pair["mse"].as_f64().unwrap(),
            pair["finite_context_blind_mse_lower_bound"]
                .as_f64()
                .unwrap()
        ));
    }
    text.push_str("\nThese are finite source-estimate bounds, not biological irreducibility. Noise, library, time, cell and perturbation context can all contribute. A repeat in another element is not a gene-specific noise baseline. Overlapping TERT comparisons are not independent samples.\n\n## Source forensics and validation\n\nThe source-anchored HBB audit maps 62 database equality descriptions to original deletion rows. They produce 53 distinct altered target sequences. Raw duplicate-label inconsistency disappears after the predefined tag filter: 35 qualified rows have 35 distinct sequences. Related representation issues have prior public reports; global novelty is unestablished.\n\nFull predictions, fits, fold/guard identities and source-derived features are in study.json.gz. Bootstrap draws and spatial controls are in resampling.json. Their tails are conditional robustness diagnostics, not calibrated biological-significance probabilities. Ninety-nine controls and six-comparison Bonferroni adjustment cannot resolve a tail below 0.05.\n\nRun the independent Python implementation and comparator before claiming numerical reproduction. Preserve failed attempts; exact-tie rank metrics can be unstable under floating reduction order even when continuous predictions/losses agree. Read the source literature/novelty ledger before labeling an observation new.\n");
    Ok(text)
}

pub fn audit(root: &Path, artifacts: &super::Artifacts, lean: &Path) -> Result<()> {
    ensure!(
        artifacts
            .files
            .keys()
            .map(String::as_str)
            .collect::<Vec<_>>()
            == FILES,
        "incomplete/unknown study artifact inventory"
    );
    for (name, expected) in &artifacts.files {
        let path = root.join(name).canonicalize()?;
        ensure!(path.starts_with(root), "study artifact escapes directory");
        ensure!(
            genome_data::sha256(&fs::read(path)?) == *expected,
            "study checksum mismatch: {name}"
        );
    }
    let read =
        |name: &str| -> Result<Value> { Ok(serde_json::from_slice(&fs::read(root.join(name))?)?) };
    let inputs = read("inputs.json")?;
    let summary = read("summary.json")?;
    let protocol_bytes = fs::read(root.join("protocol.json"))?;
    let protocol = Protocol::from_bytes(&protocol_bytes)?;
    ensure!(
        summary["protocol_sha256"] == protocol.sha256()
            && inputs["protocol_sha256"] == protocol.sha256(),
        "study protocol commitments differ"
    );
    let revisions: Vec<Revision> = serde_json::from_slice(&fs::read(root.join("ledger.json"))?)?;
    Ledger::verify(&revisions)?;
    for (event, payload) in [
        ("freeze_protocol_and_inputs", inputs.clone()),
        ("validate_metadata_selected_cohort", read("cohort.json")?),
        ("validate_prior_motif_assets", read("motifs.json")?),
        ("execute_frozen_evaluation", summary.clone()),
        ("compare_source_contexts", read("context.json")?),
        ("retain_source_forensics", read("forensics.json")?),
    ] {
        let matches: Vec<_> = revisions.iter().filter(|r| r.kind == event).collect();
        ensure!(
            matches.len() == 1 && matches[0].payload == payload,
            "study payload differs from unique ledger event: {event}"
        );
    }
    let resampling = fs::read(root.join("resampling.json"))?;
    let resampling_event = revisions
        .iter()
        .find(|r| r.kind == "execute_fixed_resampling")
        .context("missing resampling commitment")?;
    ensure!(
        resampling_event.payload["sha256"] == genome_data::sha256(&resampling),
        "resampling commitment differs"
    );
    let file = fs::File::open(root.join("study.json.gz"))?;
    let mut bytes = Vec::new();
    GzDecoder::new(file)
        .take(512 * 1024 * 1024 + 1)
        .read_to_end(&mut bytes)?;
    ensure!(
        bytes.len() <= 512 * 1024 * 1024,
        "decompressed study exceeds budget"
    );
    let full: Value = serde_json::from_slice(&bytes)?;
    for key in [
        "protocol_sha256",
        "input_records_sha256",
        "motif_manifest_sha256",
        "inclusion",
        "feature_schema",
        "metric_arithmetic",
        "limitations",
    ] {
        ensure!(
            full[key] == summary[key],
            "compressed study differs from summary: {key}"
        );
    }
    for analysis in summary["analyses"]
        .as_array()
        .context("missing summary families")?
    {
        let full_analysis = full["analyses"]
            .as_array()
            .context("missing full families")?
            .iter()
            .find(|a| a["family"] == analysis["family"])
            .context("missing full evaluation family")?;
        ensure!(
            full_analysis["folds"] == analysis["fold_diagnostics"],
            "fold diagnostics differ across artifacts"
        );
        for model in analysis["models"]
            .as_array()
            .context("missing summary models")?
        {
            let full_model = full_analysis["models"]
                .as_array()
                .context("missing full models")?
                .iter()
                .find(|m| m["model"] == model["model"])
                .context("missing fitted model")?;
            for key in [
                "macro_assay_position_mse",
                "macro_assay_position_mae",
                "pooled_row_metrics",
                "per_assay_metrics",
            ] {
                ensure!(
                    full_model[key] == model[key],
                    "model metrics differ across artifacts: {key}"
                );
            }
            ensure!(
                full_model["predictions"]
                    .as_array()
                    .context("missing full predictions")?
                    .len() as u64
                    == model["prediction_count"]
                        .as_u64()
                        .context("missing prediction count")?,
                "compressed prediction count differs from summary"
            );
        }
    }
    let controls: Value = serde_json::from_slice(&resampling)?;
    ensure!(
        full["resampling"] == controls["bootstrap"]
            && full["spatial_null"] == controls["spatial_null"],
        "resampling artifacts differ"
    );
    let proof =
        genome_proof::check_regulatory(lean, &fs::read_to_string(root.join("proof.lean"))?)?;
    let formal = revisions
        .iter()
        .find(|r| r.kind == "check_scoped_method_contracts")
        .context("missing scoped proof event")?;
    ensure!(
        formal.payload == serde_json::to_value(proof)?,
        "scoped method proof differs from ledger"
    );
    println!("Study hashes, compressed results, summary/ledger bindings, protocol and scoped proof verified.");
    println!(
        "This audit does not substitute for independent numerical or experimental replication."
    );
    Ok(())
}
