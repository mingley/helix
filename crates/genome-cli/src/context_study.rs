//! Published specificity is a distinct endpoint from raw reporter activity.
use anyhow::{ensure, Context, Result};
use flate2::{read::GzDecoder, Compression, GzBuilder};
use genome_context_model::{
    evaluate, paired_group_bootstrap, ContextEvaluation, GroupSplit, LossWeighting, MatchedPanel,
    MeasuredOutcome, ModelComparison, ModelKind, OutcomeMatrix, SequenceRecord,
};
use genome_knowledge::{Ledger, Revision};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::{Read, Write},
    path::{Path, PathBuf},
    process::Command,
};

const PROTOCOL_SHA: &str = "8f0576845fecb5cfcf9bfa354fab446d849244448a3489905a6ce33bdeb21cec";
const ENDPOINT: &str = "published_element_specificity_score";
const FILES: [&str; 11] = [
    "features.json.gz",
    "findings.md",
    "grouping.json",
    "inputs.json",
    "ledger.json",
    "protocol.json",
    "reference.json",
    "resampling.json.gz",
    "source-manifest.json",
    "study.json.gz",
    "summary.json",
];

#[derive(clap::Args)]
pub struct ContextArgs {
    #[arg(long, default_value = "experiments/context-specificity-v2.json")]
    protocol: PathBuf,
    #[arg(long, default_value = "datasets/context-panel")]
    dataset: PathBuf,
    #[arg(long, default_value = "datasets/context-panel/groups.json")]
    groups: PathBuf,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Deserialize)]
struct SourcePanel {
    endpoint: String,
    contexts: Vec<String>,
    records: Vec<SourceRecord>,
}
#[derive(Deserialize)]
struct SourceRecord {
    id: String,
    sequence: String,
    observations: BTreeMap<String, f64>,
}

fn load(path: &Path) -> Result<Value> {
    Ok(serde_json::from_slice(&fs::read(path)?)?)
}
fn compressed<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let file = fs::OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(path)?;
    let mut out = GzBuilder::new()
        .mtime(0)
        .write(file, Compression::default());
    serde_json::to_writer(&mut out, value)?;
    out.write_all(b"\n")?;
    out.finish()?.sync_all()?;
    Ok(())
}
fn decompress(path: &Path) -> Result<Value> {
    let mut bytes = Vec::new();
    GzDecoder::new(fs::File::open(path)?)
        .take(1024 * 1024 * 1024 + 1)
        .read_to_end(&mut bytes)?;
    ensure!(
        bytes.len() <= 1024 * 1024 * 1024,
        "context artifact exceeds decompression budget"
    );
    Ok(serde_json::from_slice(&bytes)?)
}
fn implementation_hashes() -> Result<BTreeMap<String, String>> {
    fn walk(root: &Path, hashes: &mut BTreeMap<String, String>) -> Result<()> {
        for entry in fs::read_dir(root)? {
            let entry = entry?;
            let path = entry.path();
            if entry.file_type()?.is_dir() {
                walk(&path, hashes)?;
            } else if path.extension().is_some_and(|v| v == "rs")
                || path.file_name().is_some_and(|v| v == "Cargo.toml")
            {
                hashes.insert(
                    path.to_string_lossy().into_owned(),
                    genome_data::sha256(&fs::read(path)?),
                );
            }
        }
        Ok(())
    }
    let mut hashes = BTreeMap::new();
    walk(Path::new("crates"), &mut hashes)?;
    for name in [
        "Cargo.toml",
        "Cargo.lock",
        "scripts/group_context_panel.py",
        "scripts/fetch_context_panel.py",
        "scripts/twobit_reference.py",
    ] {
        hashes.insert(name.into(), genome_data::sha256(&fs::read(name)?));
    }
    Ok(hashes)
}
fn summarize(evaluation: &ContextEvaluation) -> Value {
    json!({"schema_version":1,"endpoint":ENDPOINT,"panel_metadata_sha256":evaluation.panel_metadata_sha256,
        "outcomes_sha256":evaluation.outcomes_sha256,"sequence_count":evaluation.sequence_count,
        "group_count":evaluation.group_count,"measurement_count":evaluation.measurement_count,
        "alpha":evaluation.alpha,"contexts":evaluation.contexts,"model_order":evaluation.model_order,
        "models":evaluation.models,"balanced_design_checks":evaluation.balanced_design_checks,
        "folds":evaluation.folds,
        "limitations":evaluation.interpretation_limits,
        "target_processing":"Source papers describe whole-panel cell z-scores and across-cell centering, including future held-out elements. The frozen workbook violates unweighted centering and its historical formula remains unresolved. No prospective or training-only target preprocessing is claimed.",
        "formal_scope":"No proof of fitted coefficients or biological mechanism is claimed; numerical invariants and independent computation are distinct checks."})
}

pub fn run(args: &ContextArgs) -> Result<()> {
    ensure!(
        !args.out.exists(),
        "output exists; preserve previous attempts"
    );
    let protocol_bytes = fs::read(&args.protocol)?;
    ensure!(
        genome_data::sha256(&protocol_bytes) == PROTOCOL_SHA,
        "frozen context-specificity v2 protocol differs"
    );
    let root = args.dataset.canonicalize()?;
    let manifest_bytes = fs::read(root.join("manifest.json"))?;
    let manifest: Value = serde_json::from_slice(&manifest_bytes)?;
    ensure!(
        manifest["protocol_sha256"] == PROTOCOL_SHA && manifest["endpoint"] == ENDPOINT,
        "source endpoint/protocol mismatch"
    );
    for (name, expected) in manifest["files"]
        .as_object()
        .context("missing source inventory")?
    {
        let path = root.join(name).canonicalize()?;
        ensure!(
            path.starts_with(&root) && expected == &genome_data::sha256(&fs::read(path)?),
            "source artifact mismatch: {name}"
        );
    }
    let specification: Value = serde_json::from_slice(&protocol_bytes)?;
    for (key, name) in [
        ("design_workbook", "raw/Table-S7-joint-design.xlsx"),
        ("folds_workbook", "raw/Table-S9-folds-performance.xlsx"),
    ] {
        ensure!(
            specification["source"][key]["sha256"]
                == genome_data::sha256(&fs::read(root.join(name))?),
            "raw workbook differs from frozen source: {name}"
        );
    }
    let panel_bytes = fs::read(root.join("panel.json"))?;
    ensure!(
        manifest["files"]["panel.json"] == genome_data::sha256(&panel_bytes),
        "panel missing from source commitments"
    );
    let source: SourcePanel = serde_json::from_slice(&panel_bytes)?;
    ensure!(
        source.contexts == ["HepG2", "K562", "WTC11"],
        "frozen contexts differ"
    );
    ensure!(
        source.endpoint == ENDPOINT,
        "raw activity and published specificity cannot be substituted"
    );
    let reference = load(&root.join("reference-validation.json"))?;
    ensure!(
        manifest["files"]["reference-validation.json"]
            == genome_data::sha256(&fs::read(root.join("reference-validation.json"))?),
        "reference validation missing from source commitments"
    );
    ensure!(
        reference["reference_verified"] == true
            && reference["records_checked"].as_u64() == Some(source.records.len() as u64),
        "incomplete genomic reference validation"
    );
    let grouping = load(&args.groups)?;
    ensure!(
        grouping["protocol_sha256"] == PROTOCOL_SHA
            && grouping["panel_sha256"] == genome_data::sha256(&panel_bytes),
        "grouping not bound to panel/protocol"
    );
    let mut groups = BTreeMap::new();
    for row in grouping["records"]
        .as_array()
        .context("missing grouping records")?
    {
        let id = row["id"].as_str().context("missing grouped ID")?;
        ensure!(
            groups
                .insert(
                    id.to_owned(),
                    row["group"].as_str().context("missing group")?.to_owned()
                )
                .is_none(),
            "duplicate grouping ID"
        );
    }
    ensure!(
        groups.len() == source.records.len(),
        "grouping coverage differs"
    );
    let mut group_sizes = BTreeMap::<String, usize>::new();
    for group in groups.values() {
        *group_sizes.entry(group.clone()).or_default() += 1;
    }
    ensure!(
        group_sizes.len() >= 100
            && group_sizes
                .values()
                .all(|n| n * 5 <= source.records.len() * 2),
        "component feasibility differs from frozen protocol"
    );
    let records = source
        .records
        .iter()
        .map(|r| {
            Ok(SequenceRecord {
                id: r.id.clone(),
                sequence: r.sequence.clone(),
                group: groups
                    .get(&r.id)
                    .context("ungrouped source sequence")?
                    .clone(),
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let panel = MatchedPanel::new(records, source.contexts)?;
    let measurements = source
        .records
        .iter()
        .flat_map(|r| {
            r.observations.iter().map(move |(c, y)| MeasuredOutcome {
                sequence_id: r.id.clone(),
                context: c.clone(),
                value: *y,
            })
        })
        .collect::<Vec<_>>();
    let outcomes = OutcomeMatrix::from_records(&panel, &measurements)?;
    let folds = grouping["folds"]
        .as_array()
        .context("missing folds")?
        .iter()
        .map(|f| {
            let strings = |key: &str| -> Result<Vec<String>> {
                f[key]
                    .as_array()
                    .context("missing fold inventory")?
                    .iter()
                    .map(|v| Ok(v.as_str().context("invalid fold ID")?.to_owned()))
                    .collect()
            };
            let split =
                GroupSplit::with_guard(&panel, &strings("test_groups")?, &strings("guard_ids")?)?;
            let guard_groups: BTreeSet<_> = strings("guard_groups")?.into_iter().collect();
            let expected_guards: Vec<_> = panel
                .sequences()
                .iter()
                .filter(|record| guard_groups.contains(&record.group))
                .map(|record| record.id.clone())
                .collect();
            ensure!(
                expected_guards == strings("guard_ids")?,
                "guard must exclude every member of the entire component"
            );
            ensure!(
                split.training_sequence_ids == strings("train_ids")?
                    && split.test_sequence_ids == strings("test_ids")?
                    && split.excluded_guard_sequence_ids == strings("guard_ids")?,
                "engine/grouping fold memberships differ"
            );
            ensure!(
                split.training_sequence_ids.len() * 2 >= panel.sequences().len(),
                "frozen minimum training fraction failed"
            );
            Ok(split)
        })
        .collect::<Result<Vec<_>>>()?;
    ensure!(folds.len() == 5, "frozen study requires five folds");
    let inputs = json!({"protocol_sha256":PROTOCOL_SHA,"endpoint":ENDPOINT,
        "source_manifest_sha256":genome_data::sha256(&manifest_bytes),
        "panel_sha256":genome_data::sha256(&panel_bytes),"grouping_sha256":genome_data::sha256(&fs::read(&args.groups)?),
        "implementation_files":implementation_hashes()?,
        "software":{"commit":Command::new("git").args(["rev-parse","HEAD"]).output().ok().filter(|o|o.status.success()).map(|o|String::from_utf8_lossy(&o.stdout).trim().to_owned()),
        "working_tree_dirty":Command::new("git").args(["status","--porcelain"]).output().ok().map(|o|!o.status.success()||!o.stdout.is_empty())}});
    let mut ledger = Ledger::default();
    ledger.append("freeze_specificity_protocol_and_inputs", &inputs)?;
    ledger.append("validate_published_source_endpoint", &manifest)?;
    ledger.append("validate_genomic_reference", &reference)?;
    ledger.append("freeze_sequence_only_grouping", &grouping)?;
    eprintln!(
        "Fitting five fixed models on {} elements, {} groups and three known contexts",
        panel.sequences().len(),
        panel.groups().len()
    );
    let evaluation = evaluate(&panel, &folds, &outcomes, 10.0)?;
    let checks = &evaluation.balanced_design_checks;
    ensure!(
        checks.max_abs_mean_independent_beta_minus_shared_beta <= 1e-7
            && checks.max_abs_mean_independent_prediction_minus_shared_prediction <= 1e-7
            && checks.max_abs_per_sequence_mse_contrast_identity_residual <= 1e-8,
        "balanced model invariants failed"
    );
    let comparisons = [
        ModelComparison {
            baseline: ModelKind::SharedDna,
            candidate: ModelKind::AdditiveContextDna,
        },
        ModelComparison {
            baseline: ModelKind::ContextMean,
            candidate: ModelKind::AdditiveContextDna,
        },
        ModelComparison {
            baseline: ModelKind::AdditiveContextDna,
            candidate: ModelKind::IndependentContextDna,
        },
    ];
    let resampling = paired_group_bootstrap(
        &evaluation.per_group_errors,
        &comparisons,
        2000,
        0x2026102702,
        [0.025, 0.975],
        LossWeighting::EqualSequences,
    )?;
    let mut summary = summarize(&evaluation);
    let zero_mse = measurements
        .iter()
        .map(|m| m.value * m.value / measurements.len() as f64)
        .sum::<f64>();
    let zero_mae = measurements
        .iter()
        .map(|m| m.value.abs() / measurements.len() as f64)
        .sum::<f64>();
    let zero_sum = outcomes
        .values()
        .iter()
        .map(|r| r.iter().sum::<f64>().abs())
        .fold(0.0_f64, f64::max);
    summary["fixed_zero_prediction_diagnostic"] = json!({"mse":zero_mse,"mae":zero_mae,"maximum_source_element_sum_absolute":zero_sum,"not_a_selected_fitted_model":true});
    summary["source_semantics"] = json!({"status":if zero_sum>1e-8 {"unresolved_source_formula_mismatch"} else {"row_sum_check_passed_not_historical_provenance"},
        "full_scientific_semantic_validation":false,"raw_source_values_preserved":true,
        "source_formula_zero_sum_tolerance":1e-8,
        "interpretation":"Published supplementary columns are benchmarked unchanged. A nonzero across-cell sum contradicts the stated unweighted centering formula; no projection or endpoint repair is permitted. Numerical reproduction is a separate claim."});
    summary["planned_comparisons"] = serde_json::to_value(&resampling.comparisons)?;
    ledger.append("execute_fixed_specificity_evaluation", &summary)?;
    ledger.append("commit_fixed_component_bootstrap",json!({"sha256":genome_data::sha256(&super::pretty(&serde_json::to_value(&resampling)?)?),"encoding":"canonical serde_json Value pretty JSON","replicates":2000,"conditional_not_biological_significance":true}))?;
    Ledger::verify(ledger.revisions())?;
    fs::create_dir_all(args.out.parent().unwrap_or(Path::new(".")))?;
    fs::create_dir(&args.out)?;
    let features = json!({"schema_version":1,"panel_metadata_sha256":panel.metadata_sha256(),
        "feature_names":genome_context_model::feature_names(),
        "records":panel.sequences().iter().zip(panel.features()).map(|(record,features)|json!({
            "id":record.id,"length":features.length,
            "gc_count":record.sequence.as_bytes().iter().filter(|base|matches!(base,b'G'|b'C')).count(),
            "canonical_trimer_counts":features.canonical_trimer_counts,"vector":features.vector()
        })).collect::<Vec<_>>()});
    compressed(&args.out.join("features.json.gz"), &features)?;
    compressed(&args.out.join("study.json.gz"), &evaluation)?;
    compressed(&args.out.join("resampling.json.gz"), &resampling)?;
    let mut narrative=String::from("# Published cell-specificity prediction\n\nThe target is source ESS, not raw log2 reporter activity. Source normalization used the full selected panel, including future held-out elements. Model fits alone are grouped and training-only. This is not a prospective target-preprocessing test, an unseen-cell test or a new biological mechanism.\n\n| Model | Equal-element MSE | Equal-element MAE |\n| --- | ---: | ---: |\n");
    for model in &evaluation.models {
        narrative.push_str(&format!(
            "| {:?} | {:.9} | {:.9} |\n",
            model.model, model.sequence_weighted_mse, model.sequence_weighted_mae
        ));
    }
    narrative.push_str(&format!("\nThe maximum absolute sum of three published scores is {zero_sum:.12}. This source check must be reconciled with the stated centering formula before claiming full semantic validation; original scores are retained unchanged. Under ideal unweighted centering, shared activity is removed, but that formula must not be assumed for this workbook. Mean/shared baselines do not establish genomic inactivity. Primary comparison is independent-context slopes versus additive context/shared-DNA slopes. All source, grouping, model and bootstrap outputs are retained. Independent numerical comparison is required separately; an integrity audit does not establish empirical truth or reproduction.\n"));
    for (name, bytes) in [
        ("protocol.json", protocol_bytes),
        ("source-manifest.json", manifest_bytes),
        (
            "reference.json",
            fs::read(root.join("reference-validation.json"))?,
        ),
        ("grouping.json", fs::read(&args.groups)?),
        ("inputs.json", super::pretty(&inputs)?),
        ("summary.json", super::pretty(&summary)?),
        ("ledger.json", super::pretty(&ledger.revisions())?),
        ("findings.md", narrative.into_bytes()),
    ] {
        super::write_new(&args.out.join(name), &bytes)?;
    }
    let files = FILES
        .iter()
        .map(|name| {
            Ok((
                (*name).to_owned(),
                genome_data::sha256(&fs::read(args.out.join(name))?),
            ))
        })
        .collect::<Result<BTreeMap<_, _>>>()?;
    super::write_new(
        &args.out.join("artifacts.json"),
        &super::pretty(&super::Artifacts {
            schema_version: 4,
            files,
        })?,
    )?;
    for model in &evaluation.models {
        println!("{:?}: MSE {:.9}", model.model, model.sequence_weighted_mse);
    }
    println!(
        "Saved complete context-specificity results to {}",
        args.out.display()
    );
    Ok(())
}

pub fn audit(root: &Path, artifacts: &super::Artifacts) -> Result<()> {
    ensure!(
        artifacts
            .files
            .keys()
            .map(String::as_str)
            .collect::<Vec<_>>()
            == FILES,
        "incomplete context-study inventory"
    );
    for (name, expected) in &artifacts.files {
        let path = root.join(name).canonicalize()?;
        ensure!(
            path.starts_with(root) && genome_data::sha256(&fs::read(path)?) == *expected,
            "context artifact mismatch: {name}"
        );
    }
    ensure!(
        genome_data::sha256(&fs::read(root.join("protocol.json"))?) == PROTOCOL_SHA,
        "context protocol differs"
    );
    let revisions: Vec<Revision> = serde_json::from_slice(&fs::read(root.join("ledger.json"))?)?;
    Ledger::verify(&revisions)?;
    for (event, name) in [
        ("freeze_specificity_protocol_and_inputs", "inputs.json"),
        ("validate_published_source_endpoint", "source-manifest.json"),
        ("validate_genomic_reference", "reference.json"),
        ("freeze_sequence_only_grouping", "grouping.json"),
        ("execute_fixed_specificity_evaluation", "summary.json"),
    ] {
        let matching = revisions
            .iter()
            .filter(|r| r.kind == event)
            .collect::<Vec<_>>();
        ensure!(
            matching.len() == 1 && matching[0].payload == load(&root.join(name))?,
            "context ledger binding differs: {name}"
        );
    }
    let full = decompress(&root.join("study.json.gz"))?;
    let features = decompress(&root.join("features.json.gz"))?;
    let summary = load(&root.join("summary.json"))?;
    ensure!(summary["endpoint"] == ENDPOINT, "context endpoint differs");
    let inputs = load(&root.join("inputs.json"))?;
    ensure!(
        inputs["grouping_sha256"] == genome_data::sha256(&fs::read(root.join("grouping.json"))?)
            && inputs["source_manifest_sha256"]
                == genome_data::sha256(&fs::read(root.join("source-manifest.json"))?),
        "archived input byte commitments differ"
    );
    let source_manifest = load(&root.join("source-manifest.json"))?;
    ensure!(
        source_manifest["files"]["reference-validation.json"]
            == genome_data::sha256(&fs::read(root.join("reference.json"))?),
        "archived reference bytes differ from source commitment"
    );
    let expected_contexts = json!(["HepG2", "K562", "WTC11"]);
    ensure!(
        full["contexts"] == expected_contexts,
        "archived context names differ"
    );
    let mut observations = BTreeMap::<String, BTreeMap<String, f64>>::new();
    for prediction in full["predictions"]
        .as_array()
        .context("missing observed predictions")?
    {
        let id = prediction["sequence_id"]
            .as_str()
            .context("missing observed ID")?
            .to_owned();
        let context = prediction["context"]
            .as_str()
            .context("missing observed context")?
            .to_owned();
        let observed = prediction["observed"]
            .as_f64()
            .context("invalid observation")?;
        ensure!(
            observed.is_finite()
                && observations
                    .entry(id)
                    .or_default()
                    .insert(context, observed)
                    .is_none(),
            "duplicate/nonfinite observation"
        );
    }
    let mut maximum_sum = 0.0_f64;
    for values in observations.values() {
        ensure!(
            values.keys().map(String::as_str).collect::<Vec<_>>() == ["HepG2", "K562", "WTC11"],
            "incomplete observed contexts"
        );
        maximum_sum = maximum_sum.max(values.values().sum::<f64>().abs());
    }
    ensure!(
        (summary["fixed_zero_prediction_diagnostic"]["maximum_source_element_sum_absolute"]
            .as_f64()
            .context("missing row-sum audit")?
            - maximum_sum)
            .abs()
            <= 1e-8,
        "source row-sum audit differs"
    );
    ensure!(
        summary["source_semantics"]["full_scientific_semantic_validation"] == false
            && summary["source_semantics"]["raw_source_values_preserved"] == true,
        "source semantic uncertainty must remain explicit"
    );
    if maximum_sum > 1e-8 {
        ensure!(
            summary["source_semantics"]["status"] == "unresolved_source_formula_mismatch",
            "source formula mismatch was concealed"
        );
    }
    ensure!(
        features["panel_metadata_sha256"] == full["panel_metadata_sha256"]
            && features["records"]
                .as_array()
                .context("missing features")?
                .len() as u64
                == full["sequence_count"]
                    .as_u64()
                    .context("missing sequence count")?,
        "context feature coverage differs"
    );
    for key in [
        "panel_metadata_sha256",
        "outcomes_sha256",
        "sequence_count",
        "group_count",
        "measurement_count",
        "alpha",
        "contexts",
        "model_order",
        "models",
        "folds",
        "balanced_design_checks",
    ] {
        ensure!(
            summary[key] == full[key],
            "context summary/full mismatch: {key}"
        );
    }
    let resampling = decompress(&root.join("resampling.json.gz"))?;
    let record = revisions
        .iter()
        .find(|r| r.kind == "commit_fixed_component_bootstrap")
        .context("missing bootstrap commitment")?;
    ensure!(
        record.payload["sha256"] == genome_data::sha256(&super::pretty(&resampling)?),
        "context bootstrap commitment differs"
    );
    ensure!(
        summary["planned_comparisons"] == resampling["comparisons"],
        "context contrast summary differs"
    );
    println!("Context-study hashes, source/ledger bindings and complete result summaries verified. No independent numerical or biological validation is implied.");
    Ok(())
}
