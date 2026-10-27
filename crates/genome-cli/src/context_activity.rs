//! Frozen v3 activity endpoint. No specificity outcomes or fitted responses are reused.
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

const PROTOCOL_SHA: &str = "64e14037fb1f42ddc1ecaca1d23bb8096fabac1f1e892d105eaa20f525fbbb20";
const RECOVERY_SHA: &str = "00175e41cf9fa0cbf177ce8fc7ef143eb508aec8b11d058f7acb4458284dc868";
const DNA_SHA: &str = "fce6b2a1aaf5aa4d530a2e298dccb840bc7c594a05b6831fdb0c6c0f9e2a9be3";
const GROUPS_SHA: &str = "333069c4ab41ed8ea6976b0bc50cd08de9a33118ab60d683533b7ee32a6f517f";
const RAW_SHA: &str = "9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e";
const ENDPOINT: &str = "published_joint_log2_rna_dna_activity";
const CONTEXTS: [&str; 3] = ["HepG2", "K562", "WTC11"];
const SOURCE_FILES: [&str; 10] = [
    "activity.json.gz",
    "archive-validation.json",
    "extraction-receipt.json",
    "recovery-artifacts.json",
    "recovery-inputs.json",
    "recovery-summary.json",
    "source-binding.json",
    "source-export.json.gz",
    "source-export-receipt.json",
    "source-ledger.jsonl.gz",
];
const FILES: [&str; 17] = [
    "activity.json.gz",
    "archive-validation.json",
    "common-activity.json.gz",
    "extraction-receipt.json",
    "features.json.gz",
    "findings.md",
    "grouping.json",
    "inputs.json",
    "ledger.json",
    "protocol.json",
    "recovery-summary.json",
    "resampling.json.gz",
    "source-binding.json",
    "source-export.json.gz",
    "source-ledger.jsonl.gz",
    "source-manifest.json",
    "study.json.gz",
];
// Summary is a separate constant to keep the complete, sorted inventory explicit.
fn inventory_names(manifest: &Value) -> Result<Vec<String>> {
    let files = manifest["files"]
        .as_object()
        .context("missing source inventory")?;
    ensure!(
        SOURCE_FILES.iter().all(|name| files.contains_key(*name)),
        "source inventory lacks required v3 inputs"
    );
    let mut names: BTreeSet<String> = FILES.iter().map(|name| (*name).into()).collect();
    names.insert("summary.json".into());
    for name in files.keys() {
        ensure!(
            !name.is_empty()
                && !name.contains('/')
                && !name.contains('\\')
                && name != "."
                && name != ".."
                && (SOURCE_FILES.contains(&name.as_str()) || !names.contains(name))
                && name != "artifacts.json"
                && name != "failure.json",
            "unsafe/colliding source receipt name: {name}"
        );
        names.insert(name.clone());
    }
    Ok(names.into_iter().collect())
}

#[derive(clap::Args)]
pub struct ContextActivityArgs {
    #[arg(long, default_value = "experiments/context-activity-v3.json")]
    protocol: PathBuf,
    #[arg(long, default_value = "datasets/context-activity")]
    dataset: PathBuf,
    #[arg(
        long,
        default_value = "research/results/context-specificity-v2/core/grouping.json"
    )]
    groups: PathBuf,
    #[arg(long, default_value = "datasets/context-panel/panel.json")]
    dna_panel: PathBuf,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Deserialize)]
struct DnaPanel {
    records: Vec<DnaRecord>,
}
#[derive(Deserialize)]
struct DnaRecord {
    id: String,
    sequence: String,
}
#[derive(Deserialize)]
struct Activity {
    protocol_sha256: String,
    endpoint: String,
    contexts: Vec<String>,
    records: Vec<ActivityRecord>,
}
#[derive(Deserialize)]
struct ActivityRecord {
    id: String,
    observations: BTreeMap<String, f64>,
}

fn load(path: &Path) -> Result<Value> {
    Ok(serde_json::from_slice(&fs::read(path)?)?)
}
fn uncompressed(path: &Path) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    GzDecoder::new(fs::File::open(path)?)
        .take(1024 * 1024 * 1024 + 1)
        .read_to_end(&mut bytes)?;
    ensure!(
        bytes.len() <= 1024 * 1024 * 1024,
        "activity artifact exceeds decompression budget"
    );
    Ok(bytes)
}
fn load_gzip(path: &Path) -> Result<Value> {
    Ok(serde_json::from_slice(&uncompressed(path)?)?)
}
fn compressed(path: &Path, value: &impl Serialize) -> Result<()> {
    let file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)?;
    let mut writer = GzBuilder::new()
        .mtime(0)
        .write(file, Compression::default());
    serde_json::to_writer(&mut writer, value)?;
    writer.write_all(b"\n")?;
    writer.finish()?.sync_all()?;
    Ok(())
}
fn canonical_sha(value: &impl Serialize) -> Result<String> {
    Ok(genome_data::sha256(&super::pretty(&serde_json::to_value(
        value,
    )?)?))
}
fn checked_hash(path: &Path, expected: &str) -> Result<Vec<u8>> {
    let bytes = fs::read(path).with_context(|| format!("read {}", path.display()))?;
    ensure!(
        genome_data::sha256(&bytes) == expected,
        "checksum mismatch: {}",
        path.display()
    );
    Ok(bytes)
}
fn source_file_inventory(root: &Path, manifest: &Value) -> Result<()> {
    let files = manifest["files"]
        .as_object()
        .context("missing source inventory")?;
    inventory_names(manifest)?;
    let canonical = root.canonicalize()?;
    for (name, expected) in files {
        let path = root.join(name).canonicalize()?;
        ensure!(
            path.starts_with(&canonical),
            "source inventory path escapes dataset"
        );
        checked_hash(&path, expected.as_str().context("invalid source hash")?)?;
    }
    Ok(())
}
fn source_gates(manifest: &Value) -> Result<()> {
    ensure!(
        manifest["protocol_sha256"] == PROTOCOL_SHA && manifest["endpoint"] == ENDPOINT,
        "source endpoint/protocol mismatch"
    );
    ensure!(
        manifest["source"]["sha256"] == RAW_SHA && manifest["source"]["bytes"] == 53_557_805,
        "raw activity source differs"
    );
    for gate in [
        "archive_integrity_verified",
        "empirical_activity_schema_bound",
        "coverage_complete",
        "design_identity_verified",
        "replicate_aggregation_status_recorded",
    ] {
        ensure!(
            manifest["gates"][gate] == true,
            "activity source gate failed: {gate}"
        );
    }
    ensure!(
        manifest["schema_binding"]
            .as_object()
            .is_some_and(|value| !value.is_empty())
            && !manifest["aggregation_status"].is_null(),
        "schema or aggregation status is absent"
    );
    Ok(())
}
fn check_source_receipts(root: &Path, manifest: &Value) -> Result<()> {
    ensure!(
        load(&root.join("source-binding.json"))? == manifest["schema_binding"],
        "source binding differs from manifest"
    );
    ensure!(
        load(&root.join("recovery-summary.json"))? == manifest["aggregation_status"],
        "aggregation status differs from manifest"
    );
    let archive = load(&root.join("archive-validation.json"))?;
    ensure!(
        archive["complete"] == true
            && archive["all_members_crc_verified"] == true
            && archive["status"] == "archive_integrity_verified"
            && archive["protocol_sha256"] == RECOVERY_SHA,
        "archive integrity stage incomplete"
    );
    let extraction = load(&root.join("extraction-receipt.json"))?;
    ensure!(
        extraction["protocol_sha256"] == RECOVERY_SHA
            && extraction["archive_sha256"] == archive["archive_sha256"],
        "extraction/source archive binding differs"
    );
    let members = extraction["selected_members"]
        .as_array()
        .context("missing extracted member inventory")?;
    ensure!(
        members.iter().any(|member| member["sha256"] == RAW_SHA
            && member["uncompressed_bytes"] == 53_557_805
            && member["crc32_verified_by_zipfile"] == true),
        "official activity member was not integrity verified"
    );
    let exported = load(&root.join("source-export-receipt.json"))?;
    ensure!(
        exported["protocol_sha256"] == RECOVERY_SHA
            && exported["source"]["sha256"] == RAW_SHA
            && exported["schema_binding_sha256"] == manifest["files"]["source-binding.json"]
            && exported["export"]["sha256"] == manifest["files"]["source-export.json.gz"]
            && exported["ledger"]["sha256"] == manifest["files"]["source-ledger.jsonl.gz"],
        "source export/ledger receipt is stale"
    );
    let recovery = load(&root.join("recovery-inputs.json"))?;
    let recovery_artifacts = load(&root.join("recovery-artifacts.json"))?;
    ensure!(
        recovery["protocol_sha256"] == RECOVERY_SHA
            && recovery_artifacts["protocol_sha256"] == RECOVERY_SHA
            && recovery_artifacts["files"]["summary.json"]["sha256"]
                == manifest["files"]["recovery-summary.json"],
        "recovery summary input/result binding differs"
    );
    let input_hashes = recovery["files"]
        .as_object()
        .context("missing recovery input commitments")?;
    for expected in [
        json!(RAW_SHA),
        json!(DNA_SHA),
        manifest["files"]["source-export.json.gz"].clone(),
        manifest["files"]["source-binding.json"].clone(),
        manifest["files"]["source-ledger.jsonl.gz"].clone(),
        manifest["files"]["source-export-receipt.json"].clone(),
        manifest["files"]["archive-validation.json"].clone(),
    ] {
        ensure!(
            input_hashes
                .values()
                .any(|record| record["sha256"] == expected),
            "recovery arithmetic consumed a different source input"
        );
    }
    Ok(())
}
fn validate_activity(
    activity: &Activity,
    expected_ids: &BTreeSet<String>,
) -> Result<Vec<MeasuredOutcome>> {
    ensure!(
        activity.protocol_sha256 == PROTOCOL_SHA
            && activity.endpoint == ENDPOINT
            && activity
                .contexts
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>()
                == CONTEXTS,
        "activity endpoint/context identity differs"
    );
    ensure!(
        activity.records.len() == expected_ids.len(),
        "activity coverage count differs"
    );
    let mut seen = BTreeSet::new();
    let mut measurements = Vec::with_capacity(expected_ids.len() * 3);
    let mut previous: Option<&str> = None;
    for record in &activity.records {
        ensure!(
            previous.is_none_or(|id| id < record.id.as_str()),
            "activity IDs must be unique and canonically sorted"
        );
        previous = Some(&record.id);
        ensure!(
            expected_ids.contains(&record.id) && seen.insert(record.id.clone()),
            "extra/duplicate activity ID"
        );
        ensure!(
            record
                .observations
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>()
                == CONTEXTS,
            "incomplete/extra activity context"
        );
        for (context, value) in &record.observations {
            ensure!(value.is_finite(), "nonfinite published activity");
            measurements.push(MeasuredOutcome {
                sequence_id: record.id.clone(),
                context: context.clone(),
                value: *value,
            });
        }
    }
    ensure!(seen == *expected_ids, "activity ID coverage differs");
    Ok(measurements)
}
fn verify_export_values(export: &Value, activity: &Activity) -> Result<()> {
    ensure!(
        export["protocol_sha256"] == RECOVERY_SHA && export["source"]["sha256"] == RAW_SHA,
        "activity export source/protocol binding differs"
    );
    let table = &export["tables"]["all_cell_types_summary"];
    ensure!(
        table["id_header"] == "name",
        "source literal ID header differs"
    );
    let headers = table["headers"]
        .as_object()
        .context("missing source headers")?;
    for context in CONTEXTS {
        let header = format!("{context} [log2(rna/dna)]");
        ensure!(
            headers.values().filter(|value| **value == header).count() == 1,
            "source header absent/ambiguous: {header}"
        );
    }
    let required: BTreeMap<_, _> = activity
        .records
        .iter()
        .map(|record| (record.id.as_str(), record))
        .collect();
    let mut seen = BTreeSet::new();
    let mut covered = BTreeSet::new();
    for row in table["rows"]
        .as_array()
        .context("missing source export rows")?
    {
        if row["id"].is_null() || row["id"] == "" {
            continue;
        }
        let id = row["id"].as_str().context("nonliteral source ID")?;
        ensure!(seen.insert(id), "duplicate source activity ID");
        if let Some(record) = required.get(id) {
            covered.insert(id);
            for context in CONTEXTS {
                let header = format!("{context} [log2(rna/dna)]");
                let raw = row["values"][&header]
                    .as_str()
                    .context("required activity is not a raw numeric literal")?;
                let value: f64 = raw.parse().context("invalid source numeric literal")?;
                ensure!(
                    value.is_finite() && value.to_bits() == record.observations[context].to_bits(),
                    "published source value changed: {id}/{context}"
                );
            }
        }
    }
    ensure!(
        covered.len() == required.len(),
        "source export does not cover every required activity ID"
    );
    Ok(())
}
fn verify_export_and_types(root: &Path, manifest: &Value, activity: &Activity) -> Result<()> {
    let export = load_gzip(&root.join("source-export.json.gz"))?;
    ensure!(
        export["schema_binding_sha256"] == manifest["files"]["source-binding.json"],
        "export/schema binding differs"
    );
    verify_export_values(&export, activity)?;
    drop(export);
    let ids: BTreeSet<_> = activity
        .records
        .iter()
        .map(|record| record.id.as_str())
        .collect();
    for line in uncompressed(&root.join("source-ledger.jsonl.gz"))?
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.is_empty())
    {
        let record: Value = serde_json::from_slice(line)?;
        if record["record_type"] == "cell_exception"
            && record["sheet"] == "all_cell_types_summary"
            && record["id"].as_str().is_some_and(|id| ids.contains(id))
        {
            let required = CONTEXTS
                .iter()
                .any(|context| record["header"] == format!("{context} [log2(rna/dna)]"));
            ensure!(
                !required,
                "required activity is not an ordinary finite numeric source cell"
            );
        }
    }
    Ok(())
}
fn panel_and_folds(
    dna_bytes: &[u8],
    group_bytes: &[u8],
    activity: &Activity,
) -> Result<(MatchedPanel, Vec<GroupSplit>, OutcomeMatrix)> {
    let dna: DnaPanel = serde_json::from_slice(dna_bytes)?;
    let grouping: Value = serde_json::from_slice(group_bytes)?;
    let mut groups = BTreeMap::new();
    for record in grouping["records"].as_array().context("missing groups")? {
        let id = record["id"].as_str().context("missing grouped ID")?;
        ensure!(
            groups
                .insert(
                    id.to_owned(),
                    record["group"]
                        .as_str()
                        .context("missing group ID")?
                        .to_owned()
                )
                .is_none(),
            "duplicate grouped ID"
        );
    }
    ensure!(
        dna.records.len() == 55_276 && groups.len() == 55_276,
        "fixed DNA/group cohort differs"
    );
    let records = dna
        .records
        .into_iter()
        .map(|record| {
            Ok(SequenceRecord {
                id: record.id.clone(),
                sequence: record.sequence,
                group: groups
                    .get(&record.id)
                    .context("DNA ID absent from fixed groups")?
                    .clone(),
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let panel = MatchedPanel::new(
        records,
        CONTEXTS.iter().map(|cell| (*cell).into()).collect(),
    )?;
    ensure!(
        panel.groups().len() == 2761,
        "frozen component count differs"
    );
    let expected = panel
        .sequences()
        .iter()
        .map(|record| record.id.clone())
        .collect();
    let outcomes = OutcomeMatrix::from_records(&panel, &validate_activity(activity, &expected)?)?;
    let mut folds = Vec::new();
    for fold in grouping["folds"]
        .as_array()
        .context("missing frozen folds")?
    {
        let strings = |name: &str| -> Result<Vec<String>> {
            fold[name]
                .as_array()
                .context("missing fold ID inventory")?
                .iter()
                .map(|value| Ok(value.as_str().context("invalid fold ID")?.to_owned()))
                .collect()
        };
        let split =
            GroupSplit::with_guard(&panel, &strings("test_groups")?, &strings("guard_ids")?)?;
        ensure!(
            split.training_sequence_ids == strings("train_ids")?
                && split.test_sequence_ids == strings("test_ids")?
                && split.excluded_guard_sequence_ids == strings("guard_ids")?,
            "frozen fold membership differs"
        );
        ensure!(
            split.training_sequence_ids.len() * 2 >= panel.sequences().len(),
            "training support gate failed"
        );
        folds.push(split);
    }
    ensure!(folds.len() == 5, "frozen fivefold partition required");
    Ok((panel, folds, outcomes))
}

#[derive(Clone, Deserialize)]
struct NamedPredictions {
    global_mean: f64,
    context_mean: f64,
    shared_dna: f64,
    additive_context_dna: f64,
    independent_context_dna: f64,
}
impl NamedPredictions {
    fn values(&self) -> [f64; 5] {
        [
            self.global_mean,
            self.context_mean,
            self.shared_dna,
            self.additive_context_dna,
            self.independent_context_dna,
        ]
    }
}
#[derive(Clone, Deserialize)]
struct SavedPrediction {
    fold_index_0based: usize,
    sequence_id: String,
    group: String,
    context: String,
    observed: f64,
    predictions: NamedPredictions,
}
fn native_predictions(evaluation: &ContextEvaluation) -> Vec<SavedPrediction> {
    evaluation
        .predictions
        .iter()
        .map(|p| SavedPrediction {
            fold_index_0based: p.fold_index_0based,
            sequence_id: p.sequence_id.clone(),
            group: p.group.clone(),
            context: p.context.clone(),
            observed: p.observed,
            predictions: NamedPredictions {
                global_mean: p.predictions.global_mean,
                context_mean: p.predictions.context_mean,
                shared_dna: p.predictions.shared_dna,
                additive_context_dna: p.predictions.additive_context_dna,
                independent_context_dna: p.predictions.independent_context_dna,
            },
        })
        .collect()
}
fn feature_dump(panel: &MatchedPanel) -> Value {
    json!({"schema_version":1,"panel_metadata_sha256":panel.metadata_sha256(),"feature_names":genome_context_model::feature_names(),
        "records":panel.sequences().iter().zip(panel.features()).map(|(record,features)|json!({"id":record.id,"length":features.length,
            "gc_count":record.sequence.bytes().filter(|base|matches!(base,b'C'|b'G')).count(),"canonical_trimer_counts":features.canonical_trimer_counts,
            "vector":features.vector()})).collect::<Vec<_>>()})
}
fn mean(values: &[f64]) -> Result<f64> {
    ensure!(
        !values.is_empty() && values.iter().all(|value| value.is_finite()),
        "mean requires finite values"
    );
    let result = if values.iter().all(|value| *value == values[0]) {
        values[0]
    } else {
        values.iter().map(|value| value / values.len() as f64).sum()
    };
    ensure!(result.is_finite(), "mean overflow");
    Ok(result)
}
fn common_activity(predictions: &[SavedPrediction]) -> Result<Value> {
    let mut elements: BTreeMap<&str, BTreeMap<&str, &SavedPrediction>> = BTreeMap::new();
    for prediction in predictions {
        ensure!(
            prediction.observed.is_finite()
                && prediction
                    .predictions
                    .values()
                    .iter()
                    .all(|value| value.is_finite()),
            "nonfinite saved prediction"
        );
        ensure!(
            elements
                .entry(&prediction.sequence_id)
                .or_default()
                .insert(&prediction.context, prediction)
                .is_none(),
            "duplicate prediction context"
        );
    }
    let mut records = Vec::new();
    let mut maximum = 0.0_f64;
    let mut global_errors = Vec::new();
    let mut shared_errors = Vec::new();
    let mut primary = Vec::new();
    let mut common = Vec::new();
    for (id, rows) in elements {
        ensure!(
            rows.keys().copied().collect::<Vec<_>>() == CONTEXTS,
            "incomplete prediction contexts"
        );
        let ordered: Vec<_> = CONTEXTS.iter().map(|context| rows[context]).collect();
        let first = ordered[0];
        ensure!(
            ordered
                .iter()
                .all(|row| row.fold_index_0based == first.fold_index_0based
                    && row.group == first.group
                    && row.predictions.global_mean == first.predictions.global_mean
                    && row.predictions.shared_dna == first.predictions.shared_dna),
            "mean prediction/fold identity differs across cells"
        );
        let observed_mean = mean(&ordered.iter().map(|row| row.observed).collect::<Vec<_>>())?;
        let global = (observed_mean - first.predictions.global_mean).powi(2);
        let shared = (observed_mean - first.predictions.shared_dna).powi(2);
        let context_error = mean(
            &ordered
                .iter()
                .map(|row| (row.observed - row.predictions.context_mean).powi(2))
                .collect::<Vec<_>>(),
        )?;
        let additive_error = mean(
            &ordered
                .iter()
                .map(|row| (row.observed - row.predictions.additive_context_dna).powi(2))
                .collect::<Vec<_>>(),
        )?;
        let gain = context_error - additive_error;
        let common_gain = global - shared;
        let residual = gain - common_gain;
        ensure!(
            residual.abs() <= 1e-8,
            "common-activity loss identity failed"
        );
        maximum = maximum.max(residual.abs());
        global_errors.push(global);
        shared_errors.push(shared);
        primary.push(gain);
        common.push(common_gain);
        records.push(json!({"sequence_id":id,"group":first.group,"fold_index_0based":first.fold_index_0based,"observed_mean":observed_mean,
            "global_mean_prediction":first.predictions.global_mean,"shared_dna_prediction":first.predictions.shared_dna,
            "global_mean_squared_error":global,"shared_dna_squared_error":shared,"C_minus_A_cell_MSE":gain,
            "G_minus_S_mean_activity_MSE":common_gain,"identity_residual":residual}));
    }
    Ok(
        json!({"scope":"Correlated algebraic diagnostic on the source's published cell scales; not absolute expression or independent validation",
        "records":records,"global_mean_mse":mean(&global_errors)?,"shared_dna_mse":mean(&shared_errors)?,
        "mean_activity_mse_improvement":mean(&common)?,"primary_cell_mse_improvement":mean(&primary)?,
        "maximum_absolute_per_element_identity_residual":maximum}),
    )
}
fn fold_metrics(predictions: &[SavedPrediction]) -> Result<Value> {
    let mut output = Vec::new();
    for fold in 0..5 {
        let rows: Vec<_> = predictions
            .iter()
            .filter(|row| row.fold_index_0based == fold)
            .collect();
        ensure!(!rows.is_empty(), "empty fold prediction set");
        let models=ModelKind::ALL.iter().enumerate().map(|(j,model)|->Result<Value>{
            let contexts=CONTEXTS.iter().map(|context|->Result<Value>{let subset:Vec<_>=rows.iter().filter(|row|row.context==*context).collect();
                Ok(json!({"context":context,"mse":mean(&subset.iter().map(|row|(row.observed-row.predictions.values()[j]).powi(2)).collect::<Vec<_>>())?,
                    "mae":mean(&subset.iter().map(|row|(row.observed-row.predictions.values()[j]).abs()).collect::<Vec<_>>())?}))}).collect::<Result<Vec<_>>>()?;
            Ok(json!({"model":model,"mse":mean(&rows.iter().map(|row|(row.observed-row.predictions.values()[j]).powi(2)).collect::<Vec<_>>())?,
                "mae":mean(&rows.iter().map(|row|(row.observed-row.predictions.values()[j]).abs()).collect::<Vec<_>>())?,"by_context":contexts}))
        }).collect::<Result<Vec<_>>>()?;
        output.push(json!({"fold_index_0based":fold,"measurement_count":rows.len(),"sequence_count":rows.len()/3,"models":models}));
    }
    Ok(json!(output))
}
fn check_balanced(checks: &Value) -> Result<()> {
    for (name, tolerance) in [
        ("max_abs_additive_offset_prediction_residual", 1e-7),
        ("max_abs_mean_independent_beta_minus_shared_beta", 1e-7),
        (
            "max_abs_mean_independent_prediction_minus_shared_prediction",
            1e-7,
        ),
        ("max_abs_per_sequence_mse_contrast_identity_residual", 1e-8),
    ] {
        let value = checks[name]
            .as_f64()
            .context("missing balanced-design check")?;
        ensure!(
            value.is_finite() && value >= 0.0 && value <= tolerance,
            "balanced-design identity failed: {name}"
        );
    }
    Ok(())
}
fn implementation_hashes() -> Result<BTreeMap<String, String>> {
    fn walk(path: &Path, output: &mut BTreeMap<String, String>) -> Result<()> {
        for entry in fs::read_dir(path)? {
            let entry = entry?;
            let path = entry.path();
            if entry.file_type()?.is_dir() {
                walk(&path, output)?;
            } else if path.extension().is_some_and(|extension| extension == "rs")
                || path.file_name().is_some_and(|name| name == "Cargo.toml")
            {
                output.insert(
                    path.to_string_lossy().into_owned(),
                    genome_data::sha256(&fs::read(path)?),
                );
            }
        }
        Ok(())
    }
    let mut output = BTreeMap::new();
    walk(Path::new("crates"), &mut output)?;
    for path in [
        "Cargo.toml",
        "Cargo.lock",
        "scripts/build_context_activity.py",
    ] {
        output.insert(path.into(), genome_data::sha256(&fs::read(path)?));
    }
    Ok(output)
}

pub fn run(args: &ContextActivityArgs) -> Result<()> {
    ensure!(
        !args.out.exists(),
        "output exists; preserve previous attempts"
    );
    if let Some(parent) = args
        .out
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }
    fs::create_dir(&args.out)?;
    let result = run_inner(args);
    if let Err(error) = &result {
        super::write_new(
            &args.out.join("failure.json"),
            &super::pretty(&json!({"protocol_sha256":PROTOCOL_SHA,
        "status":"failed_attempt_not_a_successful_result","error":format!("{error:#}"),"prior_attempts_preserved":true}))?,
        )?;
    }
    result
}
fn run_inner(args: &ContextActivityArgs) -> Result<()> {
    let protocol_bytes = checked_hash(&args.protocol, PROTOCOL_SHA)?;
    let manifest_bytes = fs::read(args.dataset.join("manifest.json"))?;
    let manifest: Value = serde_json::from_slice(&manifest_bytes)?;
    source_gates(&manifest)?;
    source_file_inventory(&args.dataset, &manifest)?;
    check_source_receipts(&args.dataset, &manifest)?;
    let raw_path = args.dataset.join(
        manifest["source"]["raw_path"]
            .as_str()
            .context("missing raw source path")?,
    );
    let raw = checked_hash(&raw_path, RAW_SHA)?;
    ensure!(raw.len() == 53_557_805, "raw source byte count differs");
    drop(raw);
    let dna_bytes = checked_hash(&args.dna_panel, DNA_SHA)?;
    let group_bytes = checked_hash(&args.groups, GROUPS_SHA)?;
    let activity: Activity =
        serde_json::from_slice(&uncompressed(&args.dataset.join("activity.json.gz"))?)?;
    let (panel, folds, outcomes) = panel_and_folds(&dna_bytes, &group_bytes, &activity)?;
    verify_export_and_types(&args.dataset, &manifest, &activity)?;
    let inputs = json!({"protocol_sha256":PROTOCOL_SHA,"endpoint":ENDPOINT,"source_manifest_sha256":genome_data::sha256(&manifest_bytes),
        "source_dataset_path":args.dataset,"raw_source_sha256":RAW_SHA,"raw_source_bytes":53_557_805,
        "dna_panel_path":args.dna_panel,"dna_panel_sha256":DNA_SHA,"grouping_sha256":GROUPS_SHA,"activity_sha256":manifest["files"]["activity.json.gz"],
        "panel_metadata_sha256":panel.metadata_sha256(),"outcomes_sha256":outcomes.outcomes_sha256(),"implementation_files":implementation_hashes()?,
        "software":{"commit":Command::new("git").args(["rev-parse","HEAD"]).output().ok().filter(|result|result.status.success()).map(|result|String::from_utf8_lossy(&result.stdout).trim().to_owned()),
            "working_tree_dirty":Command::new("git").args(["status","--porcelain"]).output().ok().map(|result|!result.status.success()||!result.stdout.is_empty())}});
    let mut ledger = Ledger::default();
    ledger.append("freeze_activity_protocol_and_inputs", &inputs)?;
    ledger.append("validate_activity_source_gates", &manifest)?;
    for (name, bytes) in [
        ("protocol.json", protocol_bytes),
        ("source-manifest.json", manifest_bytes),
        ("grouping.json", group_bytes),
        ("inputs.json", super::pretty(&inputs)?),
    ] {
        super::write_new(&args.out.join(name), &bytes)?;
    }
    for name in manifest["files"].as_object().unwrap().keys() {
        super::write_new(&args.out.join(name), &fs::read(args.dataset.join(name))?)?;
    }
    eprintln!(
        "Fitting fresh activity-v3 responses: {} elements, three known cells, five fixed models",
        panel.sequences().len()
    );
    let evaluation = evaluate(&panel, &folds, &outcomes, 10.0)?;
    let checks = serde_json::to_value(&evaluation.balanced_design_checks)?;
    check_balanced(&checks)?;
    let predictions = native_predictions(&evaluation);
    let common = common_activity(&predictions)?;
    let fold_summary = fold_metrics(&predictions)?;
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
    close(
        resampling.comparisons[1].observed_mse_improvement,
        common["mean_activity_mse_improvement"]
            .as_f64()
            .context("missing common mean gain")?,
        1e-8,
        "primary/common mean gain",
    )?;
    let mut common_summary = common.clone();
    common_summary.as_object_mut().unwrap().remove("records");
    let summary = json!({"schema_version":1,"endpoint":ENDPOINT,"protocol_sha256":PROTOCOL_SHA,"panel_metadata_sha256":evaluation.panel_metadata_sha256,
        "outcomes_sha256":evaluation.outcomes_sha256,"sequence_count":evaluation.sequence_count,"group_count":evaluation.group_count,"measurement_count":evaluation.measurement_count,
        "contexts":evaluation.contexts,"alpha":10.0,"model_order":evaluation.model_order,"models":evaluation.models,"folds":evaluation.folds,"fold_metrics":fold_summary,
        "balanced_design_checks":checks,"common_activity":common_summary,"planned_comparisons":resampling.comparisons,
        "primary_comparison_index_0based":1,"source_gate_passed":true,"aggregation_status":manifest["aggregation_status"],"source_schema_binding":manifest["schema_binding"],
        "target_processing":"Published library-normalized log2(RNA/DNA) aggregate values unchanged; source normalization and QC precede grouped splitting; no ESS outcome or response fit reused.",
        "limits":["Selected same-study reporter panel; not absolute transcripts per cell, unseen cells, independent biological replication or a causal mechanism.",
            "Common across-cell mean is measured on the published cell scales, not an absolute expression calibration.",
            "Conditional component bootstrap omits refitting/source-normalization uncertainty; no calibrated p-value or genome-wide prevalence inference."]});
    ledger.append("execute_fixed_activity_evaluation", &summary)?;
    ledger.append("commit_activity_numeric_artifacts",json!({"encoding":"canonical serde_json Value pretty JSON",
        "study_sha256":canonical_sha(&evaluation)?,"resampling_sha256":canonical_sha(&resampling)?,"common_activity_sha256":canonical_sha(&common)?}))?;
    Ledger::verify(ledger.revisions())?;
    compressed(&args.out.join("study.json.gz"), &evaluation)?;
    compressed(&args.out.join("features.json.gz"), &feature_dump(&panel))?;
    compressed(&args.out.join("resampling.json.gz"), &resampling)?;
    compressed(&args.out.join("common-activity.json.gz"), &common)?;
    super::write_new(&args.out.join("summary.json"), &super::pretty(&summary)?)?;
    super::write_new(
        &args.out.join("ledger.json"),
        &super::pretty(&ledger.revisions())?,
    )?;
    let mut findings=String::from("# Published reporter activity, frozen v3\n\nFresh activity responses on the fixed prior DNA/groups; published log2(RNA/DNA) aggregates are unchanged. Primary comparison is context mean versus additive shared-DNA/context model. This is the same selected experiment, with known cells; no independent biological replication or absolute transcription calibration.\n\n| Model | Equal-element MSE | Equal-element MAE |\n| --- | ---: | ---: |\n");
    for model in &evaluation.models {
        findings.push_str(&format!(
            "| {:?} | {:.9} | {:.9} |\n",
            model.model, model.sequence_weighted_mse, model.sequence_weighted_mae
        ));
    }
    findings.push_str("\nAll five models, three comparisons, shared component draws, fold losses, source-gate receipts and mean-activity identity are retained. Replicate aggregation and normalization history are reported separately; numerical reproduction does not settle them.\n");
    super::write_new(&args.out.join("findings.md"), findings.as_bytes())?;
    let files = inventory_names(&manifest)?
        .into_iter()
        .map(|name| {
            Ok((
                name.to_owned(),
                genome_data::sha256(&fs::read(args.out.join(name))?),
            ))
        })
        .collect::<Result<_>>()?;
    super::write_new(
        &args.out.join("artifacts.json"),
        &super::pretty(&super::Artifacts {
            schema_version: 5,
            files,
        })?,
    )?;
    println!(
        "Saved activity-v3 results to {}; primary C-minus-A MSE gain {:.9}",
        args.out.display(),
        resampling.comparisons[1].observed_mse_improvement
    );
    Ok(())
}

fn ledger_binding<'a>(
    revisions: &'a [Revision],
    kind: &str,
    payload: &Value,
) -> Result<&'a Revision> {
    let matches: Vec<_> = revisions
        .iter()
        .filter(|revision| revision.kind == kind)
        .collect();
    ensure!(
        matches.len() == 1 && matches[0].payload == *payload,
        "activity ledger binding differs: {kind}"
    );
    Ok(matches[0])
}
fn close(left: f64, right: f64, tolerance: f64, context: &str) -> Result<()> {
    ensure!(
        left.is_finite() && right.is_finite() && (left - right).abs() <= tolerance,
        "activity numeric binding differs: {context}"
    );
    Ok(())
}
fn check_prediction_inputs(
    predictions: &[SavedPrediction],
    panel: &MatchedPanel,
    folds: &[GroupSplit],
    outcomes: &OutcomeMatrix,
) -> Result<()> {
    ensure!(
        predictions.len() == panel.sequences().len() * 3,
        "prediction coverage differs"
    );
    let mut seen = BTreeSet::new();
    let fold_by_id: BTreeMap<_, _> = folds
        .iter()
        .enumerate()
        .flat_map(|(fold, split)| {
            split
                .test_sequence_ids
                .iter()
                .map(move |id| (id.as_str(), fold))
        })
        .collect();
    ensure!(
        fold_by_id.len() == panel.sequences().len(),
        "fold coverage differs"
    );
    for prediction in predictions {
        let i = panel.sequence_index(&prediction.sequence_id)?;
        let j = panel.context_index(&prediction.context)?;
        ensure!(
            seen.insert((i, j))
                && prediction.observed.to_bits() == outcomes.values()[i][j].to_bits()
                && prediction.group == panel.sequences()[i].group
                && fold_by_id[prediction.sequence_id.as_str()] == prediction.fold_index_0based,
            "saved observation/group/fold differs from frozen source"
        );
    }
    Ok(())
}

#[derive(Deserialize)]
struct SavedRidge {
    training_rows: usize,
    feature_count: usize,
    alpha: f64,
    intercept: f64,
    training_feature_means: Vec<f64>,
    training_feature_scales: Vec<f64>,
    standardized_coefficients: Vec<f64>,
}
impl SavedRidge {
    fn effect(&self, row: &[f64]) -> Result<f64> {
        ensure!(
            self.feature_count == 33
                && row.len() == 33
                && self.training_feature_means.len() == 33
                && self.training_feature_scales.len() == 33
                && self.standardized_coefficients.len() == 33
                && self.alpha == 10.0
                && self.intercept.is_finite()
                && self
                    .training_feature_means
                    .iter()
                    .chain(&self.standardized_coefficients)
                    .all(|value| value.is_finite())
                && self
                    .training_feature_scales
                    .iter()
                    .all(|value| value.is_finite() && *value > 0.0),
            "invalid saved ridge fit"
        );
        let value = row
            .iter()
            .enumerate()
            .map(|(j, x)| {
                (x - self.training_feature_means[j]) / self.training_feature_scales[j]
                    * self.standardized_coefficients[j]
            })
            .sum::<f64>();
        ensure!(value.is_finite(), "saved ridge prediction overflow");
        Ok(value)
    }
}
#[derive(Deserialize)]
struct SavedFit {
    training_sequence_ids: Vec<String>,
    contexts: Vec<String>,
    training_sequences: usize,
    training_measurements: usize,
    alpha: f64,
    feature_names: Vec<String>,
    global_mean: f64,
    training_context_means: BTreeMap<String, f64>,
    shared_dna_fit: SavedRidge,
    independent_context_fits: BTreeMap<String, SavedRidge>,
}

fn check_saved_fits(
    full: &Value,
    predictions: &[SavedPrediction],
    panel: &MatchedPanel,
    folds: &[GroupSplit],
    outcomes: &OutcomeMatrix,
) -> Result<()> {
    let saved = full["folds"].as_array().context("missing saved folds")?;
    ensure!(saved.len() == folds.len(), "saved fold count differs");
    let mut fits = Vec::new();
    for (k, fold) in saved.iter().enumerate() {
        ensure!(
            fold["fold_index_0based"] == k && fold["split"] == serde_json::to_value(&folds[k])?,
            "saved training/test/guard split differs"
        );
        let fit: SavedFit = serde_json::from_value(fold["fit"].clone())?;
        ensure!(
            fit.training_sequence_ids == folds[k].training_sequence_ids
                && fit.contexts == panel.contexts()
                && fit.training_sequences == folds[k].training_indices().len()
                && fit.training_measurements == fit.training_sequences * 3
                && fit.alpha == 10.0
                && fit.feature_names == genome_context_model::feature_names()
                && fit
                    .independent_context_fits
                    .keys()
                    .map(String::as_str)
                    .collect::<Vec<_>>()
                    == CONTEXTS
                && fit
                    .training_context_means
                    .keys()
                    .map(String::as_str)
                    .collect::<Vec<_>>()
                    == CONTEXTS,
            "saved fit shape/feature identity differs"
        );
        let training = folds[k].training_indices();
        let target_means = training
            .iter()
            .map(|i| mean(&outcomes.values()[*i]))
            .collect::<Result<Vec<_>>>()?;
        close(
            fit.global_mean,
            mean(&target_means)?,
            1e-7,
            "training global response mean",
        )?;
        close(
            fit.shared_dna_fit.intercept,
            fit.global_mean,
            1e-7,
            "shared response intercept",
        )?;
        for (j, context) in CONTEXTS.iter().enumerate() {
            let expected = mean(
                &training
                    .iter()
                    .map(|i| outcomes.values()[*i][j])
                    .collect::<Vec<_>>(),
            )?;
            close(
                fit.training_context_means[*context],
                expected,
                1e-7,
                "training context response mean",
            )?;
            let specific = &fit.independent_context_fits[*context];
            close(
                specific.intercept,
                expected,
                1e-7,
                "independent context intercept",
            )?;
            ensure!(
                specific.training_rows == training.len()
                    && specific.training_feature_means == fit.shared_dna_fit.training_feature_means
                    && specific.training_feature_scales
                        == fit.shared_dna_fit.training_feature_scales,
                "training scaler was not shared"
            );
        }
        ensure!(
            fit.shared_dna_fit.training_rows == training.len(),
            "shared training row count differs"
        );
        let features: Vec<_> = training
            .iter()
            .map(|i| panel.features()[*i].vector())
            .collect();
        // Validate vector shapes before indexing any saved feature statistic.
        fit.shared_dna_fit.effect(&features[0])?;
        for specific in fit.independent_context_fits.values() {
            specific.effect(&features[0])?;
        }
        for j in 0..33 {
            let expected_mean = mean(&features.iter().map(|row| row[j]).collect::<Vec<_>>())?;
            let sd = mean(
                &features
                    .iter()
                    .map(|row| (row[j] - expected_mean).powi(2))
                    .collect::<Vec<_>>(),
            )?
            .sqrt();
            close(
                fit.shared_dna_fit.training_feature_means[j],
                expected_mean,
                1e-7,
                "training feature mean",
            )?;
            close(
                fit.shared_dna_fit.training_feature_scales[j],
                if sd <= 1e-12 { 1.0 } else { sd },
                1e-7,
                "training population feature SD",
            )?;
            close(
                fit.shared_dna_fit.standardized_coefficients[j],
                mean(
                    &fit.independent_context_fits
                        .values()
                        .map(|ridge| ridge.standardized_coefficients[j])
                        .collect::<Vec<_>>(),
                )?,
                1e-7,
                "shared/independent coefficient linearity",
            )?;
        }
        fits.push(fit);
    }
    for prediction in predictions {
        let fit = &fits[prediction.fold_index_0based];
        let row = panel.features()[panel.sequence_index(&prediction.sequence_id)?].vector();
        let shared_effect = fit.shared_dna_fit.effect(&row)?;
        let specific = &fit.independent_context_fits[&prediction.context];
        let expected = [
            fit.global_mean,
            fit.training_context_means[&prediction.context],
            fit.global_mean + shared_effect,
            fit.training_context_means[&prediction.context] + shared_effect,
            specific.intercept + specific.effect(&row)?,
        ];
        for (observed, expected) in prediction.predictions.values().iter().zip(expected) {
            close(*observed, expected, 1e-7, "saved fit/prediction")?;
        }
    }
    Ok(())
}

pub fn audit(root: &Path, artifacts: &super::Artifacts) -> Result<()> {
    let manifest = load(&root.join("source-manifest.json"))?;
    ensure!(
        artifacts.files.keys().cloned().collect::<Vec<_>>() == inventory_names(&manifest)?,
        "incomplete activity artifact inventory"
    );
    let canonical = root.canonicalize()?;
    for (name, expected) in &artifacts.files {
        let path = root.join(name).canonicalize()?;
        ensure!(
            path.starts_with(&canonical),
            "artifact path escapes run directory"
        );
        checked_hash(&path, expected)?;
    }
    checked_hash(&root.join("protocol.json"), PROTOCOL_SHA)?;
    let inputs = load(&root.join("inputs.json"))?;
    source_gates(&manifest)?;
    source_file_inventory(root, &manifest)?;
    check_source_receipts(root, &manifest)?;
    ensure!(
        inputs["protocol_sha256"] == PROTOCOL_SHA
            && inputs["endpoint"] == ENDPOINT
            && inputs["source_manifest_sha256"]
                == genome_data::sha256(&fs::read(root.join("source-manifest.json"))?)
            && inputs["grouping_sha256"] == GROUPS_SHA
            && inputs["dna_panel_sha256"] == DNA_SHA
            && inputs["raw_source_sha256"] == RAW_SHA
            && inputs["activity_sha256"] == manifest["files"]["activity.json.gz"],
        "activity input/source commitments differ"
    );
    let dataset = PathBuf::from(
        inputs["source_dataset_path"]
            .as_str()
            .context("missing recorded source dataset path")?,
    );
    let raw = checked_hash(
        &dataset.join(
            manifest["source"]["raw_path"]
                .as_str()
                .context("missing raw source path")?,
        ),
        RAW_SHA,
    )?;
    ensure!(raw.len() == 53_557_805, "raw source bytes differ");
    drop(raw);
    let dna = checked_hash(
        Path::new(
            inputs["dna_panel_path"]
                .as_str()
                .context("missing recorded DNA path")?,
        ),
        DNA_SHA,
    )?;
    let groups = checked_hash(&root.join("grouping.json"), GROUPS_SHA)?;
    let activity: Activity =
        serde_json::from_slice(&uncompressed(&root.join("activity.json.gz"))?)?;
    let (panel, folds, outcomes) = panel_and_folds(&dna, &groups, &activity)?;
    verify_export_and_types(root, &manifest, &activity)?;
    let revisions: Vec<Revision> = serde_json::from_slice(&fs::read(root.join("ledger.json"))?)?;
    Ledger::verify(&revisions)?;
    let summary = load(&root.join("summary.json"))?;
    ledger_binding(&revisions, "freeze_activity_protocol_and_inputs", &inputs)?;
    ledger_binding(&revisions, "validate_activity_source_gates", &manifest)?;
    ledger_binding(&revisions, "execute_fixed_activity_evaluation", &summary)?;
    ensure!(
        summary["endpoint"] == ENDPOINT
            && summary["protocol_sha256"] == PROTOCOL_SHA
            && summary["source_gate_passed"] == true
            && summary["primary_comparison_index_0based"] == 1
            && summary["aggregation_status"] == manifest["aggregation_status"]
            && summary["source_schema_binding"] == manifest["schema_binding"],
        "activity summary semantic scope differs"
    );
    let full = load_gzip(&root.join("study.json.gz"))?;
    let resampling = load_gzip(&root.join("resampling.json.gz"))?;
    let common = load_gzip(&root.join("common-activity.json.gz"))?;
    ensure!(
        load_gzip(&root.join("features.json.gz"))? == feature_dump(&panel),
        "saved DNA feature fractions/counts differ"
    );
    let commitments = json!({"encoding":"canonical serde_json Value pretty JSON","study_sha256":canonical_sha(&full)?,
        "resampling_sha256":canonical_sha(&resampling)?,"common_activity_sha256":canonical_sha(&common)?});
    ledger_binding(
        &revisions,
        "commit_activity_numeric_artifacts",
        &commitments,
    )?;
    for key in [
        "panel_metadata_sha256",
        "outcomes_sha256",
        "sequence_count",
        "group_count",
        "measurement_count",
        "contexts",
        "alpha",
        "model_order",
        "models",
        "folds",
        "balanced_design_checks",
    ] {
        ensure!(
            summary[key] == full[key],
            "activity summary/full binding differs: {key}"
        );
    }
    ensure!(
        full["panel_metadata_sha256"] == panel.metadata_sha256()
            && full["outcomes_sha256"] == outcomes.outcomes_sha256()
            && inputs["panel_metadata_sha256"] == panel.metadata_sha256()
            && inputs["outcomes_sha256"] == outcomes.outcomes_sha256(),
        "typed source/feature/outcome commitment differs"
    );
    ensure!(
        full["sequence_count"] == 55_276
            && full["group_count"] == 2761
            && full["measurement_count"] == 165_828
            && full["alpha"] == 10.0
            && full["contexts"] == json!(CONTEXTS)
            && full["model_order"] == json!(ModelKind::ALL),
        "frozen activity model/count contract differs"
    );
    check_balanced(&full["balanced_design_checks"])?;
    let predictions: Vec<SavedPrediction> = serde_json::from_value(full["predictions"].clone())?;
    check_prediction_inputs(&predictions, &panel, &folds, &outcomes)?;
    check_saved_fits(&full, &predictions, &panel, &folds, &outcomes)?;
    let recalculated = common_activity(&predictions)?;
    ensure!(
        recalculated == common,
        "common activity diagnostic differs from held-out predictions"
    );
    let mut common_summary = recalculated;
    common_summary.as_object_mut().unwrap().remove("records");
    ensure!(
        summary["common_activity"] == common_summary
            && summary["fold_metrics"] == fold_metrics(&predictions)?,
        "activity mean/fold summary differs"
    );
    ensure!(
        summary["planned_comparisons"] == resampling["comparisons"]
            && resampling["replicates"] == 2000
            && resampling["seed"] == 0x2026102702_u64
            && resampling["weighting"] == "equal_sequences"
            && resampling["interval_probabilities"] == json!([0.025, 0.975])
            && resampling["model_order"] == json!(ModelKind::ALL),
        "bootstrap contract differs"
    );
    let models = full["models"]
        .as_array()
        .context("missing model summaries")?;
    ensure!(models.len() == 5, "frozen five model summaries required");
    for (j, model) in models.iter().enumerate() {
        let mse = mean(
            &predictions
                .iter()
                .map(|row| (row.observed - row.predictions.values()[j]).powi(2))
                .collect::<Vec<_>>(),
        )?;
        let mae = mean(
            &predictions
                .iter()
                .map(|row| (row.observed - row.predictions.values()[j]).abs())
                .collect::<Vec<_>>(),
        )?;
        close(
            model["sequence_weighted_mse"]
                .as_f64()
                .context("missing model MSE")?,
            mse,
            1e-8,
            "modelMSE/predictions",
        )?;
        close(
            model["sequence_weighted_mae"]
                .as_f64()
                .context("missing model MAE")?,
            mae,
            1e-8,
            "modelMAE/predictions",
        )?;
        close(
            resampling["observed_model_mse"][j]
                .as_f64()
                .context("missing bootstrap model MSE")?,
            mse,
            1e-8,
            "bootstrap/model MSE",
        )?;
        close(
            resampling["observed_model_mae"][j]
                .as_f64()
                .context("missing bootstrap model MAE")?,
            mae,
            1e-8,
            "bootstrap/model MAE",
        )?;
    }
    close(
        resampling["comparisons"][1]["observed_mse_improvement"]
            .as_f64()
            .context("missing primary improvement")?,
        common["mean_activity_mse_improvement"]
            .as_f64()
            .context("missing common gain")?,
        1e-8,
        "primary/common gain",
    )?;
    println!("Activity-v3 artifact/source/ledger commitments, frozen folds, exact observed targets, loss summaries and common-activity identity verified. This audit does not independently refit models or establish biological validation.");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn activity() -> Activity {
        Activity {
            protocol_sha256: PROTOCOL_SHA.into(),
            endpoint: ENDPOINT.into(),
            contexts: CONTEXTS.iter().map(|cell| (*cell).into()).collect(),
            records: vec![
                ActivityRecord {
                    id: "a".into(),
                    observations: BTreeMap::from([
                        ("HepG2".into(), 1.25),
                        ("K562".into(), -2.5),
                        ("WTC11".into(), 0.0),
                    ]),
                },
                ActivityRecord {
                    id: "b".into(),
                    observations: BTreeMap::from([
                        ("HepG2".into(), 2.0),
                        ("K562".into(), 3.0),
                        ("WTC11".into(), 4.0),
                    ]),
                },
            ],
        }
    }
    fn export() -> Value {
        json!({"protocol_sha256":RECOVERY_SHA,"source":{"sha256":RAW_SHA},"tables":{"all_cell_types_summary":{
        "id_header":"name","headers":{"1":"name","2":"HepG2 [log2(rna/dna)]","3":"K562 [log2(rna/dna)]","4":"WTC11 [log2(rna/dna)]"},
        "rows":[{"id":"a","values":{"HepG2 [log2(rna/dna)]":"1.25","K562 [log2(rna/dna)]":"-2.5","WTC11 [log2(rna/dna)]":"0"}},
            {"id":"b","values":{"HepG2 [log2(rna/dna)]":"2","K562 [log2(rna/dna)]":"3","WTC11 [log2(rna/dna)]":"4"}}]}}})
    }
    #[test]
    fn activity_adapter_rejects_missing_extra_duplicate_and_nonfinite_targets() {
        let ids = BTreeSet::from(["a".into(), "b".into()]);
        assert_eq!(validate_activity(&activity(), &ids).unwrap().len(), 6);
        let mut missing = activity();
        missing.records[0].observations.remove("WTC11");
        assert!(validate_activity(&missing, &ids).is_err());
        let mut extra = activity();
        extra.records[1].id = "other".into();
        assert!(validate_activity(&extra, &ids).is_err());
        let mut duplicate = activity();
        duplicate.records[1].id = "a".into();
        assert!(validate_activity(&duplicate, &ids).is_err());
        let mut nonfinite = activity();
        nonfinite.records[0]
            .observations
            .insert("K562".into(), f64::NAN);
        assert!(validate_activity(&nonfinite, &ids).is_err());
        let mut wrong = activity();
        wrong.contexts[0] = "different cell".into();
        assert!(validate_activity(&wrong, &ids).is_err());
    }
    #[test]
    fn source_export_is_exactly_bound_without_transforming_published_values() {
        verify_export_values(&export(), &activity()).unwrap();
        let mut changed = export();
        changed["tables"]["all_cell_types_summary"]["rows"][0]["values"]["HepG2 [log2(rna/dna)]"] =
            json!("1.25000001");
        assert!(verify_export_values(&changed, &activity()).is_err());
        changed = export();
        changed["tables"]["all_cell_types_summary"]["rows"][0]["values"]["HepG2 [log2(rna/dna)]"] =
            Value::Null;
        assert!(verify_export_values(&changed, &activity()).is_err());
    }
    #[test]
    fn source_gate_failure_blocks_incomplete_schema_and_aggregation_status() {
        let mut manifest = json!({"protocol_sha256":PROTOCOL_SHA,"endpoint":ENDPOINT,"source":{"sha256":RAW_SHA,"bytes":53557805},
            "gates":{"archive_integrity_verified":true,"empirical_activity_schema_bound":true,"coverage_complete":true,"design_identity_verified":true,"replicate_aggregation_status_recorded":true},
            "schema_binding":{"literal_header_bound":true},"aggregation_status":{"status":"reported_uncertainty_not_a_favorable_pass"}});
        source_gates(&manifest).unwrap();
        manifest["gates"]["coverage_complete"] = json!(false);
        assert!(source_gates(&manifest).is_err());
    }
    #[test]
    fn common_activity_identity_and_tampered_observation_binding_are_checked() {
        let sequences = ["AAAAAA", "CCCCCC", "ACACAC", "AGAGAG", "ATATAT"];
        let panel = MatchedPanel::new(
            sequences
                .iter()
                .enumerate()
                .map(|(i, sequence)| SequenceRecord {
                    id: format!("s{i}"),
                    group: format!("g{i}"),
                    sequence: (*sequence).into(),
                })
                .collect(),
            CONTEXTS.iter().map(|cell| (*cell).into()).collect(),
        )
        .unwrap();
        let rows = panel
            .sequences()
            .iter()
            .enumerate()
            .flat_map(|(i, sequence)| {
                CONTEXTS
                    .iter()
                    .enumerate()
                    .map(move |(j, context)| MeasuredOutcome {
                        sequence_id: sequence.id.clone(),
                        context: (*context).into(),
                        value: (i * i + 2 * j) as f64 - (i * j) as f64,
                    })
            })
            .collect::<Vec<_>>();
        let outcomes = OutcomeMatrix::from_records(&panel, &rows).unwrap();
        let folds = GroupSplit::partition(
            &panel,
            &(0..5).map(|i| vec![format!("g{i}")]).collect::<Vec<_>>(),
        )
        .unwrap();
        let evaluation = evaluate(&panel, &folds, &outcomes, 10.0).unwrap();
        let predictions = native_predictions(&evaluation);
        check_prediction_inputs(&predictions, &panel, &folds, &outcomes).unwrap();
        let mut saved = serde_json::to_value(&evaluation).unwrap();
        check_saved_fits(&saved, &predictions, &panel, &folds, &outcomes).unwrap();
        saved["folds"][0]["fit"]["shared_dna_fit"]["standardized_coefficients"][0] = json!(999.0);
        assert!(check_saved_fits(&saved, &predictions, &panel, &folds, &outcomes).is_err());
        let diagnostic = common_activity(&predictions).unwrap();
        assert!(
            diagnostic["maximum_absolute_per_element_identity_residual"]
                .as_f64()
                .unwrap()
                < 1e-8
        );
        let mut tampered = predictions;
        tampered[0].observed += 1.0;
        assert!(check_prediction_inputs(&tampered, &panel, &folds, &outcomes).is_err());
    }
    #[test]
    fn canonical_numeric_commitment_and_ledger_tamper_binding() {
        #[derive(Serialize)]
        struct DifferentOrder {
            z: u64,
            a: u64,
        }
        let payload = DifferentOrder { z: 7, a: 3 };
        assert_eq!(
            canonical_sha(&payload).unwrap(),
            canonical_sha(&json!({"a":3,"z":7})).unwrap()
        );
        let mut ledger = Ledger::default();
        ledger.append("numeric", json!({"result":3})).unwrap();
        ledger_binding(ledger.revisions(), "numeric", &json!({"result":3})).unwrap();
        assert!(ledger_binding(ledger.revisions(), "numeric", &json!({"result":4})).is_err());
    }

    #[test]
    fn dynamic_source_receipts_are_bound_and_cannot_overwrite_result_files() {
        let mut files = serde_json::Map::new();
        for name in SOURCE_FILES {
            files.insert(name.into(), json!("0".repeat(64)));
        }
        files.insert(
            "additional-source-receipt.json".into(),
            json!("1".repeat(64)),
        );
        let mut manifest = json!({"files":files});
        assert!(inventory_names(&manifest)
            .unwrap()
            .contains(&"additional-source-receipt.json".into()));
        manifest["files"]["summary.json"] = json!("1".repeat(64));
        assert!(inventory_names(&manifest).is_err());
    }
}
