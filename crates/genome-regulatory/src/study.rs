//! Frozen, retrospective regulatory-prediction experiment.
use crate::{
    math::{
        regression_metrics, regression_metrics_with_observed_ranks, RegressionMetrics, RidgeDesign,
        StandardizedRidge,
    },
    motif::{MotifEffect, MotifSet},
};
use anyhow::{ensure, Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::Path,
};

const FROZEN_PROTOCOL_SHA256: &str =
    "c8e97fabef2afe8ccaa5a01a7ef7bb952f9c9f15857e5c540af112a66e78cb58";
const HBB: &str = "urn:mavedb:00000018-a-1";
pub const MODEL_NAMES: [&str; 6] = [
    "mean",
    "chemistry",
    "coordinate",
    "delta3mer",
    "motif",
    "full",
];
const COMPARISONS: [(&str, usize, usize); 3] = [
    ("delta3mer_vs_chemistry", 3, 1),
    ("motif_vs_chemistry", 4, 1),
    ("full_vs_coordinate", 5, 2),
];
const MUTATIONS: [&str; 12] = [
    "A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G",
];

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct InputRecord {
    pub id: String,
    pub assay_id: String,
    pub position_1based: usize,
    pub reference: String,
    pub alternate: String,
    pub target_sequence: String,
    pub score: Option<f64>,
    /// Original score-column decimal, retained for exact observed rank ties.
    pub score_raw: String,
    pub unique_tags: Option<f64>,
    pub dna_sequences: Option<f64>,
    pub rna_sequences: Option<f64>,
    pub cell_context: String,
}

pub struct Protocol {
    json: Value,
    sha256: String,
}
impl Protocol {
    pub fn open(path: &Path) -> Result<Self> {
        Self::from_bytes(&fs::read(path)?)
    }
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let sha256 = digest(bytes);
        ensure!(
            sha256 == FROZEN_PROTOCOL_SHA256,
            "protocol differs from the final pre-acquisition frozen v1 definition"
        );
        let json: Value = serde_json::from_slice(bytes)?;
        ensure!(
            json["schema_version"] == 1
                && json["status"] == "frozen"
                && json["protocol_id"] == "regulatory-prediction-v1",
            "unsupported regulatory protocol"
        );
        Ok(Self { json, sha256 })
    }
    pub fn as_value(&self) -> &Value {
        &self.json
    }
    pub fn sha256(&self) -> &str {
        &self.sha256
    }
    fn primary_ids(&self) -> Vec<String> {
        self.json["cohort"]["assay_roles"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|row| row["role"] == "primary")
            .map(|row| row["assay_id"].as_str().unwrap().to_owned())
            .collect()
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Prediction {
    pub id: String,
    pub assay_id: String,
    pub position_1based: usize,
    pub reference: String,
    pub alternate: String,
    pub fold_index_0based: usize,
    pub observed: f64,
    pub observed_source_decimal: String,
    pub predicted: f64,
    pub residual_observed_minus_predicted: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct FoldDiagnostic {
    pub family: String,
    pub fold_index_0based: usize,
    pub target_assay: String,
    pub test_interval_0based_half_open: Option<[usize; 2]>,
    pub guard_bp: usize,
    pub train_ids: Vec<String>,
    pub test_ids: Vec<String>,
    pub excluded_guard_ids: Vec<String>,
    pub training_assays: Vec<String>,
    pub training_cell_contexts: Vec<String>,
    pub test_cell_contexts: Vec<String>,
    pub minimum_same_assay_train_test_center_distance: Option<usize>,
    pub shared_train_test_positions: usize,
    pub test_rows_with_overlapping_same_assay_training_feature_support: usize,
    pub strict_same_assay_feature_support_disjoint: bool,
    pub contextual_rule: String,
}

struct Plan {
    diagnostic: FoldDiagnostic,
    training: Vec<usize>,
    test: Vec<usize>,
}
struct PreparedPlan {
    plan: Plan,
    designs: Vec<RidgeDesign>,
    /// Identical UNSCALED feature rows share one evaluated prediction per fit.
    test_groups: Vec<Vec<Vec<usize>>>,
}

#[derive(Clone, Debug, Serialize)]
pub struct FittedFold {
    pub fold_index_0based: usize,
    pub fit: StandardizedRidge,
}

#[derive(Clone, Debug, Serialize)]
pub struct PositionError {
    pub assay_id: String,
    pub position_1based: usize,
    pub allele_count: usize,
    pub observed_mean: f64,
    pub observed_mean_exact_numerator: String,
    pub observed_mean_exact_denominator: String,
    pub predicted_mean: f64,
    pub mean_squared_error: f64,
    pub mean_absolute_error: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct AssayMetrics {
    pub assay_id: String,
    pub cell_context: String,
    pub target_rows: usize,
    pub target_positions: usize,
    pub position_weighted_mse: f64,
    pub position_weighted_mae: f64,
    pub row_metrics: RegressionMetrics,
    pub position_mean_metrics: RegressionMetrics,
    pub relative_mse_versus_training_mean: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct ModelResult {
    pub model: String,
    pub feature_names: Vec<String>,
    pub fits: Vec<FittedFold>,
    pub predictions: Vec<Prediction>,
    pub per_position_errors: Vec<PositionError>,
    pub per_assay_metrics: Vec<AssayMetrics>,
    pub pooled_row_metrics: RegressionMetrics,
    pub macro_assay_position_mse: f64,
    pub macro_assay_position_mae: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct AnalysisResult {
    pub family: String,
    pub folds: Vec<FoldDiagnostic>,
    pub models: Vec<ModelResult>,
    pub interpretation: String,
}

#[derive(Serialize)]
pub struct StudyResult {
    pub schema_version: u32,
    pub protocol_sha256: String,
    pub protocol: Value,
    pub input_records_sha256: String,
    pub motif_manifest_sha256: String,
    pub metric_arithmetic: Value,
    pub inclusion: Value,
    pub feature_schema: BTreeMap<String, Vec<String>>,
    pub motif_features: Vec<Value>,
    pub analyses: Vec<AnalysisResult>,
    pub resampling: Value,
    pub spatial_null: Value,
    pub limitations: Vec<String>,
}

fn prepare_inputs(input: &[InputRecord], protocol: &Protocol) -> Result<(Vec<InputRecord>, Value)> {
    ensure!(
        !input.is_empty() && input.len() <= 30_000,
        "empty input or regulatory row budget exceeded"
    );
    let roles: BTreeMap<String, String> = protocol.json["cohort"]["assay_roles"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            (
                row["assay_id"].as_str().unwrap().into(),
                row["role"].as_str().unwrap().into(),
            )
        })
        .collect();
    let mut ids = BTreeSet::new();
    let mut variants = BTreeSet::new();
    let mut assay_context = BTreeMap::<String, (String, String)>::new();
    let mut included = Vec::new();
    let mut exclusions = Vec::new();
    let mut exclusion_counts = BTreeMap::<String, usize>::new();
    let mut assay_input_counts = BTreeMap::<String, usize>::new();
    for record in input {
        ensure!(
            roles.contains_key(&record.assay_id),
            "input assay is not in frozen protocol"
        );
        ensure!(
            !record.id.is_empty() && ids.insert(record.id.clone()),
            "empty or duplicate regulatory record ID"
        );
        ensure!(
            !record.cell_context.trim().is_empty(),
            "missing assay cell context"
        );
        let seq = record.target_sequence.as_bytes();
        ensure!(
            seq.len() >= 11 && seq.len() <= 5000 && seq.iter().all(|b| b"ACGT".contains(b)),
            "unsupported target sequence/alphabet/length"
        );
        ensure!(
            (1..=seq.len()).contains(&record.position_1based),
            "target variant position out of bounds"
        );
        let reference = record.reference.as_bytes();
        let alternate = record.alternate.as_bytes();
        ensure!(
            reference.len() == 1
                && alternate.len() == 1
                && b"ACGT".contains(&reference[0])
                && b"ACGT".contains(&alternate[0])
                && reference[0] != alternate[0],
            "record is not an unambiguous unequal SNV"
        );
        ensure!(
            seq[record.position_1based - 1] == reference[0],
            "regulatory variant reference mismatch"
        );
        ensure!(
            variants.insert((
                record.assay_id.clone(),
                record.position_1based,
                record.alternate.clone()
            )),
            "duplicate regulatory SNV"
        );
        if let Some((sequence, cell)) = assay_context.get(&record.assay_id) {
            ensure!(
                sequence == &record.target_sequence && cell == &record.cell_context,
                "assay has inconsistent reference/cell context"
            );
        } else {
            assay_context.insert(
                record.assay_id.clone(),
                (record.target_sequence.clone(), record.cell_context.clone()),
            );
        }
        for value in [
            record.score,
            record.unique_tags,
            record.dna_sequences,
            record.rna_sequences,
        ]
        .into_iter()
        .flatten()
        {
            ensure!(value.is_finite(),"nonfinite regulatory numeric input must be explicitly represented as missing with source-status accounting");
        }
        *assay_input_counts
            .entry(record.assay_id.clone())
            .or_default() += 1;
        let mut reasons = Vec::new();
        if record.score.is_none() {
            reasons.push("missing_or_invalid_score");
        }
        if record
            .unique_tags
            .is_none_or(|tags| tags < 10.0 || tags.fract() != 0.0)
        {
            reasons.push("missing_invalid_or_fewer_than_10_unique_tags");
        }
        if record
            .dna_sequences
            .is_none_or(|count| count <= 0.0 || count.fract() != 0.0)
        {
            reasons.push("missing_invalid_or_nonpositive_dna_count");
        }
        if record
            .rna_sequences
            .is_none_or(|count| count <= 0.0 || count.fract() != 0.0)
        {
            reasons.push("missing_invalid_or_nonpositive_rna_count");
        }
        if reasons.is_empty() {
            included.push(record.clone());
        } else {
            for reason in &reasons {
                *exclusion_counts.entry((*reason).into()).or_default() += 1;
            }
            exclusions.push(json!({"id":record.id,"assay_id":record.assay_id,"position_1based":record.position_1based,"reasons":reasons}));
        }
    }
    included.sort_by(|a, b| {
        (&a.assay_id, a.position_1based, &a.alternate, &a.id).cmp(&(
            &b.assay_id,
            b.position_1based,
            &b.alternate,
            &b.id,
        ))
    });
    let mut per_assay = Vec::new();
    for (assay, (sequence, cell)) in assay_context {
        let rows: Vec<_> = included.iter().filter(|r| r.assay_id == assay).collect();
        let mut position_counts = BTreeMap::<usize, usize>::new();
        for row in &rows {
            *position_counts.entry(row.position_1based).or_default() += 1;
        }
        ensure!(
            position_counts.values().all(|&n| n <= 3),
            "more than three alternate alleles at a position"
        );
        let all_input_positions: BTreeSet<_> = input
            .iter()
            .filter(|r| r.assay_id == assay)
            .map(|r| r.position_1based)
            .collect();
        per_assay.push(json!({"assay_id":assay,"role":roles[&assay],"cell_context":cell,"target_length":sequence.len(),
            "input_supported_snv_rows":assay_input_counts[&assay],"included_rows":rows.len(),"included_positions":position_counts.len(),
            "positions_with_one_allele":position_counts.values().filter(|&&n|n==1).count(),
            "positions_with_two_alleles":position_counts.values().filter(|&&n|n==2).count(),
            "positions_with_three_alleles":position_counts.values().filter(|&&n|n==3).count(),
            "input_positions_with_no_eligible_alleles":all_input_positions.len()-position_counts.len(),
            "included_ids":rows.iter().map(|r|r.id.as_str()).collect::<Vec<_>>()}));
    }
    for assay in protocol.primary_ids().into_iter().chain([HBB.to_owned()]) {
        ensure!(
            included.iter().any(|r| r.assay_id == assay),
            "required evaluation assay has no eligible rows: {assay}"
        );
    }
    let report = json!({"input_supported_snv_rows":input.len(),"included_rows":included.len(),"excluded_rows":exclusions.len(),
        "exclusion_reason_counts":exclusion_counts,"exclusions":exclusions,"per_assay":per_assay,
        "incomplete_position_policy":"Each quality-passing allele retained; missing siblings do not discard valid rows",
        "source_scope":"Input records are already reference-validated SNVs; caller must retain raw non-SNV/source-status accounting and source manifest commitments"});
    Ok((included, report))
}

type FeatureRows = Vec<[Vec<f64>; 6]>;
fn features(record: &InputRecord, motifs: &MotifSet) -> Result<([Vec<f64>; 6], Vec<MotifEffect>)> {
    let sequence = record.target_sequence.as_bytes();
    let position = record.position_1based - 1;
    let mutation = format!("{}>{}", record.reference, record.alternate);
    let mutation_index = MUTATIONS
        .iter()
        .position(|&m| m == mutation)
        .context("unknown substitution identity")?;
    let mut chemistry = vec![0.0; 12];
    chemistry[mutation_index] = 1.0;
    let x = 2.0 * position as f64 / (sequence.len() - 1) as f64 - 1.0;
    let coordinate = vec![x, x * x];
    let mut delta3mer = vec![0.0; 64];
    let alternate = record.alternate.as_bytes()[0];
    for start in position.saturating_sub(2)..=position.min(sequence.len() - 3) {
        let original = &sequence[start..start + 3];
        let mut changed = [original[0], original[1], original[2]];
        changed[position - start] = alternate;
        delta3mer[kmer_index(original)?] -= 1.0;
        delta3mer[kmer_index(&changed)?] += 1.0;
    }
    let effects: Vec<_> = motifs
        .motifs
        .iter()
        .map(|motif| {
            motif
                .effect(
                    &record.target_sequence,
                    record.position_1based,
                    &record.alternate,
                )?
                .context("no complete motif context")
        })
        .collect::<Result<_>>()?;
    let motif_deltas: Vec<_> = effects
        .iter()
        .map(|effect| effect.delta_max_log2_odds)
        .collect();
    let concat = |parts: &[&[f64]]| {
        parts
            .iter()
            .flat_map(|part| part.iter().copied())
            .collect::<Vec<_>>()
    };
    Ok((
        [
            Vec::new(),
            chemistry.clone(),
            concat(&[&chemistry, &coordinate]),
            concat(&[&chemistry, &delta3mer]),
            concat(&[&chemistry, &motif_deltas]),
            concat(&[&chemistry, &coordinate, &delta3mer, &motif_deltas]),
        ],
        effects,
    ))
}

fn feature_schema() -> BTreeMap<String, Vec<String>> {
    let chemistry: Vec<_> = MUTATIONS.iter().map(|m| format!("mutation:{m}")).collect();
    let coordinate = vec!["position:x".into(), "position:x_squared".into()];
    let mut kmer = Vec::new();
    for a in b"ACGT" {
        for b in b"ACGT" {
            for c in b"ACGT" {
                kmer.push(format!(
                    "delta3mer:{}{}{}",
                    *a as char, *b as char, *c as char
                ));
            }
        }
    }
    let motif: Vec<String> = ["MA0035.4", "MA0493.3", "MA0060.3", "MA0108.3"]
        .iter()
        .map(|id| format!("motif_delta:{id}"))
        .collect();
    let concat = |parts: &[&[String]]| {
        parts
            .iter()
            .flat_map(|part| part.iter().cloned())
            .collect::<Vec<_>>()
    };
    [
        ("mean", Vec::new()),
        ("chemistry", chemistry.clone()),
        ("coordinate", concat(&[&chemistry, &coordinate])),
        ("delta3mer", concat(&[&chemistry, &kmer])),
        ("motif", concat(&[&chemistry, &motif])),
        ("full", concat(&[&chemistry, &coordinate, &kmer, &motif])),
    ]
    .into_iter()
    .map(|(id, names)| (id.into(), names))
    .collect()
}

fn kmer_index(kmer: &[u8]) -> Result<usize> {
    let mut index = 0;
    for base in kmer {
        index = index * 4
            + b"ACGT"
                .iter()
                .position(|b| b == base)
                .context("invalid kmer base")?;
    }
    Ok(index)
}

fn make_plan(
    specification: (&str, usize, &str, Option<[usize; 2]>, usize),
    training: Vec<usize>,
    test: Vec<usize>,
    excluded: Vec<usize>,
    records: &[InputRecord],
) -> Result<Plan> {
    let (family, fold, target, interval, guard) = specification;
    ensure!(
        !training.is_empty() && !test.is_empty(),
        "empty training/test fold: {family}/{target}/{fold}"
    );
    let train_set: BTreeSet<_> = training.iter().copied().collect();
    ensure!(
        test.iter().all(|i| !train_set.contains(i)),
        "training/test row overlap"
    );
    let mut min_distance = None;
    let mut shared_positions = BTreeSet::new();
    let mut overlapping_test_rows = 0;
    for &test_index in &test {
        let b = &records[test_index];
        let mut overlaps = false;
        for &train_index in &training {
            let a = &records[train_index];
            if a.assay_id == b.assay_id {
                let distance = a.position_1based.abs_diff(b.position_1based);
                min_distance =
                    Some(min_distance.map_or(distance, |prior: usize| prior.min(distance)));
                if distance == 0 {
                    shared_positions.insert(b.position_1based);
                }
                if distance <= 20 {
                    overlaps = true;
                }
            }
        }
        if overlaps {
            overlapping_test_rows += 1;
        }
    }
    if family == "hbb_guarded_position_cv" {
        ensure!(
            min_distance.is_some_and(|d| d > 20) && overlapping_test_rows == 0,
            "guard does not isolate radius10 sequence supports"
        );
    }
    let strings = |indices: &[usize]| {
        indices
            .iter()
            .map(|&i| records[i].id.clone())
            .collect::<Vec<_>>()
    };
    let unique = |indices: &[usize], cell: bool| {
        indices
            .iter()
            .map(|&i| {
                if cell {
                    records[i].cell_context.clone()
                } else {
                    records[i].assay_id.clone()
                }
            })
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>()
    };
    Ok(Plan{diagnostic:FoldDiagnostic{
        family:family.into(),fold_index_0based:fold,target_assay:target.into(),
        test_interval_0based_half_open:interval,guard_bp:guard,
        train_ids:strings(&training),test_ids:strings(&test),excluded_guard_ids:strings(&excluded),
        training_assays:unique(&training,false),training_cell_contexts:unique(&training,true),test_cell_contexts:unique(&test,true),
        minimum_same_assay_train_test_center_distance:min_distance,
        shared_train_test_positions:shared_positions.len(),test_rows_with_overlapping_same_assay_training_feature_support:overlapping_test_rows,
        strict_same_assay_feature_support_disjoint:overlapping_test_rows==0,
        contextual_rule:match family{
            "primary_leave_element_out"=>"Training contains other primary-role loci only; validation-only/HBB excluded. Cell contexts can differ and are disclosed, not encoded as covariates.",
            "hbb_to_external_transfer"=>"All training is the previously exposed HBB source; target is a different primary locus/context; no test normalization.",
            "row_random_diagnostic"=>"Deliberately leaky HBB row split; never a primary or tuning criterion.",
            _=>"Same-target positions blocked; centers >20bp apart so radius10 feature supports are disjoint; assay/source dependence can remain."
        }.into()
    },training,test})
}

fn plans(records: &[InputRecord], protocol: &Protocol) -> Result<Vec<Plan>> {
    let mut plans = Vec::new();
    let hbb: Vec<_> = records
        .iter()
        .enumerate()
        .filter(|(_, r)| r.assay_id == HBB)
        .map(|(i, _)| i)
        .collect();
    let length = records[hbb[0]].target_sequence.len();
    for fold in 0..5 {
        let start = fold * length / 5;
        let end = (fold + 1) * length / 5;
        let test: Vec<_> = hbb
            .iter()
            .copied()
            .filter(|&i| (start..end).contains(&(records[i].position_1based - 1)))
            .collect();
        let centers: BTreeSet<_> = test.iter().map(|&i| records[i].position_1based).collect();
        let mut train = Vec::new();
        let mut excluded = Vec::new();
        for &i in &hbb {
            if centers.contains(&records[i].position_1based) {
                continue;
            }
            if centers
                .iter()
                .all(|&p| records[i].position_1based.abs_diff(p) > 20)
            {
                train.push(i);
            } else {
                excluded.push(i);
            }
        }
        plans.push(make_plan(
            ("hbb_guarded_position_cv", fold, HBB, Some([start, end]), 20),
            train,
            test,
            excluded,
            records,
        )?);
    }
    let primary: BTreeSet<_> = protocol.primary_ids().into_iter().collect();
    for (fold, target) in primary.iter().enumerate() {
        let test: Vec<_> = records
            .iter()
            .enumerate()
            .filter(|(_, r)| &r.assay_id == target)
            .map(|(i, _)| i)
            .collect();
        let train: Vec<_> = records
            .iter()
            .enumerate()
            .filter(|(_, r)| primary.contains(&r.assay_id) && &r.assay_id != target)
            .map(|(i, _)| i)
            .collect();
        ensure!(
            train.iter().all(|&i| records[i].assay_id != HBB
                && !protocol.json["cohort"]["assay_roles"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(
                        |a| a["assay_id"] == records[i].assay_id && a["role"] == "validation_only"
                    )),
            "context/role training leakage"
        );
        plans.push(make_plan(
            ("primary_leave_element_out", fold, target, None, 0),
            train,
            test.clone(),
            Vec::new(),
            records,
        )?);
        plans.push(make_plan(
            ("hbb_to_external_transfer", fold, target, None, 0),
            hbb.clone(),
            test,
            Vec::new(),
            records,
        )?);
    }
    for fold in 0..5 {
        let (test, train): (Vec<_>, Vec<_>) = hbb
            .iter()
            .copied()
            .partition(|&i| row_fold(&records[i]) == fold);
        plans.push(make_plan(
            ("row_random_diagnostic", fold, HBB, None, 0),
            train,
            test,
            Vec::new(),
            records,
        )?);
    }
    Ok(plans)
}

fn row_fold(record: &InputRecord) -> usize {
    let input = format!(
        "{}|{}|{}>{}|row-random-v1",
        record.assay_id, record.position_1based, record.reference, record.alternate
    );
    let hash = Sha256::digest(input.as_bytes());
    (u64::from_le_bytes(hash[..8].try_into().unwrap()) % 5) as usize
}

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

/// Execute exactly the frozen experiment. All feature choices precede fitting;
/// cached designs depend on training X only and every null fit recomputes ymean.
pub fn run(input: &[InputRecord], protocol: &Protocol, motifs: &MotifSet) -> Result<StudyResult> {
    let (records, inclusion) = prepare_inputs(input, protocol)?;
    let schema = feature_schema();
    let mut feature_rows = Vec::with_capacity(records.len());
    let mut motif_features = Vec::with_capacity(records.len());
    for record in &records {
        let (row, effects) = features(record, motifs)?;
        for (model, names) in MODEL_NAMES.iter().enumerate() {
            ensure!(
                row[model].len() == schema[*names].len(),
                "feature vocabulary/dimension mismatch"
            );
        }
        if let Some(score) = record.score {
            ensure!(
                record.score_raw.trim().parse::<f64>().ok() == Some(score),
                "source score decimal disagrees with parsed outcome"
            );
            crate::exact::decimal(&record.score_raw)?;
        }
        motif_features.push(json!({"id":record.id,"assay_id":record.assay_id,"position_1based":record.position_1based,"effects":effects}));
        feature_rows.push(row);
    }
    let labels: Vec<_> = records.iter().map(|r| r.score.unwrap()).collect();
    let fold_plans = plans(&records, protocol)?;
    let mut prepared = Vec::new();
    for plan in fold_plans {
        let mut designs = Vec::new();
        let mut test_groups = Vec::new();
        for model in 0..6 {
            let rows: Vec<_> = plan
                .training
                .iter()
                .map(|&i| feature_rows[i][model].clone())
                .collect();
            designs.push(RidgeDesign::prepare(&rows, 10.0)?);
            let mut groups = BTreeMap::<Vec<u64>, Vec<usize>>::new();
            for &index in &plan.test {
                let key = feature_rows[index][model]
                    .iter()
                    .map(|&x| if x == 0.0 { 0 } else { x.to_bits() })
                    .collect();
                groups.entry(key).or_default().push(index);
            }
            test_groups.push(groups.into_values().collect());
        }
        prepared.push(PreparedPlan {
            plan,
            designs,
            test_groups,
        });
    }
    let analyses = evaluate_real(&prepared, &feature_rows, &labels, &records, &schema)?;
    let resampling = resample_losses(&analyses)?;
    let spatial_null = spatial_control(
        &prepared,
        &feature_rows,
        &labels,
        &records,
        protocol,
        &analyses,
    )?;
    let mut canonical_input = input.to_vec();
    canonical_input.sort_by(|a, b| {
        (&a.assay_id, a.position_1based, &a.alternate, &a.id).cmp(&(
            &b.assay_id,
            b.position_1based,
            &b.alternate,
            &b.id,
        ))
    });
    Ok(StudyResult{
        schema_version:1,protocol_sha256:protocol.sha256.clone(),protocol:protocol.json.clone(),
        input_records_sha256:digest(&serde_json::to_vec(&canonical_input)?),
        motif_manifest_sha256:motifs.manifest_sha256().into(),
        metric_arithmetic:json!({
            "observed_position_means":"Exact rational arithmetic on original source CSV decimals; exact fractions determine observed position rank ties, and are converted to f64 only for numeric display and nonrank metrics",
            "predicted_values":"Actual f64 ridge predictions without quantization, epsilon ties or rounding; identical unscaled feature rows receive one evaluated prediction per fit",
            "constant_predictions":"Identical per-position predictions preserve the exact value before position-rank calculation",
            "unchanged_scientific_protocol":"Numerical correctness fixes preserve the frozen models/features/folds/alpha/controls; no outcome-driven parameter selection"
        }),
        inclusion,feature_schema:schema,motif_features,analyses,resampling,spatial_null,
        limitations:vec![
            "Retrospective selected public promoters from one study; neither random sampling nor independent laboratory replication is established".into(),
            "HBB outcomes were explored previously. Models, fixed alpha/features, final protocol and new cohort selection were fixed before the new score acquisition; no pristine-HBB or external preregistration claim".into(),
            "Training rows have equal weight, while errors average alleles per position and then loci; cell contexts and assay scales may differ and are never normalized using test outcomes".into(),
            "Blocked HBB centers have disjoint radius10 sequence supports, but source/library dependence and biology can remain. Row-random diagnostic is deliberately leaky and not a model-selection criterion".into(),
            "JASPAR log-odds deltas are sequence-similarity features, not TF binding or expression measurements; missing TBP taxonomy is retained in the motif source metadata".into(),
            "Bootstrap intervals condition on fixed OOF predictions. Artificial spatial rotations test a chosen null transformation and are not calibrated biological-significance p-values".into(),
            "Every planned model/control and negative result is retained; no discovery or causal mechanism is declared from thresholds".into()
        ],
    })
}

fn predict_groups(
    prepared: &PreparedPlan,
    fit: &StandardizedRidge,
    model: usize,
    features: &FeatureRows,
    total_rows: usize,
) -> Result<Vec<f64>> {
    let mut values = vec![0.0; total_rows];
    for group in &prepared.test_groups[model] {
        let first = group[0];
        let value = fit.predict(&features[first][model])?;
        for &index in group {
            values[index] = value;
        }
    }
    Ok(values)
}

fn evaluate_real(
    prepared: &[PreparedPlan],
    features: &FeatureRows,
    labels: &[f64],
    records: &[InputRecord],
    schema: &BTreeMap<String, Vec<String>>,
) -> Result<Vec<AnalysisResult>> {
    let families = [
        "hbb_guarded_position_cv",
        "primary_leave_element_out",
        "hbb_to_external_transfer",
        "row_random_diagnostic",
    ];
    let mut analyses = Vec::new();
    for family in families {
        let selected: Vec<_> = prepared
            .iter()
            .filter(|p| p.plan.diagnostic.family == family)
            .collect();
        let mut models = Vec::new();
        for (model, name) in MODEL_NAMES.iter().enumerate() {
            let mut predictions = Vec::new();
            let mut fits = Vec::new();
            let mut seen = BTreeSet::new();
            for prepared in &selected {
                let outcomes: Vec<_> = prepared.plan.training.iter().map(|&i| labels[i]).collect();
                let fit = prepared.designs[model].fit(&outcomes)?;
                let fitted_values = predict_groups(prepared, &fit, model, features, records.len())?;
                for &index in &prepared.plan.test {
                    let row = &records[index];
                    ensure!(
                        seen.insert(row.id.as_str()),
                        "target predicted more than once in evaluation family"
                    );
                    let predicted = fitted_values[index];
                    predictions.push(Prediction {
                        id: row.id.clone(),
                        assay_id: row.assay_id.clone(),
                        position_1based: row.position_1based,
                        reference: row.reference.clone(),
                        alternate: row.alternate.clone(),
                        fold_index_0based: prepared.plan.diagnostic.fold_index_0based,
                        observed: labels[index],
                        observed_source_decimal: row.score_raw.clone(),
                        predicted,
                        residual_observed_minus_predicted: labels[index] - predicted,
                    });
                }
                fits.push(FittedFold {
                    fold_index_0based: prepared.plan.diagnostic.fold_index_0based,
                    fit,
                });
            }
            predictions.sort_by(|a, b| {
                (&a.assay_id, a.position_1based, &a.alternate, &a.id).cmp(&(
                    &b.assay_id,
                    b.position_1based,
                    &b.alternate,
                    &b.id,
                ))
            });
            models.push(summarize_model(
                name,
                schema[*name].clone(),
                fits,
                predictions,
                records,
            )?);
        }
        let baseline: BTreeMap<_, _> = models[0]
            .per_assay_metrics
            .iter()
            .map(|m| (m.assay_id.clone(), m.position_weighted_mse))
            .collect();
        for model in &mut models {
            for assay in &mut model.per_assay_metrics {
                let denominator = baseline[&assay.assay_id];
                assay.relative_mse_versus_training_mean =
                    (denominator > 0.0).then(|| assay.position_weighted_mse / denominator);
            }
        }
        analyses.push(AnalysisResult{
            family:family.into(),folds:selected.iter().map(|p|p.plan.diagnostic.clone()).collect(),models,
            interpretation:match family{
                "hbb_guarded_position_cv"=>"Exploratory within-HBB OOF predictions on previously exposed outcomes; guard20 isolates radius10 sequence supports.",
                "primary_leave_element_out"=>"Nine held-out primary loci; training excludes that locus, all validation-only records and HBB. Report heterogeneous cell contexts and one-study dependence.",
                "hbb_to_external_transfer"=>"Secondary transport from already explored HBB to nine other elements without target normalization or tuning.",
                _=>"Explicitly leaky row diagnostic; shared positions/context are disclosed. Differences from blocked CV also reflect split difficulty and cannot be attributed solely to leakage."
            }.into()
        });
    }
    Ok(analyses)
}

fn summarize_model(
    name: &str,
    feature_names: Vec<String>,
    fits: Vec<FittedFold>,
    predictions: Vec<Prediction>,
    records: &[InputRecord],
) -> Result<ModelResult> {
    let mut positions = BTreeMap::<(String, usize), Vec<&Prediction>>::new();
    let mut assays = BTreeMap::<String, Vec<&Prediction>>::new();
    for prediction in &predictions {
        positions
            .entry((prediction.assay_id.clone(), prediction.position_1based))
            .or_default()
            .push(prediction);
        assays
            .entry(prediction.assay_id.clone())
            .or_default()
            .push(prediction);
    }
    let mut exact_position_means = BTreeMap::new();
    for (key, rows) in &positions {
        let values: Vec<_> = rows
            .iter()
            .map(|r| crate::exact::decimal(&r.observed_source_decimal))
            .collect::<Result<_>>()?;
        exact_position_means.insert(key.clone(), crate::exact::mean(&values)?);
    }
    let per_position_errors: Vec<_> = positions
        .into_iter()
        .map(|((assay_id, position_1based), rows)| {
            let n = rows.len() as f64;
            let exact_mean = &exact_position_means[&(assay_id.clone(), position_1based)];
            Ok(PositionError {
                assay_id,
                position_1based,
                allele_count: rows.len(),
                // Preserve true constants across groups with different allele
                // counts; otherwise one-ULP summation differences create fake
                // rank variation in a constant mean predictor.
                observed_mean: crate::exact::to_f64(exact_mean)?,
                observed_mean_exact_numerator: exact_mean.numer().to_string(),
                observed_mean_exact_denominator: exact_mean.denom().to_string(),
                predicted_mean: if rows.iter().all(|r| r.predicted == rows[0].predicted) {
                    rows[0].predicted
                } else {
                    rows.iter().map(|r| r.predicted / n).sum()
                },
                mean_squared_error: rows
                    .iter()
                    .map(|r| (r.observed - r.predicted).powi(2) / n)
                    .sum(),
                mean_absolute_error: rows
                    .iter()
                    .map(|r| (r.observed - r.predicted).abs() / n)
                    .sum(),
            })
        })
        .collect::<Result<_>>()?;
    let mut per_assay_metrics = Vec::new();
    for (assay, rows) in assays {
        let positions: Vec<_> = per_position_errors
            .iter()
            .filter(|p| p.assay_id == assay)
            .collect();
        let actual: Vec<_> = rows.iter().map(|r| r.observed).collect();
        let predicted: Vec<_> = rows.iter().map(|r| r.predicted).collect();
        let position_actual: Vec<_> = positions.iter().map(|p| p.observed_mean).collect();
        let position_predicted: Vec<_> = positions.iter().map(|p| p.predicted_mean).collect();
        let exact_means: Vec<_> = positions
            .iter()
            .map(|p| exact_position_means[&(assay.clone(), p.position_1based)].clone())
            .collect();
        let observed_ranks = crate::exact::ranks(&exact_means);
        let n = positions.len() as f64;
        per_assay_metrics.push(AssayMetrics {
            assay_id: assay.clone(),
            cell_context: records
                .iter()
                .find(|r| r.assay_id == assay)
                .unwrap()
                .cell_context
                .clone(),
            target_rows: rows.len(),
            target_positions: positions.len(),
            position_weighted_mse: positions.iter().map(|p| p.mean_squared_error / n).sum(),
            position_weighted_mae: positions.iter().map(|p| p.mean_absolute_error / n).sum(),
            row_metrics: regression_metrics(&actual, &predicted)?,
            position_mean_metrics: regression_metrics_with_observed_ranks(
                &position_actual,
                &position_predicted,
                &observed_ranks,
            )?,
            relative_mse_versus_training_mean: None,
        });
    }
    let actual: Vec<_> = predictions.iter().map(|r| r.observed).collect();
    let predicted: Vec<_> = predictions.iter().map(|r| r.predicted).collect();
    let n = per_assay_metrics.len() as f64;
    ensure!(n > 0.0, "evaluation has zero test assays");
    Ok(ModelResult {
        model: name.into(),
        feature_names,
        fits,
        predictions,
        per_position_errors,
        pooled_row_metrics: regression_metrics(&actual, &predicted)?,
        macro_assay_position_mse: per_assay_metrics
            .iter()
            .map(|a| a.position_weighted_mse / n)
            .sum(),
        macro_assay_position_mae: per_assay_metrics
            .iter()
            .map(|a| a.position_weighted_mae / n)
            .sum(),
        per_assay_metrics,
    })
}

fn analysis<'a>(analyses: &'a [AnalysisResult], family: &str) -> &'a AnalysisResult {
    analyses.iter().find(|a| a.family == family).unwrap()
}

fn position_differences(
    result: &AnalysisResult,
) -> Result<Vec<genome_research_stats::PositionLosses>> {
    let by_model: Vec<BTreeMap<_, _>> = result
        .models
        .iter()
        .map(|m| {
            m.per_position_errors
                .iter()
                .map(|p| {
                    (
                        (p.assay_id.clone(), p.position_1based),
                        p.mean_squared_error,
                    )
                })
                .collect()
        })
        .collect();
    let mut output = Vec::new();
    for key in by_model[0].keys() {
        ensure!(key.0 == HBB, "HBB bootstrap received a different assay");
        let differences = COMPARISONS
            .map(|(_, candidate, baseline)| by_model[baseline][key] - by_model[candidate][key]);
        output.push(genome_research_stats::PositionLosses {
            position_1based: key.1,
            loss_differences: differences,
        });
    }
    Ok(output)
}

fn element_differences(result: &AnalysisResult) -> Vec<[f64; 3]> {
    let by_model: Vec<BTreeMap<_, _>> = result
        .models
        .iter()
        .map(|m| {
            m.per_assay_metrics
                .iter()
                .map(|p| (p.assay_id.clone(), p.position_weighted_mse))
                .collect()
        })
        .collect();
    by_model[0]
        .keys()
        .map(|assay| {
            COMPARISONS.map(|(_, candidate, baseline)| {
                by_model[baseline][assay] - by_model[candidate][assay]
            })
        })
        .collect()
}

fn resample_losses(analyses: &[AnalysisResult]) -> Result<Value> {
    let positions = position_differences(analysis(analyses, "hbb_guarded_position_cv"))?;
    let elements = element_differences(analysis(analyses, "primary_leave_element_out"));
    ensure!(
        elements.len() == 9,
        "primary bootstrap requires nine frozen loci"
    );
    Ok(json!({
        "comparison_order":COMPARISONS.map(|(id,_,_)|id),
        "hbb_position_block":genome_research_stats::paired_block_bootstrap(&positions,21,2000,2_026_102_701)?,
        "primary_element":genome_research_stats::paired_element_bootstrap(&elements,2000,2_026_102_702)?,
        "relative_mse_definition":"Per-assay model MSE divided by that assay's training-mean-baseline MSE; null when denominator is zero",
        "interpretation":"Conditional descriptive resampling of frozen OOF losses; no refitting, no guarantee of population coverage, shared studies and overlapping OOF training fits remain"
    }))
}

fn control_mses(
    prepared: &[PreparedPlan],
    features: &FeatureRows,
    labels: &[f64],
    records: &[InputRecord],
) -> Result<[[f64; 6]; 2]> {
    let mut families = [[0.0; 6]; 2];
    for (family_index, family) in ["hbb_guarded_position_cv", "primary_leave_element_out"]
        .iter()
        .enumerate()
    {
        for model in 0..6 {
            let mut losses = BTreeMap::<(String, usize), (f64, usize)>::new();
            for prepared in prepared
                .iter()
                .filter(|p| p.plan.diagnostic.family == *family)
            {
                let outcomes: Vec<_> = prepared.plan.training.iter().map(|&i| labels[i]).collect();
                let fit = prepared.designs[model].fit(&outcomes)?;
                let fitted_values = predict_groups(prepared, &fit, model, features, records.len())?;
                for &index in &prepared.plan.test {
                    let prediction = fitted_values[index];
                    let entry = losses
                        .entry((
                            records[index].assay_id.clone(),
                            records[index].position_1based,
                        ))
                        .or_default();
                    entry.0 += (labels[index] - prediction).powi(2);
                    entry.1 += 1;
                }
            }
            let mut assay_losses = BTreeMap::<String, (f64, usize)>::new();
            for ((assay, _), (sum, count)) in losses {
                let entry = assay_losses.entry(assay).or_default();
                entry.0 += sum / count as f64;
                entry.1 += 1;
            }
            let n = assay_losses.len() as f64;
            families[family_index][model] = assay_losses
                .values()
                .map(|(sum, count)| sum / (*count as f64) / n)
                .sum();
        }
    }
    Ok(families)
}

fn spatial_control(
    prepared: &[PreparedPlan],
    features: &FeatureRows,
    labels: &[f64],
    records: &[InputRecord],
    protocol: &Protocol,
    analyses: &[AnalysisResult],
) -> Result<Value> {
    let assays: BTreeSet<_> = protocol
        .primary_ids()
        .into_iter()
        .chain([HBB.to_owned()])
        .collect();
    let mut strata = Vec::new();
    for assay in assays {
        for mutation in MUTATIONS {
            let indices: Vec<_> = records
                .iter()
                .enumerate()
                .filter(|(_, r)| {
                    r.assay_id == assay && format!("{}>{}", r.reference, r.alternate) == mutation
                })
                .map(|(i, _)| i)
                .collect();
            strata.push((assay.clone(), mutation, indices));
        }
    }
    let mut rng = genome_research_stats::SplitMix64::new(2_026_102_703);
    let mut statistics = Vec::new();
    let mut transformations = Vec::new();
    let mut family_statistics = [Vec::<[f64; 3]>::new(), Vec::<[f64; 3]>::new()];
    for replicate in 0..99 {
        let mut rotated = labels.to_vec();
        for (assay, mutation, indices) in &strata {
            let offset = if indices.len() >= 2 {
                Some(1 + rng.index(indices.len() - 1)?)
            } else {
                None
            };
            if let Some(offset) = offset {
                for (i, &destination) in indices.iter().enumerate() {
                    rotated[destination] = labels[indices[(i + offset) % indices.len()]];
                }
            }
            transformations.push(json!({"replicate_0based":replicate,"assay_id":assay,"mutation_identity":mutation,"stratum_rows":indices.len(),"circular_offset":offset}));
        }
        let mses = control_mses(prepared, features, &rotated, records)?;
        let hbb =
            COMPARISONS.map(|(_, candidate, baseline)| mses[0][baseline] - mses[0][candidate]);
        let primary =
            COMPARISONS.map(|(_, candidate, baseline)| mses[1][baseline] - mses[1][candidate]);
        family_statistics[0].push(hbb);
        family_statistics[1].push(primary);
        statistics.push(json!({"replicate_0based":replicate,"hbb_model_mses":mses[0],"primary_macro_model_mses":mses[1],
            "hbb_loss_improvements":hbb,"primary_loss_improvements":primary}));
    }
    let mut tails = [[0.0; 3]; 2];
    let mut observed = [[0.0; 3]; 2];
    for (family_index, name) in ["hbb_guarded_position_cv", "primary_leave_element_out"]
        .iter()
        .enumerate()
    {
        let result = analysis(analyses, name);
        observed[family_index] = COMPARISONS.map(|(_, candidate, baseline)| {
            result.models[baseline].macro_assay_position_mse
                - result.models[candidate].macro_assay_position_mse
        });
        for comparison in 0..3 {
            let nulls: Vec<_> = family_statistics[family_index]
                .iter()
                .map(|values| values[comparison])
                .collect();
            tails[family_index][comparison] = genome_research_stats::empirical_rotation_tail(
                observed[family_index][comparison],
                &nulls,
            )?;
        }
    }
    Ok(json!({"replicates":99,"seed":2_026_102_703u64,
        "transformation":"One SplitMix64 stream over replicate then fullURNlex and directed-mutation order; each destination receives source(i+nonzerooffset)%stratumsize; size<2 unchanged without RNG draw",
        "comparison_order":COMPARISONS.map(|(id,_,_)|id),"model_order":MODEL_NAMES,
        "observed_guarded_hbb_loss_improvements":observed[0],"observed_primary_loss_improvements":observed[1],
        "statistics":statistics,"transformations":transformations,
        "adjusted_tails":genome_research_stats::adjust_control_tails(tails[0],tails[1])?,
        "factorization_cache":"Training X means/scales and Gram Cholesky reused because X is fixed; every rotated training ymean and RHS is recomputed",
        "interpretation":"Artificial chemistry-stratified spatial rotation diagnostic; does not preserve cross-class/within-position score correlations, absolute gaps or all library dependence. Not a calibrated biological-significance p-value; no threshold-based discovery claim."
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn protocol() -> Protocol {
        Protocol::open(
            &Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../experiments/regulatory-prediction-v1.json"),
        )
        .unwrap()
    }
    fn motifs() -> MotifSet {
        MotifSet::open(
            &Path::new(env!("CARGO_MANIFEST_DIR")).join("../../datasets/regulatory-motifs"),
        )
        .unwrap()
    }
    fn toy_rows() -> Vec<InputRecord> {
        let protocol = protocol();
        let mut rows = Vec::new();
        for assay in protocol
            .primary_ids()
            .into_iter()
            .chain([HBB.into(), "urn:mavedb:00000023-a-2".into()])
        {
            let length = if assay == HBB { 187 } else { 35 };
            let sequence = "A".repeat(length);
            for position in 1..=length {
                for alternate in ["C", "G", "T"] {
                    rows.push(InputRecord {
                        id: format!("{assay}#n.{position}A>{alternate}"),
                        assay_id: assay.clone(),
                        position_1based: position,
                        reference: "A".into(),
                        alternate: alternate.into(),
                        target_sequence: sequence.clone(),
                        score: Some((position % 7) as f64 / 10.0),
                        score_raw: ((position % 7) as f64 / 10.0).to_string(),
                        unique_tags: Some(20.0),
                        dna_sequences: Some(100.0),
                        rna_sequences: Some(200.0),
                        cell_context: "synthetic test context".into(),
                    });
                }
            }
        }
        rows
    }
    #[test]
    fn blocked_folds_have_exact_protocol_counts_disjoint_supports_and_no_role_leakage() {
        let protocol = protocol();
        let (records, _) = prepare_inputs(&toy_rows(), &protocol).unwrap();
        let plans = plans(&records, &protocol).unwrap();
        let guarded: Vec<_> = plans
            .iter()
            .filter(|p| p.diagnostic.family == "hbb_guarded_position_cv")
            .collect();
        assert_eq!(
            guarded.iter().map(|p| p.test.len() / 3).collect::<Vec<_>>(),
            vec![37, 37, 38, 37, 38]
        );
        assert_eq!(
            guarded
                .iter()
                .map(|p| p.training.len() / 3)
                .collect::<Vec<_>>(),
            vec![130, 110, 109, 110, 129]
        );
        for plan in guarded {
            assert_eq!(
                plan.diagnostic
                    .minimum_same_assay_train_test_center_distance,
                Some(21)
            );
            assert_eq!(plan.diagnostic.shared_train_test_positions, 0);
            assert!(plan.diagnostic.strict_same_assay_feature_support_disjoint);
            for &train in &plan.training {
                for &test in &plan.test {
                    assert!(
                        records[train]
                            .position_1based
                            .abs_diff(records[test].position_1based)
                            > 20
                    );
                }
            }
        }
        for plan in plans
            .iter()
            .filter(|p| p.diagnostic.family == "primary_leave_element_out")
        {
            assert_eq!(plan.diagnostic.training_assays.len(), 8);
            assert!(!plan.diagnostic.training_assays.contains(&HBB.to_owned()));
            assert!(!plan
                .diagnostic
                .training_assays
                .contains(&"urn:mavedb:00000023-a-2".to_owned()));
            assert!(!plan
                .diagnostic
                .training_assays
                .contains(&plan.diagnostic.target_assay));
        }
        assert!(plans
            .iter()
            .filter(|p| p.diagnostic.family == "row_random_diagnostic")
            .all(|p| p.diagnostic.shared_train_test_positions > 0));
    }
    #[test]
    fn eligibility_keeps_valid_siblings_and_rejects_nonfinite_or_duplicate_inputs() {
        let protocol = protocol();
        let mut rows = toy_rows();
        rows[0].unique_tags = Some(9.0);
        let (included, report) = prepare_inputs(&rows, &protocol).unwrap();
        assert_eq!(included.len(), rows.len() - 1);
        assert!(included.iter().any(
            |r| r.assay_id == rows[0].assay_id && r.position_1based == rows[0].position_1based
        ));
        let item = report["per_assay"]
            .as_array()
            .unwrap()
            .iter()
            .find(|r| r["assay_id"] == rows[0].assay_id)
            .unwrap();
        assert_eq!(item["positions_with_two_alleles"], 1);
        let mut bad = rows.clone();
        bad[1].score = Some(f64::NAN);
        assert!(prepare_inputs(&bad, &protocol).is_err());
        let mut bad = rows.clone();
        bad.push(rows[1].clone());
        assert!(prepare_inputs(&bad, &protocol).is_err());
        let mut bad = rows;
        bad[1].reference = "C".into();
        assert!(prepare_inputs(&bad, &protocol).is_err());
    }
    #[test]
    fn fixed_sequence_features_ignore_outcomes_counts_and_cell_labels() {
        let motifs = motifs();
        let mut row = toy_rows().remove(0);
        let (original, _) = features(&row, &motifs).unwrap();
        assert_eq!(
            original.iter().map(Vec::len).collect::<Vec<_>>(),
            vec![0, 12, 14, 76, 16, 82]
        );
        assert_eq!(original[3][12], -1.0); // AAA count delta at the first target base.
        assert_eq!(original[3][12 + 16], 1.0); // CAA gain.
        row.score = Some(999.0);
        row.unique_tags = Some(9999.0);
        row.cell_context = "changed outcome metadata".into();
        let (changed, _) = features(&row, &motifs).unwrap();
        assert_eq!(original, changed);
        assert_eq!(feature_schema()["full"].len(), 82);
    }
    #[test]
    fn cached_design_refits_y_instead_of_reusing_real_outcome_means() {
        let design = RidgeDesign::prepare(&[vec![-1.0], vec![1.0]], 2.0).unwrap();
        let first = design.fit(&[1.0, 3.0]).unwrap();
        let second = design.fit(&[13.0, 11.0]).unwrap();
        assert_eq!(first.intercept, 2.0);
        assert_eq!(second.intercept, 12.0);
        assert!(
            first.standardized_coefficients[0] > 0.0 && second.standardized_coefficients[0] < 0.0
        );
        assert_eq!(first.training_feature_means, second.training_feature_means);
    }
    #[test]
    fn final_protocol_hash_rejects_silent_amendment() {
        let original = fs::read(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../experiments/regulatory-prediction-v1.json"),
        )
        .unwrap();
        Protocol::from_bytes(&original).unwrap();
        let mut altered = original;
        altered.push(b'\n');
        assert!(Protocol::from_bytes(&altered).is_err());
    }
    #[test]
    fn position_aggregation_preserves_constant_predictions_with_partial_alleles() {
        let records = toy_rows();
        let assay = records[0].assay_id.clone();
        let mut predictions = Vec::new();
        for position in 1..=3 {
            for allele in 0..position {
                predictions.push(Prediction {
                    id: format!("{position}/{allele}"),
                    assay_id: assay.clone(),
                    position_1based: position,
                    reference: "A".into(),
                    alternate: ["C", "G", "T"][allele].into(),
                    fold_index_0based: 0,
                    observed: position as f64 / 10.0,
                    observed_source_decimal: (position as f64 / 10.0).to_string(),
                    predicted: 0.1,
                    residual_observed_minus_predicted: position as f64 / 10.0 - 0.1,
                });
            }
        }
        let result =
            summarize_model("mean", Vec::new(), Vec::new(), predictions, &records).unwrap();
        assert!(result
            .per_position_errors
            .iter()
            .all(|p| p.predicted_mean == 0.1));
        assert!(result.per_assay_metrics[0]
            .position_mean_metrics
            .pearson
            .is_none());
        assert!(result.per_assay_metrics[0]
            .position_mean_metrics
            .spearman
            .is_none());
    }
}
