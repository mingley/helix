//! Post-hoc audit of evidence strength, ascertainment and split limitations.
//! This does not constitute new held-out evaluation or biological validation.
use anyhow::{ensure, Context, Result};
use genome_data::{sha256, Dataset, LabelRecord};
use genome_inference::{ConsequenceModel, Predictor};
use genome_types::{ClinicalLabel, Consequence, FeatureRecord};
use serde::Serialize;
use serde_json::Value;
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::Path,
};

#[derive(Debug, Serialize)]
pub struct AuditReport {
    pub schema_version: u32,
    pub dataset_manifest_sha256: String,
    pub source_sha256: String,
    pub analysis_status: String,
    pub source_pool: SourcePool,
    pub selected: PartitionSummary,
    pub training: PartitionSummary,
    pub heldout: PartitionSummary,
    pub partition_independence: PartitionIndependence,
    pub source_submissions: SubmissionSummary,
    pub model_gaps: Vec<ModelGap>,
    pub post_hoc_comparison: PostHocComparison,
    pub limitations: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct SourcePool {
    pub total_records: usize,
    pub classification_counts: BTreeMap<String, usize>,
    pub review_status_counts: BTreeMap<String, usize>,
    /// Stars summarize review process, never measurement independence.
    pub review_star_counts: BTreeMap<String, usize>,
    pub selected_records: usize,
    pub excluded_records: usize,
    pub exclusion_reason_counts: BTreeMap<String, usize>,
    pub exclusions: Vec<ExcludedRecord>,
}

#[derive(Debug, Serialize)]
pub struct ExcludedRecord {
    pub accession: String,
    pub title: String,
    pub classification: String,
    pub reason: String,
}

#[derive(Debug, Serialize)]
pub struct PartitionSummary {
    pub record_count: usize,
    pub unique_loci: usize,
    pub binary_class_counts: BTreeMap<String, usize>,
    pub original_classification_counts: BTreeMap<String, usize>,
    pub review_status_counts: BTreeMap<String, usize>,
    pub review_star_counts: BTreeMap<String, usize>,
    pub consequence_counts: BTreeMap<String, usize>,
    pub label_by_consequence: BTreeMap<String, BTreeMap<String, usize>>,
    /// A record may name several traits; these counts are not exclusive.
    pub trait_name_counts: BTreeMap<String, usize>,
    pub records_without_specific_condition: usize,
}

#[derive(Debug, Serialize)]
pub struct PartitionIndependence {
    pub exact_locus_overlap: usize,
    pub shared_gene: bool,
    pub shared_source_database: bool,
    pub independent_submitters_established: bool,
    pub independent_cohorts_established: bool,
    pub minimum_cross_partition_distance_bp: Option<u64>,
    pub cross_partition_pairs_within_5_bp: usize,
    pub cross_partition_pairs_within_50_bp: usize,
    pub nearest_cross_partition_pairs: Vec<NearbyPair>,
    pub interpretation: String,
}

#[derive(Debug, Serialize)]
pub struct NearbyPair {
    pub training_id: String,
    pub heldout_id: String,
    pub training_position_zero_based: u64,
    pub heldout_position_zero_based: u64,
    pub distance_bp: u64,
}

#[derive(Debug, Serialize)]
pub struct SubmissionSummary {
    pub unique_scv_accessions: usize,
    pub scv_accessions_repeated_across_records: usize,
    pub records_without_scv_accessions: usize,
    pub submitter_organization_identities_available: bool,
    pub experimental_cohort_identities_available: bool,
    pub interpretation: String,
}

#[derive(Debug, Serialize)]
pub struct ModelGap {
    pub consequence: String,
    pub benign_training_records: usize,
    pub pathogenic_training_records: usize,
    pub interpretation: String,
}

#[derive(Debug, Serialize)]
pub struct PostHocComparison {
    pub status: String,
    pub model: String,
    pub disagreement_count: usize,
    pub disagreements: Vec<ClassificationDisagreement>,
    pub interpretation: String,
}

#[derive(Debug, Serialize)]
pub struct ClassificationDisagreement {
    pub accession: String,
    pub title: String,
    pub source_url: String,
    pub variant_key: String,
    pub source_protein_change: String,
    pub clinical_classification: String,
    pub review_status: String,
    pub trait_names: Vec<String>,
    pub source_consequence: String,
    pub probability_pathogenic: f64,
    pub predicted_label: ClinicalLabel,
    pub observed_binary_label: ClinicalLabel,
}

/// Open only the existing validated, hash-pinned dataset. All labels are
/// deliberately visible in this descriptive audit: never use this as a fresh
/// model-selection benchmark after exploring it.
pub fn audit_dataset(directory: &Path) -> Result<AuditReport> {
    let dataset = Dataset::open(directory)?;
    let root = directory.canonicalize()?;
    let source_path = root.join("source.json").canonicalize()?;
    ensure!(
        source_path.starts_with(&root),
        "source escapes dataset directory"
    );
    let source_bytes = fs::read(source_path)?;
    ensure!(
        sha256(&source_bytes) == dataset.manifest.files["source.json"],
        "source changed after dataset validation"
    );
    let source: Value = serde_json::from_slice(&source_bytes)?;
    let result = source["result"]
        .as_object()
        .context("missing ESummary result")?;
    let uids = result
        .get("uids")
        .and_then(Value::as_array)
        .context("missing ESummary UID list")?;
    let mut seen_uids = BTreeSet::new();
    let mut by_accession = BTreeMap::new();
    let mut records = Vec::new();
    for uid in uids {
        let uid = uid.as_str().context("ESummary UID is not a string")?;
        ensure!(seen_uids.insert(uid), "duplicate ESummary UID");
        let record = result.get(uid).context("ESummary UID has no record")?;
        let accession = string(record, "accession_version")?;
        ensure!(
            by_accession.insert(accession.to_owned(), record).is_none(),
            "duplicate ESummary versioned accession"
        );
        records.push(record);
    }
    let labels: Vec<_> = dataset
        .training
        .iter()
        .map(|t| &t.observation)
        .chain(dataset.heldout())
        .collect();
    let selected_accessions: BTreeSet<_> = labels.iter().map(|l| l.id.as_str()).collect();
    let feature_map: BTreeMap<_, _> = dataset
        .features
        .iter()
        .map(|f| (f.id.as_str(), f))
        .collect();
    let mut classification_counts = BTreeMap::new();
    let mut review_status_counts = BTreeMap::new();
    let mut review_star_counts = BTreeMap::new();
    let mut exclusion_reason_counts = BTreeMap::new();
    let mut exclusions = Vec::new();
    let mut scv_counts = BTreeMap::new();
    let mut records_without_scv_accessions = 0;
    for record in &records {
        let accession = string(record, "accession_version")?;
        let classification = record["germline_classification"]["description"]
            .as_str()
            .context("missing germline classification")?;
        let review = record["germline_classification"]["review_status"]
            .as_str()
            .context("missing review status")?;
        increment(&mut classification_counts, classification);
        increment(&mut review_status_counts, review);
        increment(&mut review_star_counts, review_stars(review));
        if !selected_accessions.contains(accession) {
            let reason = exclusion_reason(record)?;
            increment(&mut exclusion_reason_counts, &reason);
            exclusions.push(ExcludedRecord {
                accession: accession.into(),
                title: string(record, "title")?.into(),
                classification: classification.into(),
                reason,
            });
        }
        let scvs = record["supporting_submissions"]["scv"]
            .as_array()
            .context("missing supporting submission accessions")?;
        if scvs.is_empty() {
            records_without_scv_accessions += 1;
        }
        let mut record_scvs = BTreeSet::new();
        for scv in scvs {
            let scv = scv.as_str().context("invalid SCV accession")?;
            if record_scvs.insert(scv) {
                increment(&mut scv_counts, scv);
            }
        }
    }
    ensure!(
        records.len() == labels.len() + exclusions.len(),
        "selected/excluded records do not exhaust source pool"
    );
    let training_labels: Vec<_> = dataset.training.iter().map(|t| &t.observation).collect();
    let heldout_labels: Vec<_> = dataset.heldout().iter().collect();
    let selected = partition(&labels, &feature_map, &by_accession)?;
    let training = partition(&training_labels, &feature_map, &by_accession)?;
    let heldout = partition(&heldout_labels, &feature_map, &by_accession)?;
    let partition_independence = split_audit(&training_labels, &heldout_labels, &feature_map);
    let model = ConsequenceModel::fit(&dataset.training)?;
    let mut model_gaps = Vec::new();
    for (group, counts) in &model.groups {
        if counts.pathogenic > 0 && counts.pathogenic < counts.total {
            model_gaps.push(ModelGap {
                consequence: consequence_name(*group),
                benign_training_records: (counts.total - counts.pathogenic) as usize,
                pathogenic_training_records: counts.pathogenic as usize,
                interpretation: "The same consequence group has different clinical assertions for different variants. This challenges group-level sufficiency; it is not a same-variant/context contradiction.".into(),
            });
        }
    }
    let mut disagreements = Vec::new();
    for label in dataset.heldout() {
        let feature = feature_map[label.id.as_str()];
        let prediction = model.predict(feature)?;
        if prediction.probability_pathogenic.label() != label.label {
            let record = by_accession[&label.id];
            disagreements.push(ClassificationDisagreement {
                accession: label.id.clone(),
                title: string(record, "title")?.into(),
                source_url: label.source_url.clone(),
                variant_key: feature.variant.key(),
                source_protein_change: string(record, "protein_change")?.into(),
                clinical_classification: label.classification.clone(),
                review_status: label.review_status.clone(),
                trait_names: trait_names(record)?,
                source_consequence: consequence_name(feature.consequence),
                probability_pathogenic: prediction.probability_pathogenic.value(),
                predicted_label: prediction.probability_pathogenic.label(),
                observed_binary_label: label.label,
            });
        }
    }
    Ok(AuditReport {
        schema_version: 1,
        dataset_manifest_sha256: dataset.manifest_sha256.clone(),
        source_sha256: dataset.manifest.files["source.json"].clone(),
        analysis_status: "descriptive post-hoc exploration of an already exposed demonstration dataset".into(),
        source_pool: SourcePool {
            total_records: records.len(), classification_counts, review_status_counts,
            review_star_counts, selected_records: labels.len(), excluded_records: exclusions.len(),
            exclusion_reason_counts, exclusions,
        },
        selected, training, heldout, partition_independence,
        source_submissions: SubmissionSummary {
            unique_scv_accessions: scv_counts.len(),
            scv_accessions_repeated_across_records: scv_counts.values().filter(|&&n| n > 1).count(),
            records_without_scv_accessions,
            submitter_organization_identities_available: false,
            experimental_cohort_identities_available: false,
            interpretation: "SCV accessions identify submissions, not necessarily distinct organizations, studies, patients, or experiments. Multiple submitters/no conflicts is a review category, not proof of replication.".into(),
        },
        model_gaps,
        post_hoc_comparison: PostHocComparison {
            status: "post_hoc_on_previously_exposed_holdout".into(),
            model: "consequence-beta-binomial-v1".into(),
            disagreement_count: disagreements.len(), disagreements,
            interpretation: "Disagreements prioritize source/assay investigation. Neither a classifier score nor a curated clinical label establishes measured molecular function. The explored holdout cannot now be treated as untouched validation for future tuning.".into(),
        },
        limitations: vec![
            "The first 100 search hits are an ascertainment-selected recent-record sample, not a census of HBB variation.".into(),
            "Likely and definitive classifications were collapsed; uncertainty and conflicts were excluded rather than resolved.".into(),
            "Aggregate clinical labels can mix disease, inheritance, zygosity, population and assay contexts.".into(),
            "Review stars summarize assertion review; they do not establish biological truth or independent experiments.".into(),
            "A source snapshot can contain no conflicts because of selection; absence here does not establish consensus for the gene.".into(),
            "Disjoint loci do not establish independent cohorts, submitters, haplotypes, regulatory contexts or studies.".into(),
            "The source contains no measured effect sizes, assay units, individual participant identifiers or experimentally established cohort independence.".into(),
            "Sequence-derived annotation agreement can validate computation but cannot independently validate pathogenicity.".into(),
        ],
    })
}

fn increment(map: &mut BTreeMap<String, usize>, value: &str) {
    *map.entry(value.to_owned()).or_default() += 1;
}
fn string<'a>(record: &'a Value, key: &str) -> Result<&'a str> {
    record[key]
        .as_str()
        .with_context(|| format!("missing string field {key}"))
}
fn review_stars(status: &str) -> &str {
    match status {
        "practice guideline" => "4",
        "reviewed by expert panel" => "3",
        "criteria provided, multiple submitters, no conflicts" => "2",
        "criteria provided, single submitter"
        | "criteria provided, conflicting classifications" => "1",
        "no assertion criteria provided"
        | "no assertion provided"
        | "no classification provided" => "0",
        _ => "unmapped",
    }
}
fn consequence_name(group: Consequence) -> String {
    match group {
        Consequence::Missense => "missense",
        Consequence::Synonymous => "synonymous",
        Consequence::Splice => "splice",
        Consequence::Stop => "stop",
        Consequence::Intron => "intron",
        Consequence::Utr => "utr",
        Consequence::Other => "other",
    }
    .into()
}
fn exclusion_reason(record: &Value) -> Result<String> {
    let classification = record["germline_classification"]["description"]
        .as_str()
        .context("missing excluded classification")?;
    if classification == "Uncertain significance" {
        return Ok("uncertain_significance".into());
    }
    if classification.to_ascii_lowercase().contains("conflict") {
        return Ok("conflicting_classifications".into());
    }
    if !matches!(
        classification,
        "Benign" | "Likely benign" | "Pathogenic" | "Likely pathogenic"
    ) {
        return Ok("other_nonbinary_classification".into());
    }
    let variations = record["variation_set"]
        .as_array()
        .context("missing variation set")?;
    if variations.len() != 1 {
        return Ok("multiple_variant_or_haplotype_record".into());
    }
    let spdi = variations[0]["canonical_spdi"].as_str().unwrap_or_default();
    let fields: Vec<_> = spdi.split(':').collect();
    if fields.len() != 4 {
        return Ok("missing_or_invalid_spdi".into());
    }
    if fields[0] != "NC_000011.10" {
        return Ok("different_reference_sequence".into());
    }
    if fields[2].len() != 1
        || fields[3].len() != 1
        || !fields[2].bytes().all(|b| b"ACGT".contains(&b))
        || !fields[3].bytes().all(|b| b"ACGT".contains(&b))
        || fields[2] == fields[3]
    {
        return Ok("unsupported_or_identical_alleles".into());
    }
    Ok("eligible_record_not_selected_or_duplicate".into())
}
fn trait_names(record: &Value) -> Result<Vec<String>> {
    record["germline_classification"]["trait_set"]
        .as_array()
        .context("missing trait set")?
        .iter()
        .map(|t| Ok(string(t, "trait_name")?.to_owned()))
        .collect()
}
fn partition(
    labels: &[&LabelRecord],
    features: &BTreeMap<&str, &FeatureRecord>,
    source: &BTreeMap<String, &Value>,
) -> Result<PartitionSummary> {
    let mut summary = PartitionSummary {
        record_count: labels.len(),
        unique_loci: 0,
        binary_class_counts: BTreeMap::new(),
        original_classification_counts: BTreeMap::new(),
        review_status_counts: BTreeMap::new(),
        review_star_counts: BTreeMap::new(),
        consequence_counts: BTreeMap::new(),
        label_by_consequence: BTreeMap::new(),
        trait_name_counts: BTreeMap::new(),
        records_without_specific_condition: 0,
    };
    let mut loci = BTreeSet::new();
    for label in labels {
        let feature = features
            .get(label.id.as_str())
            .context("missing audited feature")?;
        let source = source.get(&label.id).context("missing source accession")?;
        let binary = if label.label.positive() {
            "pathogenic"
        } else {
            "benign"
        };
        let consequence = consequence_name(feature.consequence);
        increment(&mut summary.binary_class_counts, binary);
        increment(
            &mut summary.original_classification_counts,
            &label.classification,
        );
        increment(&mut summary.review_status_counts, &label.review_status);
        increment(
            &mut summary.review_star_counts,
            review_stars(&label.review_status),
        );
        increment(&mut summary.consequence_counts, &consequence);
        increment(
            summary.label_by_consequence.entry(consequence).or_default(),
            binary,
        );
        loci.insert(feature.variant.position().offset());
        let traits = trait_names(source)?;
        if traits.is_empty()
            || traits.iter().all(|name| {
                matches!(
                    name.to_ascii_lowercase().as_str(),
                    "not provided" | "not specified"
                )
            })
        {
            summary.records_without_specific_condition += 1;
        }
        let unique_traits: BTreeSet<_> = traits.into_iter().collect();
        for name in unique_traits {
            increment(&mut summary.trait_name_counts, &name);
        }
    }
    summary.unique_loci = loci.len();
    Ok(summary)
}
fn split_audit(
    training: &[&LabelRecord],
    heldout: &[&LabelRecord],
    features: &BTreeMap<&str, &FeatureRecord>,
) -> PartitionIndependence {
    let mut pairs = Vec::new();
    let training_loci: BTreeSet<_> = training
        .iter()
        .map(|l| features[l.id.as_str()].variant.position().offset())
        .collect();
    let heldout_loci: BTreeSet<_> = heldout
        .iter()
        .map(|l| features[l.id.as_str()].variant.position().offset())
        .collect();
    for train in training {
        for test in heldout {
            let a = features[train.id.as_str()].variant.position().offset();
            let b = features[test.id.as_str()].variant.position().offset();
            pairs.push(NearbyPair {
                training_id: train.id.clone(),
                heldout_id: test.id.clone(),
                training_position_zero_based: a,
                heldout_position_zero_based: b,
                distance_bp: a.abs_diff(b),
            });
        }
    }
    pairs.sort_by(|a, b| {
        (a.distance_bp, &a.training_id, &a.heldout_id).cmp(&(
            b.distance_bp,
            &b.training_id,
            &b.heldout_id,
        ))
    });
    let minimum = pairs.first().map(|p| p.distance_bp);
    let within5 = pairs.iter().filter(|p| p.distance_bp <= 5).count();
    let within50 = pairs.iter().filter(|p| p.distance_bp <= 50).count();
    pairs.truncate(10);
    PartitionIndependence {
        exact_locus_overlap: training_loci.intersection(&heldout_loci).count(),
        shared_gene: true, shared_source_database: true,
        independent_submitters_established: false, independent_cohorts_established: false,
        minimum_cross_partition_distance_bp: minimum,
        cross_partition_pairs_within_5_bp: within5, cross_partition_pairs_within_50_bp: within50,
        nearest_cross_partition_pairs: pairs,
        interpretation: "The position-hash split prevents exact-locus overlap but leaves nearby variants, one gene and one curation source shared. Proximity is a leakage risk to investigate, not proof that individual observations share a cohort.".into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    fn fixture() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../datasets/clinvar-hbb")
    }
    #[test]
    fn source_selection_and_review_strength_are_auditable() {
        let report = audit_dataset(&fixture()).unwrap();
        assert_eq!(report.source_pool.total_records, 100);
        assert_eq!(report.source_pool.selected_records, 54);
        assert_eq!(report.source_pool.excluded_records, 46);
        assert_eq!(
            report.source_pool.exclusion_reason_counts["uncertain_significance"],
            45
        );
        assert_eq!(
            report.source_pool.exclusion_reason_counts["multiple_variant_or_haplotype_record"],
            1
        );
        assert_eq!(report.heldout.review_star_counts["0"], 2);
        assert_eq!(report.heldout.review_star_counts["1"], 23);
        assert!(!report.heldout.review_star_counts.contains_key("3"));
        assert!(
            !report
                .partition_independence
                .independent_cohorts_established
        );
        assert_eq!(report.partition_independence.exact_locus_overlap, 0);
        assert_eq!(
            report
                .partition_independence
                .minimum_cross_partition_distance_bp,
            Some(3)
        );
        assert_eq!(
            report
                .partition_independence
                .cross_partition_pairs_within_5_bp,
            2
        );
        assert_eq!(
            report
                .partition_independence
                .cross_partition_pairs_within_50_bp,
            39
        );
    }
    #[test]
    fn a_model_disagreement_is_explicitly_post_hoc_and_not_a_measured_contradiction() {
        let report = audit_dataset(&fixture()).unwrap();
        let comparison = report.post_hoc_comparison;
        assert_eq!(comparison.status, "post_hoc_on_previously_exposed_holdout");
        assert_eq!(comparison.disagreement_count, 1);
        let example = &comparison.disagreements[0];
        assert_eq!(example.accession, "VCV004685642.1");
        assert_eq!(example.source_protein_change, "H78L");
        assert_eq!(example.clinical_classification, "Likely benign");
        assert_eq!(example.trait_names, ["not provided"]);
        assert_eq!(example.predicted_label, ClinicalLabel::Pathogenic);
        assert_eq!(example.observed_binary_label, ClinicalLabel::Benign);
    }
    #[test]
    fn altered_source_cannot_enter_audit() {
        let temp = tempfile::tempdir().unwrap();
        for entry in fs::read_dir(fixture()).unwrap() {
            let entry = entry.unwrap();
            if entry.file_type().unwrap().is_file() {
                fs::copy(entry.path(), temp.path().join(entry.file_name())).unwrap();
            }
        }
        fs::write(temp.path().join("source.json"), "{}").unwrap();
        assert!(audit_dataset(temp.path())
            .unwrap_err()
            .to_string()
            .contains("checksum mismatch"));
    }
    #[test]
    fn unknown_review_status_does_not_default_to_zero_stars() {
        assert_eq!(review_stars("an unfamiliar future status"), "unmapped");
        assert_eq!(review_stars("reviewed by expert panel"), "3");
    }
}
