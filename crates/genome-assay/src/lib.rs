//! Reference-validated measured promoter reporter effects, separate from clinical labels.
//!
//! MaveDB's n.-coordinates are positions in its submitted 187-base target, not
//! NM_000518.5 cDNA coordinates. Source equality records are retained separately:
//! they are not inferred deletions or validated neutral wild-type controls.
use anyhow::{bail, ensure, Context, Result};
use genome_mechanism::{ReferenceBundle, Region};
use genome_types::{GRCh38, GenomicPosition, Variant};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Component, Path},
};

const ASSAY_URN: &str = "urn:mavedb:00000018-a-1";
const HEADERS: [&str; 9] = [
    "accession",
    "hgvs_nt",
    "hgvs_splice",
    "hgvs_pro",
    "score",
    "p-value",
    "unique_tags",
    "dna_sequences",
    "rna_sequences",
];

#[derive(Deserialize)]
struct Manifest {
    schema_version: u32,
    assay_urn: String,
    assembly: String,
    chromosome: String,
    target_start: u64,
    target_end: u64,
    target_strand: String,
    gene_transcription_strand: String,
    expected_rows: usize,
    cell_line: String,
    score_name: String,
    score_units: String,
    replicates: Value,
    reference_fasta_sha256: String,
    source_target_sequence_sha256: String,
    data_license: Value,
    sources: Vec<Value>,
    files: BTreeMap<String, String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum NumericStatus {
    Finite,
    Missing,
    NonNumeric,
    NonFinite,
    OutsideDomain,
}

#[derive(Clone, Debug, Serialize)]
pub struct SourceNumber {
    pub raw: String,
    pub value: Option<f64>,
    pub status: NumericStatus,
}

#[derive(Clone, Copy)]
enum NumericDomain {
    AnyFinite,
    Probability,
    Count,
}

#[derive(Clone, Copy, Debug, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ReporterUnit {
    Log2ReporterExpressionCoefficient,
}

#[derive(Clone, Debug, Serialize)]
pub struct ReporterEffect {
    pub coefficient_log2: f64,
    pub unit: ReporterUnit,
    pub assay_urn: String,
    pub cell_line: String,
    pub experimental_context: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Representation {
    SingleNucleotideSubstitution,
    SourceEquality,
    UnsupportedDeletion,
    UnsupportedInsertionOrDuplication,
    UnsupportedComplexAllele,
    UnsupportedCoordinateSystem,
    UnsupportedOrMalformedDescription,
}

#[derive(Clone, Debug, Serialize)]
pub struct AssaySourceRow {
    /// One-based data row number, excluding the CSV header.
    pub source_row_1based: usize,
    pub source_accession: String,
    pub source_hgvs_nt: String,
    pub source_hgvs_splice: String,
    pub source_hgvs_pro: String,
    pub representation: Representation,
    pub target_position_1based: Option<usize>,
    pub exclusion_reason_from_snv_analysis: Option<String>,
    pub score: SourceNumber,
    pub p_value: SourceNumber,
    pub unique_tags: SourceNumber,
    pub dna_sequences: SourceNumber,
    pub rna_sequences: SourceNumber,
}

#[derive(Clone, Debug, Serialize)]
pub struct MeasuredSnvRecord {
    pub variant_key: String,
    pub variant: Variant<GRCh38>,
    pub source_accession: String,
    pub source_hgvs_nt: String,
    pub target_position_1based: usize,
    pub genomic_position_0based: u64,
    pub region: Region,
    /// None means the SNV maps but its reported score is not a finite number.
    /// Missing scores never remove a source row from the representation counts.
    pub effect: Option<ReporterEffect>,
    pub score: SourceNumber,
    pub p_value: SourceNumber,
    pub unique_tags: SourceNumber,
    pub dna_sequences: SourceNumber,
    pub rna_sequences: SourceNumber,
}

pub struct PromoterAssay {
    manifest_sha256: String,
    records: Vec<MeasuredSnvRecord>,
    source_rows: Vec<AssaySourceRow>,
    report: Value,
}

impl PromoterAssay {
    pub fn open(directory: &Path, reference: &ReferenceBundle) -> Result<Self> {
        let root = directory
            .canonicalize()
            .context("promoter assay directory")?;
        let manifest_bytes = fs::read(root.join("manifest.json"))?;
        let manifest: Manifest = serde_json::from_slice(&manifest_bytes)?;
        ensure!(
            manifest.schema_version == 1
                && manifest.assay_urn == ASSAY_URN
                && manifest.assembly == "GRCh38"
                && manifest.chromosome == "11",
            "unsupported promoter assay schema, accession, assembly, or chromosome"
        );
        ensure!(
            manifest.target_start < manifest.target_end
                && manifest.target_end - manifest.target_start == 187
                && manifest.expected_rows > 0
                && manifest.expected_rows <= 100_000,
            "invalid promoter target bounds or row budget"
        );
        ensure!(
            manifest.gene_transcription_strand == reference.manifest().transcript.strand,
            "gene strand differs from reference"
        );
        ensure!(
            manifest.reference_fasta_sha256 == reference.manifest().files["genomic.fasta"],
            "assay was mapped against a different frozen genomic reference"
        );
        for name in [
            "metadata.json",
            "scores.csv",
            "target.fasta",
            "license-scope.json",
        ] {
            ensure!(
                manifest.files.contains_key(name),
                "assay manifest does not pin {name}"
            );
        }
        let mut contents = BTreeMap::new();
        for (name, expected) in &manifest.files {
            let mut parts = Path::new(name).components();
            ensure!(
                matches!(parts.next(), Some(Component::Normal(_))) && parts.next().is_none(),
                "assay files must be basenames"
            );
            let path = root.join(name).canonicalize()?;
            ensure!(path.starts_with(&root), "assay file escapes bundle");
            let bytes = fs::read(path)?;
            ensure!(
                digest(&bytes) == *expected,
                "assay checksum mismatch: {name}"
            );
            contents.insert(name.as_str(), bytes);
        }
        ensure!(
            !manifest.sources.is_empty(),
            "assay source provenance is required"
        );
        for source in &manifest.sources {
            ensure!(
                source["url"]
                    .as_str()
                    .is_some_and(|url| url.starts_with("https://")),
                "assay source requires an HTTPS URL"
            );
            let file = source["file"]
                .as_str()
                .context("assay source has no file")?;
            ensure!(
                manifest.files.contains_key(file),
                "assay provenance names an unpinned file"
            );
            if let Some(hash) = source["sha256"].as_str() {
                ensure!(
                    manifest.files.get(file).map(String::as_str) == Some(hash),
                    "assay source commitment mismatch"
                );
            }
        }
        let metadata: Value = serde_json::from_slice(&contents["metadata.json"])?;
        let license_scope: Value = serde_json::from_slice(&contents["license-scope.json"])?;
        validate_metadata(&manifest, &metadata, &license_scope)?;
        let target = target_fasta(&contents["target.fasta"])?;
        ensure!(target.len() == 187, "assay target FASTA length mismatch");
        let source_target = metadata["targetGenes"][0]["targetSequence"]["sequence"]
            .as_str()
            .context("missing source target sequence")?;
        ensure!(
            source_target.as_bytes() == target,
            "target FASTA differs from MaveDB target"
        );
        ensure!(
            digest(&target) == manifest.source_target_sequence_sha256,
            "source target sequence commitment mismatch"
        );
        let matches = matching_intervals(
            reference.genomic_sequence().as_bytes(),
            &target,
            reference.manifest().sequence_start,
        )?;
        ensure!(
            matches.len() == 1,
            "assay target does not map uniquely in the bundled reference"
        );
        let mapping = &matches[0];
        ensure!(
            mapping.start == manifest.target_start
                && mapping.end == manifest.target_end
                && mapping.strand == manifest.target_strand,
            "assay target coordinates/orientation disagree with exact reference match"
        );
        // MaveDB's genomic coordinates are verified here as one-based closed,
        // not assumed from the style of a chromosome/start/end field.
        let source_coordinates = &metadata["extraMetadata"];
        ensure!(
            source_coordinates["reference"] == "hg38"
                && source_coordinates["chr"] == "11"
                && source_coordinates["start"].as_u64() == Some(mapping.start + 1)
                && source_coordinates["end"].as_u64() == Some(mapping.end),
            "source genomic coordinate convention does not reconcile with reference"
        );
        let mut reader = csv::ReaderBuilder::new().from_reader(contents["scores.csv"].as_slice());
        ensure!(
            reader.headers()?.iter().eq(HEADERS),
            "unsupported promoter assay CSV header"
        );
        let mut source_rows = Vec::new();
        let mut records = Vec::new();
        let mut accessions = BTreeSet::new();
        let mut representations = BTreeSet::new();
        let mut variant_keys = BTreeSet::new();
        for (index, raw) in reader.records().enumerate() {
            ensure!(
                index < manifest.expected_rows,
                "assay has more rows than declared"
            );
            let raw = raw.context("invalid promoter assay CSV row")?;
            let accession = raw.get(0).context("missing assay accession")?;
            let hgvs = raw.get(1).context("missing target-relative HGVS")?;
            ensure!(
                accession.starts_with(&format!("{ASSAY_URN}#"))
                    && accession.len() > ASSAY_URN.len() + 1,
                "invalid assay row accession"
            );
            ensure!(
                accessions.insert(accession.to_owned()),
                "duplicate assay accession"
            );
            ensure!(
                representations.insert(hgvs.to_owned()),
                "duplicate assay source representation"
            );
            let parsed = parse_hgvs(hgvs, &target)?;
            let score = source_number(&raw[4], NumericDomain::AnyFinite);
            let p_value = source_number(&raw[5], NumericDomain::Probability);
            let unique_tags = source_number(&raw[6], NumericDomain::Count);
            let dna_sequences = source_number(&raw[7], NumericDomain::Count);
            let rna_sequences = source_number(&raw[8], NumericDomain::Count);
            if let Some((position, original, alternate)) = parsed.snv {
                let variant = map_target_snv(mapping, position, original, alternate)?;
                let annotation = reference.annotate(&variant)?;
                ensure!(
                    variant_keys.insert(variant.key()),
                    "duplicate normalized assay SNV"
                );
                let effect = score.value.map(|coefficient_log2| ReporterEffect {
                    coefficient_log2,
                    unit: ReporterUnit::Log2ReporterExpressionCoefficient,
                    assay_urn: manifest.assay_urn.clone(),
                    cell_line: manifest.cell_line.clone(),
                    experimental_context: "Published plasmid MPRA; linear-model coefficient derived from reporter DNA/RNA barcode counts in random-mutagenesis libraries. Constructs can contain multiple changes and the model assumes additive effects; not a direct endogenous HBB or clinical outcome".into(),
                });
                records.push(MeasuredSnvRecord {
                    variant_key: variant.key(),
                    genomic_position_0based: variant.position().offset(),
                    variant,
                    source_accession: accession.into(),
                    source_hgvs_nt: hgvs.into(),
                    target_position_1based: position,
                    region: annotation.region,
                    effect,
                    score: score.clone(),
                    p_value: p_value.clone(),
                    unique_tags: unique_tags.clone(),
                    dna_sequences: dna_sequences.clone(),
                    rna_sequences: rna_sequences.clone(),
                });
            }
            source_rows.push(AssaySourceRow {
                source_row_1based: index + 1,
                source_accession: accession.into(),
                source_hgvs_nt: hgvs.into(),
                source_hgvs_splice: raw[2].into(),
                source_hgvs_pro: raw[3].into(),
                representation: parsed.representation,
                target_position_1based: parsed.position,
                exclusion_reason_from_snv_analysis: parsed.exclusion_reason,
                score,
                p_value,
                unique_tags,
                dna_sequences,
                rna_sequences,
            });
        }
        ensure!(
            source_rows.len() == manifest.expected_rows,
            "assay row count differs from source metadata"
        );
        let mut representation_counts = BTreeMap::<String, usize>::new();
        let mut score_status_counts = BTreeMap::<String, usize>::new();
        let mut counts_status = BTreeMap::<String, BTreeMap<String, usize>>::new();
        for row in &source_rows {
            increment(&mut representation_counts, enum_name(&row.representation)?);
            increment(&mut score_status_counts, enum_name(&row.score.status)?);
            for (name, number) in [
                ("p_value", &row.p_value),
                ("unique_tags", &row.unique_tags),
                ("dna_sequences", &row.dna_sequences),
                ("rna_sequences", &row.rna_sequences),
            ] {
                increment(
                    counts_status.entry(name.into()).or_default(),
                    enum_name(&number.status)?,
                );
            }
        }
        let positions: BTreeSet<_> = records.iter().map(|r| r.target_position_1based).collect();
        let mut by_region = BTreeMap::<String, Vec<f64>>::new();
        let mut region_record_counts = BTreeMap::<String, usize>::new();
        for record in &records {
            let region = enum_name(&record.region)?;
            increment(&mut region_record_counts, region.clone());
            if let Some(effect) = &record.effect {
                by_region
                    .entry(region)
                    .or_default()
                    .push(effect.coefficient_log2);
            }
        }
        let region_scores: BTreeMap<_, _> = by_region
            .into_iter()
            .map(|(region, values)| (region, summarize(&values)))
            .collect();
        let all_scores: Vec<_> = source_rows.iter().filter_map(|r| r.score.value).collect();
        let snv_scores: Vec<_> = records
            .iter()
            .filter_map(|r| r.effect.as_ref().map(|effect| effect.coefficient_log2))
            .collect();
        let equality_scores: Vec<_> = source_rows
            .iter()
            .filter(|r| r.representation == Representation::SourceEquality)
            .filter_map(|r| r.score.value)
            .collect();
        let mut target_regions = BTreeMap::<String, usize>::new();
        for position in mapping.start..mapping.end {
            let index = usize::try_from(position - reference.manifest().sequence_start)?;
            let base = reference.genomic_sequence().as_bytes()[index] as char;
            let alternate = if base == 'A' { 'C' } else { 'A' };
            let variant = Variant::<GRCh38>::snv(
                GenomicPosition::new("11", position)?,
                &base.to_string(),
                &alternate.to_string(),
            )?;
            let region = reference.annotate(&variant)?.region;
            increment(&mut target_regions, enum_name(&region)?);
        }
        let report = json!({
            "schema_version":1,
            "assay_urn":manifest.assay_urn,
            "source_metadata_sha256":manifest.files["metadata.json"],
            "source_scores_sha256":manifest.files["scores.csv"],
            "reference_manifest_sha256":reference.manifest_sha256(),
            "data_license":manifest.data_license,
            "data_license_scope":license_scope,
            "endpoint":{
                "name":"published reporter-expression linear-model coefficient",
                "unit":"log2_reporter_expression_coefficient",
                "source_unit_description":manifest.score_units,
                "cell_line":manifest.cell_line,
                "experimental_context":"plasmid MPRA",
                "replicate_context":manifest.replicates,
                "method_limits":[
                    "Random-mutagenesis constructs can carry multiple changes; reported per-variant effects are additive linear-model coefficients rather than isolated-variant endogenous measurements",
                    "The primary study reports inflated coefficient p-values; supplied p-values are retained without assigning new significance labels",
                    "Three cultures share one constructed library; culture/barcode counts do not establish independent study replication"
                ],
                "interpretation":"Measured reporter endpoint derived by the source study's regression model; not endogenous HBB expression, protein function, a pathogenicity probability, or a calibrated clinical score"
            },
            "mapping":{
                "assembly":manifest.assembly,"chromosome":manifest.chromosome,
                "target_interval_zero_based_half_open":[mapping.start,mapping.end],
                "target_reference_alignment_strand":mapping.strand,
                "hbb_gene_transcription_strand":manifest.gene_transcription_strand,
                "target_sequence_bases":target.len(),
                "exact_sequence_matches_in_bundled_interval":matches.len(),
                "genome_wide_uniqueness_established":false,
                "hgvs_coordinate_system":"MaveDB n. indexes its own submitted target, one-based; these are not RefSeq cDNA or chromosome coordinates",
                "source_coordinate_convention":"one-based inclusive start/end, reconciled by full-target sequence equality",
                "target_region_base_counts":target_regions,
                "orientation_limit":"Target-to-reference alignment orientation does not by itself establish physical reporter construct orientation"
            },
            "rows":{
                "all_source_rows":source_rows.len(),"expected_source_rows":manifest.expected_rows,
                "representation_counts":representation_counts,
                "reference_validated_snv_rows":records.len(),
                "excluded_from_snv_analysis":source_rows.len()-records.len(),
                "source_equality_rows":source_rows.iter().filter(|row| row.representation==Representation::SourceEquality).count(),
                "numeric_finite_score_rows_all_source":all_scores.len(),
                "numeric_finite_score_rows_snvs":snv_scores.len(),
                "score_status_counts":score_status_counts,
                "provided_numeric_field_status_counts":counts_status,
                "source_equality_interpretation":"Literal n.<position>= descriptions are retained exactly and excluded from SNV summaries; not reinterpreted as deletions or validated neutral wild-type controls"
            },
            "coverage":{
                "unique_target_positions_with_snvs":positions.len(),
                "target_positions":target.len(),
                "potential_single_base_substitutions":target.len()*3,
                "observed_unique_reference_validated_snvs":variant_keys.len(),
                "complete_target_snv_enumeration":variant_keys.len()==target.len()*3,
                "snv_record_counts_by_reference_region":region_record_counts
            },
            "descriptive_scores":{
                "all_source_rows_including_equality":summarize(&all_scores),
                "reference_validated_snvs":summarize(&snv_scores),
                "source_equality_records_not_a_neutral_control":summarize(&equality_scores),
                "snvs_by_reference_region":region_scores,
                "quantile_method":"Linear interpolation at rank p*(n-1), equivalent to R type 7; descriptive only"
            },
            "analysis_status":"Descriptive external functional-assay integration; no model fitting, threshold tuning, association testing, or new untouched-holdout claim",
            "limitations":[
                "All source rows are preserved and scores keep their supplied units; unsupported representations and missing/nonfinite scores are separately accounted for",
                "Scores from variants at the same position share local sequence and assay context; unique tags, cultures and independent studies are different units",
                "Reporter coefficients in HEL 92.1.7 cells do not establish endogenous expression, clinical pathogenicity or universal cell-context effects",
                "The assay is an experimental endpoint beyond annotation/clinical curation; database and clinical evidence ancestry may overlap",
                "Several rows were seen in availability checks; this is not an untouched or independently preregistered validation set",
                "No transcription-factor motif, causal pathway, effect significance threshold or clinical label is inferred from these descriptive summaries"
            ]
        });
        Ok(Self {
            manifest_sha256: digest(&manifest_bytes),
            records,
            source_rows,
            report,
        })
    }
    pub fn report(&self) -> &Value {
        &self.report
    }
    pub fn manifest_sha256(&self) -> &str {
        &self.manifest_sha256
    }
    pub fn records(&self) -> &[MeasuredSnvRecord] {
        &self.records
    }
    pub fn source_rows(&self) -> &[AssaySourceRow] {
        &self.source_rows
    }
}

fn validate_metadata(manifest: &Manifest, metadata: &Value, license_scope: &Value) -> Result<()> {
    ensure!(
        metadata["urn"].as_str() == Some(manifest.assay_urn.as_str()),
        "score-set accession differs from metadata"
    );
    ensure!(
        metadata["numVariants"].as_u64() == Some(manifest.expected_rows as u64),
        "row count differs from MaveDB metadata"
    );
    ensure!(
        metadata["license"] == manifest.data_license
            && license_scope["data_license"] == manifest.data_license,
        "assay data license commitments disagree"
    );
    ensure!(
        manifest.data_license["shortName"] == "CC0"
            && manifest.data_license["version"] == "1.0"
            && manifest.data_license["link"]
                == "https://creativecommons.org/publicdomain/zero/1.0/",
        "this adapter requires the verified CC0 source grant"
    );
    ensure!(
        license_scope["scope_source_doi"] == "10.1186/s13059-019-1845-6"
            && license_scope["scope_excerpt"]
                .as_str()
                .is_some_and(|s| s.contains("downloaded file")),
        "data-license scope evidence is missing"
    );
    ensure!(
        manifest.cell_line == "HEL 92.1.7"
            && manifest.score_name == "score"
            && manifest.score_units == "log2 reporter-expression effect; linear-model coefficient",
        "unsupported assay endpoint/context declaration"
    );
    let targets = metadata["targetGenes"]
        .as_array()
        .context("no assay targets")?;
    ensure!(
        targets.len() == 1
            && targets[0]["category"] == "regulatory"
            && targets[0]["targetSequence"]["sequenceType"] == "dna"
            && targets[0]["targetSequence"]["taxonomy"]["code"] == 9606,
        "unsupported target or organism"
    );
    let method = metadata["methodText"]
        .as_str()
        .context("assay method description missing")?;
    ensure!(
        method.contains("Log2 variant expression effect")
            && method.contains("linear model")
            && method.contains("unique barcode tags"),
        "source does not describe the expected endpoint"
    );
    ensure!(
        metadata["experiment"]["shortDescription"]
            .as_str()
            .is_some_and(|s| s.contains("HEL 92.1.7")),
        "source cell context differs from assay declaration"
    );
    Ok(())
}

#[derive(Clone, Debug)]
struct TargetMapping {
    start: u64,
    end: u64,
    strand: String,
}

fn matching_intervals(
    reference: &[u8],
    target: &[u8],
    reference_start: u64,
) -> Result<Vec<TargetMapping>> {
    ensure!(!target.is_empty(), "empty assay target");
    let reverse: Vec<_> = target
        .iter()
        .rev()
        .map(|&base| complement(base))
        .collect::<Result<_>>()?;
    let mut matches = Vec::new();
    if target.len() > reference.len() {
        return Ok(matches);
    }
    for (offset, window) in reference.windows(target.len()).enumerate() {
        for (sequence, strand) in [(target, "+"), (reverse.as_slice(), "-")] {
            if window == sequence {
                matches.push(TargetMapping {
                    start: reference_start + offset as u64,
                    end: reference_start + offset as u64 + target.len() as u64,
                    strand: strand.into(),
                });
            }
        }
    }
    Ok(matches)
}

struct ParsedHgvs {
    representation: Representation,
    position: Option<usize>,
    snv: Option<(usize, u8, u8)>,
    exclusion_reason: Option<String>,
}

fn parse_hgvs(description: &str, target: &[u8]) -> Result<ParsedHgvs> {
    let excluded = |representation, reason: &str, position| ParsedHgvs {
        representation,
        position,
        snv: None,
        exclusion_reason: Some(reason.into()),
    };
    let Some(body) = description.strip_prefix("n.") else {
        return Ok(excluded(Representation::UnsupportedCoordinateSystem, "Only positive n. coordinates relative to the source target are supported; no implicit cDNA/genomic reinterpretation", None));
    };
    let digits = body.bytes().take_while(u8::is_ascii_digit).count();
    if digits == 0 {
        let category = if body.contains(';') || body.contains('[') {
            Representation::UnsupportedComplexAllele
        } else {
            Representation::UnsupportedOrMalformedDescription
        };
        return Ok(excluded(
            category,
            "Compound, offset, negative or malformed target coordinates are not mapped as SNVs",
            None,
        ));
    }
    let position: usize = body[..digits]
        .parse()
        .context("target-relative coordinate overflow")?;
    ensure!(
        (1..=target.len()).contains(&position),
        "assay target coordinate outside reference sequence"
    );
    if digits > 1 && body.starts_with('0') {
        return Ok(excluded(
            Representation::UnsupportedOrMalformedDescription,
            "Noncanonical zero-padded HGVS coordinate is retained without normalization",
            Some(position),
        ));
    }
    let suffix = &body[digits..];
    if suffix == "=" {
        return Ok(excluded(Representation::SourceEquality, "Source equality description; not an SNV, inferred deletion, or validated neutral control", Some(position)));
    }
    let bytes = suffix.as_bytes();
    if bytes.len() == 3
        && bytes[1] == b'>'
        && b"ACGT".contains(&bytes[0])
        && b"ACGT".contains(&bytes[2])
    {
        ensure!(
            target[position - 1] == bytes[0],
            "assay HGVS reference allele disagrees with target sequence"
        );
        ensure!(
            bytes[0] != bytes[2],
            "substitution notation has equal alleles; equality is not silently normalized"
        );
        return Ok(ParsedHgvs {
            representation: Representation::SingleNucleotideSubstitution,
            position: Some(position),
            snv: Some((position, bytes[0], bytes[2])),
            exclusion_reason: None,
        });
    }
    let (representation, reason) = if suffix.contains("del") && !suffix.contains("ins") {
        (
            Representation::UnsupportedDeletion,
            "Deletion representation is retained but not interpreted as an SNV",
        )
    } else if suffix.contains("dup") || (suffix.contains("ins") && !suffix.contains("del")) {
        (
            Representation::UnsupportedInsertionOrDuplication,
            "Insertion/duplication representation is retained but not interpreted as an SNV",
        )
    } else if suffix.contains("delins") || suffix.contains(';') || suffix.contains('_') {
        (
            Representation::UnsupportedComplexAllele,
            "Complex or multi-position allele is retained but not interpreted as an SNV",
        )
    } else {
        (Representation::UnsupportedOrMalformedDescription, "Unsupported or malformed allele description is retained without genomic reinterpretation")
    };
    Ok(excluded(representation, reason, Some(position)))
}

fn map_target_snv(
    mapping: &TargetMapping,
    position: usize,
    reference: u8,
    alternate: u8,
) -> Result<Variant<GRCh38>> {
    ensure!(
        position > 0 && position as u64 <= mapping.end - mapping.start,
        "target position outside mapped interval"
    );
    let (offset, reference, alternate) = match mapping.strand.as_str() {
        "+" => (mapping.start + position as u64 - 1, reference, alternate),
        "-" => (
            mapping.end - position as u64,
            complement(reference)?,
            complement(alternate)?,
        ),
        _ => bail!("unsupported target alignment strand"),
    };
    Variant::snv(
        GenomicPosition::new("11", offset)?,
        &(reference as char).to_string(),
        &(alternate as char).to_string(),
    )
}

fn source_number(raw: &str, domain: NumericDomain) -> SourceNumber {
    let normalized = raw.trim();
    let missing = normalized.is_empty()
        || ["NA", "N/A", "NULL", "."]
            .iter()
            .any(|token| normalized.eq_ignore_ascii_case(token));
    let (value, status) = if missing {
        (None, NumericStatus::Missing)
    } else {
        match normalized.parse::<f64>() {
            Err(_) => (None, NumericStatus::NonNumeric),
            Ok(value) if !value.is_finite() => (None, NumericStatus::NonFinite),
            Ok(value) => {
                let in_domain = match domain {
                    NumericDomain::AnyFinite => true,
                    NumericDomain::Probability => (0.0..=1.0).contains(&value),
                    NumericDomain::Count => {
                        value >= 0.0 && value.fract() == 0.0 && value <= (1u64 << 53) as f64
                    }
                };
                if in_domain {
                    (Some(value), NumericStatus::Finite)
                } else {
                    (None, NumericStatus::OutsideDomain)
                }
            }
        }
    };
    SourceNumber {
        raw: raw.into(),
        value,
        status,
    }
}

fn summarize(values: &[f64]) -> Option<Value> {
    if values.is_empty() {
        return None;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let quantile = |p: f64| {
        let rank = p * (sorted.len() - 1) as f64;
        let low = rank.floor() as usize;
        let high = rank.ceil() as usize;
        sorted[low] * (1.0 - (rank - low as f64)) + sorted[high] * (rank - low as f64)
    };
    Some(json!({
        "finite_count":sorted.len(),"minimum":sorted[0],"q25":quantile(0.25),
        "median":quantile(0.5),"q75":quantile(0.75),"maximum":sorted[sorted.len()-1],
        "mean":sorted.iter().map(|value| value / sorted.len() as f64).sum::<f64>()
    }))
}

fn target_fasta(bytes: &[u8]) -> Result<Vec<u8>> {
    let text = std::str::from_utf8(bytes)?;
    let mut lines = text.lines();
    ensure!(
        lines.next().is_some_and(|line| line.starts_with('>')),
        "assay FASTA header missing"
    );
    let mut sequence = Vec::new();
    for line in lines {
        ensure!(
            !line.starts_with('>'),
            "multiple assay target FASTA records"
        );
        sequence.extend_from_slice(line.trim().as_bytes());
    }
    ensure!(
        !sequence.is_empty() && sequence.iter().all(|base| b"ACGT".contains(base)),
        "assay target must contain unambiguous uppercase DNA"
    );
    Ok(sequence)
}

fn complement(base: u8) -> Result<u8> {
    match base {
        b'A' => Ok(b'T'),
        b'T' => Ok(b'A'),
        b'C' => Ok(b'G'),
        b'G' => Ok(b'C'),
        _ => bail!("ambiguous target DNA"),
    }
}
fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn enum_name(value: &impl Serialize) -> Result<String> {
    serde_json::to_value(value)?
        .as_str()
        .map(str::to_owned)
        .context("expected string enum")
}
fn increment(map: &mut BTreeMap<String, usize>, key: String) {
    *map.entry(key).or_default() += 1;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../datasets/hbb-promoter-assay")
    }
    fn reference() -> ReferenceBundle {
        ReferenceBundle::open(
            &Path::new(env!("CARGO_MANIFEST_DIR")).join("../../datasets/hbb-reference"),
        )
        .unwrap()
    }
    fn copied_fixture() -> tempfile::TempDir {
        let temp = tempfile::tempdir().unwrap();
        for entry in fs::read_dir(fixture()).unwrap() {
            let entry = entry.unwrap();
            if entry.file_type().unwrap().is_file() {
                fs::copy(entry.path(), temp.path().join(entry.file_name())).unwrap();
            }
        }
        temp
    }
    fn recommit(root: &Path, name: &str, bytes: &[u8]) {
        fs::write(root.join(name), bytes).unwrap();
        let path = root.join("manifest.json");
        let mut manifest: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        let hash = digest(bytes);
        manifest["files"][name] = json!(hash);
        for source in manifest["sources"].as_array_mut().unwrap() {
            if source["file"] == name && source.get("sha256").is_some() {
                source["sha256"] = json!(hash);
            }
        }
        fs::write(path, serde_json::to_vec(&manifest).unwrap()).unwrap();
    }
    fn mutate_rows(root: &Path, mutate: impl FnOnce(&mut Vec<Vec<String>>)) {
        let mut reader = csv::Reader::from_path(root.join("scores.csv")).unwrap();
        let mut rows: Vec<Vec<String>> = reader
            .records()
            .map(|record| record.unwrap().iter().map(str::to_owned).collect())
            .collect();
        mutate(&mut rows);
        let mut writer = csv::Writer::from_writer(Vec::new());
        writer.write_record(HEADERS).unwrap();
        for row in rows {
            writer.write_record(row).unwrap();
        }
        recommit(root, "scores.csv", &writer.into_inner().unwrap());
    }
    fn load_error(root: &Path) -> String {
        PromoterAssay::open(root, &reference())
            .err()
            .unwrap()
            .to_string()
    }

    #[test]
    fn frozen_assay_reconciles_all_source_rows_and_measured_snvs() {
        let assay = PromoterAssay::open(&fixture(), &reference()).unwrap();
        assert_eq!(assay.records().len(), 561);
        assert_eq!(assay.source_rows().len(), 623);
        assert_eq!(assay.report()["rows"]["source_equality_rows"], 62);
        assert_eq!(
            assay.report()["rows"]["numeric_finite_score_rows_all_source"],
            623
        );
        assert_eq!(
            assay.report()["mapping"]["target_reference_alignment_strand"],
            "+"
        );
        assert_eq!(
            assay.report()["mapping"]["hbb_gene_transcription_strand"],
            "-"
        );
        assert_eq!(
            assay.report()["mapping"]["target_region_base_counts"]["five_prime_utr"],
            50
        );
        assert_eq!(
            assay.report()["mapping"]["target_region_base_counts"]["upstream"],
            137
        );
        assert_eq!(
            assay.report()["coverage"]["snv_record_counts_by_reference_region"]["five_prime_utr"],
            150
        );
        assert_eq!(
            assay.report()["coverage"]["snv_record_counts_by_reference_region"]["upstream"],
            411
        );
        assert_eq!(
            assay.report()["coverage"]["complete_target_snv_enumeration"],
            true
        );
        assert!(assay.records().iter().all(|record| record
            .effect
            .as_ref()
            .unwrap()
            .coefficient_log2
            .is_finite()));
        let specific = assay
            .records()
            .iter()
            .find(|record| record.source_hgvs_nt == "n.43G>A")
            .unwrap();
        assert_eq!(specific.genomic_position_0based, 5_227_063);
        assert_eq!(specific.region, Region::FivePrimeUtr);
        assert_eq!(specific.variant.reference(), 'G');
        assert_eq!(specific.variant.alternate(), 'A');
    }

    #[test]
    fn target_snv_mapping_preserves_reference_orientation_not_gene_strand() {
        let forward = TargetMapping {
            start: 100,
            end: 103,
            strand: "+".into(),
        };
        let reverse = TargetMapping {
            start: 100,
            end: 103,
            strand: "-".into(),
        };
        let a = map_target_snv(&forward, 1, b'A', b'C').unwrap();
        let b = map_target_snv(&reverse, 1, b'A', b'C').unwrap();
        assert_eq!(a.key(), "GRCh38:11:100:A>C");
        assert_eq!(b.key(), "GRCh38:11:102:T>G");
        let matches = matching_intervals(b"AAACGCCC", b"ACG", 100).unwrap();
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].start, 102);
        assert_eq!(matches[0].strand, "+");
        assert_eq!(matching_intervals(b"ACGT", b"ACGT", 0).unwrap().len(), 2); // palindromic orientation is ambiguous
    }

    #[test]
    fn strict_target_relative_hgvs_does_not_reinterpret_equality_or_complex_alleles() {
        let parsed = parse_hgvs("n.2C>T", b"ACGT").unwrap();
        assert_eq!(parsed.snv, Some((2, b'C', b'T')));
        let equality = parse_hgvs("n.2=", b"ACGT").unwrap();
        assert_eq!(equality.representation, Representation::SourceEquality);
        assert!(equality.snv.is_none());
        for (description, expected) in [
            ("n.2del", Representation::UnsupportedDeletion),
            ("n.1_2delinsTT", Representation::UnsupportedComplexAllele),
            ("n.[1A>T;2C>G]", Representation::UnsupportedComplexAllele),
            ("c.2C>T", Representation::UnsupportedCoordinateSystem),
            ("n.02C>T", Representation::UnsupportedOrMalformedDescription),
        ] {
            let parsed = parse_hgvs(description, b"ACGT").unwrap();
            assert_eq!(parsed.representation, expected);
            assert!(parsed.snv.is_none());
        }
        assert!(parse_hgvs("n.0=", b"ACGT").is_err());
        assert!(parse_hgvs("n.5A>T", b"ACGT").is_err());
        assert!(parse_hgvs("n.2A>T", b"ACGT").is_err());
        assert!(parse_hgvs("n.2C>C", b"ACGT").is_err());
    }

    #[test]
    fn missingness_and_numeric_domains_remain_distinct() {
        assert_eq!(
            source_number("NA", NumericDomain::AnyFinite).status,
            NumericStatus::Missing
        );
        assert_eq!(
            source_number("nonsense", NumericDomain::AnyFinite).status,
            NumericStatus::NonNumeric
        );
        assert_eq!(
            source_number("NaN", NumericDomain::AnyFinite).status,
            NumericStatus::NonFinite
        );
        assert_eq!(
            source_number("inf", NumericDomain::AnyFinite).status,
            NumericStatus::NonFinite
        );
        assert_eq!(
            source_number("68.0", NumericDomain::Count).value,
            Some(68.0)
        );
        assert_eq!(
            source_number("68.5", NumericDomain::Count).status,
            NumericStatus::OutsideDomain
        );
        assert_eq!(
            source_number("-1", NumericDomain::Count).status,
            NumericStatus::OutsideDomain
        );
        assert_eq!(
            source_number("1.2", NumericDomain::Probability).status,
            NumericStatus::OutsideDomain
        );
        assert_eq!(
            source_number("-0.51", NumericDomain::AnyFinite).value,
            Some(-0.51)
        );
    }

    #[test]
    fn missing_scores_do_not_discard_rows_or_miscount_equality() {
        let temp = copied_fixture();
        mutate_rows(temp.path(), |rows| {
            rows[0][4] = "NaN".into(); // equality row
            rows[1][4] = "NA".into(); // SNV row
            rows[2][4] = "not-a-number".into(); // another SNV row
        });
        let assay = PromoterAssay::open(temp.path(), &reference()).unwrap();
        assert_eq!(assay.records().len(), 561);
        assert_eq!(assay.source_rows().len(), 623);
        assert_eq!(assay.report()["rows"]["source_equality_rows"], 62);
        assert_eq!(
            assay.report()["rows"]["numeric_finite_score_rows_all_source"],
            620
        );
        assert_eq!(
            assay.report()["rows"]["numeric_finite_score_rows_snvs"],
            559
        );
        assert_eq!(assay.report()["rows"]["score_status_counts"]["missing"], 1);
        assert_eq!(
            assay.report()["rows"]["score_status_counts"]["non_numeric"],
            1
        );
        assert_eq!(
            assay.report()["rows"]["score_status_counts"]["non_finite"],
            1
        );
    }

    #[test]
    fn unsupported_alleles_are_retained_with_explicit_exclusion() {
        let temp = copied_fixture();
        mutate_rows(temp.path(), |rows| {
            rows[1][1] = "n.3del".into();
        });
        let assay = PromoterAssay::open(temp.path(), &reference()).unwrap();
        assert_eq!(assay.source_rows().len(), 623);
        assert_eq!(assay.records().len(), 560);
        assert_eq!(
            assay.report()["rows"]["representation_counts"]["unsupported_deletion"],
            1
        );
        assert!(assay.source_rows()[1]
            .exclusion_reason_from_snv_analysis
            .is_some());
    }

    #[test]
    fn checksums_reference_mismatches_and_duplicate_accessions_fail() {
        let temp = copied_fixture();
        fs::write(temp.path().join("scores.csv"), "tampered").unwrap();
        assert!(load_error(temp.path()).contains("checksum mismatch"));
        let temp = copied_fixture();
        mutate_rows(temp.path(), |rows| {
            rows[1][1] = "n.3A>G".into();
        });
        assert!(load_error(temp.path()).contains("reference allele disagrees"));
        let temp = copied_fixture();
        mutate_rows(temp.path(), |rows| {
            rows[1][0] = rows[0][0].clone();
        });
        assert!(load_error(temp.path()).contains("duplicate assay accession"));
    }

    #[test]
    fn source_orientation_licenses_and_sequence_identity_are_checked() {
        let temp = copied_fixture();
        let path = temp.path().join("manifest.json");
        let mut manifest: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        manifest["target_strand"] = json!("-");
        fs::write(&path, serde_json::to_vec(&manifest).unwrap()).unwrap();
        assert!(load_error(temp.path()).contains("coordinates/orientation"));
        let temp = copied_fixture();
        let mut metadata: Value =
            serde_json::from_slice(&fs::read(temp.path().join("metadata.json")).unwrap()).unwrap();
        metadata["license"]["shortName"] = json!("unknown");
        recommit(
            temp.path(),
            "metadata.json",
            &serde_json::to_vec(&metadata).unwrap(),
        );
        assert!(load_error(temp.path()).contains("data license"));
        let temp = copied_fixture();
        let mut metadata: Value =
            serde_json::from_slice(&fs::read(temp.path().join("metadata.json")).unwrap()).unwrap();
        let target = metadata["targetGenes"][0]["targetSequence"]["sequence"]
            .as_str()
            .unwrap()
            .to_owned();
        metadata["targetGenes"][0]["targetSequence"]["sequence"] =
            json!(format!("A{}", &target[1..]));
        recommit(
            temp.path(),
            "metadata.json",
            &serde_json::to_vec(&metadata).unwrap(),
        );
        assert!(load_error(temp.path()).contains("FASTA differs"));
    }

    #[test]
    fn descriptive_quantiles_have_declared_interpolation_and_empty_behavior() {
        let summary = summarize(&[-2.0, 0.0, 2.0, 4.0]).unwrap();
        assert_eq!(summary["q25"], -0.5);
        assert_eq!(summary["median"], 1.0);
        assert_eq!(summary["q75"], 2.5);
        assert_eq!(summary["mean"], 1.0);
        assert!(summarize(&[]).is_none());
    }
}
