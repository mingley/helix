//! Licensed, metadata-selected regulatory assay ingestion. Source outcomes keep
//! reporter units and remain separate from clinical claims or independent studies.
use anyhow::{bail, ensure, Context, Result};
use genome_types::{GRCh38, GenomicPosition, Variant};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Component, Path},
};
const UNIT: &str = "log2_reporter_expression_coefficient";
const DOI: &str = "10.1038/s41467-019-11526-w";
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

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Mapping {
    pub assembly: String,
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub target_reference_alignment_strand: String,
    pub source_coordinate_convention: String,
    pub genome_wide_uniqueness_established: bool,
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RegulatoryAssay {
    pub assay_id: String,
    pub element_id: String,
    pub role: String,
    pub title: String,
    pub cell_context: String,
    pub assay_context: String,
    pub target_sequence: String,
    pub target_length: usize,
    pub metadata_file: String,
    pub reference_file: String,
    pub scores_file: String,
    pub expected_source_rows: usize,
    pub mapping: Mapping,
}
#[derive(Clone, Debug, Serialize)]
pub struct SourceNumber {
    pub raw: String,
    pub value: Option<f64>,
    pub status: String,
}
#[derive(Clone, Debug, Serialize)]
pub struct RegulatoryRecord {
    pub assay_id: String,
    pub element_id: String,
    pub id: String,
    pub role: String,
    pub target_position_1based: usize,
    pub reference: String,
    pub alternate: String,
    pub target_sequence: String,
    pub score: Option<f64>,
    pub score_status: String,
    pub score_raw: String,
    pub unique_tags: Option<f64>,
    pub dna_sequences: Option<f64>,
    pub rna_sequences: Option<f64>,
    pub cell_context: String,
    pub assay_context: String,
    pub source_hgvs_nt: String,
    pub variant_key: String,
    pub genomic_position_0based: u64,
    pub genomic_reference: String,
    pub genomic_alternate: String,
}
#[derive(Clone, Debug, Serialize)]
pub struct RegulatorySourceRow {
    pub assay_id: String,
    pub id: String,
    pub source_row_1based: usize,
    pub source_hgvs_nt: String,
    pub representation: String,
    pub score: SourceNumber,
    pub p_value: SourceNumber,
    pub unique_tags: SourceNumber,
    pub dna_sequences: SourceNumber,
    pub rna_sequences: SourceNumber,
    pub exclusion_reason_from_snv_representation: Option<String>,
}
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Manifest {
    pub schema_version: u32,
    pub cohort_id: String,
    pub study_doi: String,
    pub assay_count: usize,
    pub primary_element_count: usize,
    pub score_unit: String,
    pub endpoint: String,
    pub independent_study_count: usize,
    pub metadata_manifest_sha256: String,
    pub frozen_protocol_sha256: String,
    pub outcome_acquisition_gate_sha256: String,
    pub acquired_at: String,
    pub assays: Vec<RegulatoryAssay>,
    pub row_accounting: Vec<Value>,
    pub sources: Vec<Value>,
    pub files: BTreeMap<String, String>,
    pub limitations: Vec<String>,
}
pub struct RegulatoryCohort {
    pub manifest: Manifest,
    manifest_sha256: String,
    records: Vec<RegulatoryRecord>,
    source_rows: Vec<RegulatorySourceRow>,
    report: Value,
}
impl RegulatoryCohort {
    pub fn open(directory: &Path) -> Result<Self> {
        let root = directory
            .canonicalize()
            .context("regulatory cohort directory")?;
        let manifest_bytes = fs::read(root.join("manifest.json"))?;
        let manifest: Manifest = serde_json::from_slice(&manifest_bytes)?;
        ensure!(
            manifest.schema_version == 1
                && manifest.cohort_id == "kircher2019-human-promoters-v1"
                && manifest.study_doi == DOI
                && manifest.score_unit == UNIT
                && manifest.independent_study_count == 1,
            "unsupported cohort identity/endpoint"
        );
        ensure!(
            manifest.assay_count == 15
                && manifest.assay_count == manifest.assays.len()
                && manifest.primary_element_count == 9,
            "unexpected frozen cohort size"
        );
        let mut contents = BTreeMap::new();
        for (name, expected) in &manifest.files {
            let mut components = Path::new(name).components();
            ensure!(
                matches!(components.next(), Some(Component::Normal(_)))
                    && components.next().is_none(),
                "cohort paths must be basenames"
            );
            let path = root.join(name).canonicalize()?;
            ensure!(path.starts_with(&root), "cohort source escapes bundle");
            let bytes = fs::read(path)?;
            ensure!(
                digest(&bytes) == *expected,
                "cohort checksum mismatch: {name}"
            );
            contents.insert(name.as_str(), bytes);
        }
        let get = |name: &str| -> Result<&[u8]> {
            contents
                .get(name)
                .map(Vec::as_slice)
                .context("cohort manifest is missing a consumed input")
        };
        ensure!(
            digest(get("metadata-manifest.json")?) == manifest.metadata_manifest_sha256
                && digest(get("frozen-protocol.json")?) == manifest.frozen_protocol_sha256
                && digest(get("outcome-acquisition-gate.json")?)
                    == manifest.outcome_acquisition_gate_sha256,
            "cohort phase/protocol commitments disagree"
        );
        let metadata_phase: Value = serde_json::from_slice(get("metadata-manifest.json")?)?;
        ensure!(
            metadata_phase["phase"] == "metadata_only_no_score_endpoints_requested",
            "missing metadata-first acquisition phase"
        );
        let phase_files = metadata_phase["files"]
            .as_object()
            .context("metadata phase commitments missing")?;
        for (name, hash) in phase_files {
            ensure!(
                hash.as_str() == Some(digest(get(name)?).as_str()),
                "metadata changed after phase freeze"
            );
        }
        let gate: Value = serde_json::from_slice(get("outcome-acquisition-gate.json")?)?;
        let protocol: Value = serde_json::from_slice(get("frozen-protocol.json")?)?;
        ensure!(
            gate["frozen_before_new_outcome_acquisition"] == true
                && gate["protocol_sha256"] == manifest.frozen_protocol_sha256
                && gate["metadata_manifest_sha256"] == manifest.metadata_manifest_sha256,
            "acquisition gate does not bind protocol/metadata phase"
        );
        ensure!(
            protocol["status"] == "frozen"
                && protocol["frozen_before_new_outcome_acquisition"] == true,
            "protocol was not frozen before new outcome acquisition"
        );
        let declared = protocol["cohort"]["assay_roles"]
            .as_array()
            .context("frozen protocol has no role inventory")?;
        let mut expected_roles = BTreeMap::new();
        for item in declared {
            let id = item["assay_id"]
                .as_str()
                .context("invalid assay role identity")?;
            let role = item["role"].as_str().context("invalid assay role")?;
            ensure!(
                expected_roles.insert(id, role).is_none(),
                "duplicate protocol assay ID"
            );
        }
        ensure!(
            expected_roles.len() == 15,
            "protocol cohort role count differs"
        );
        let actual_roles: BTreeMap<_, _> = manifest
            .assays
            .iter()
            .map(|a| (a.assay_id.as_str(), a.role.as_str()))
            .collect();
        ensure!(
            actual_roles == expected_roles && actual_roles.len() == 15,
            "cohort roles differ from frozen protocol"
        );
        ensure!(
            gate["assay_roles"] == serde_json::to_value(&actual_roles)?,
            "gate roles differ from protocol"
        );
        let mut seen_sources = BTreeSet::new();
        for source in &manifest.sources {
            let name = source["file"]
                .as_str()
                .context("source has no retained file")?;
            ensure!(
                source["url"]
                    .as_str()
                    .is_some_and(|u| u.starts_with("https://"))
                    && manifest.files.contains_key(name),
                "invalid source provenance URL/file"
            );
            ensure!(
                source["sha256"].as_str() == manifest.files.get(name).map(String::as_str),
                "source provenance commitment differs"
            );
            seen_sources.insert(name);
        }
        let license: Value = serde_json::from_slice(get("license-scope.json")?)?;
        ensure!(
            license["scope_source_doi"] == "10.1186/s13059-019-1845-6"
                && license["scope_excerpt"]
                    .as_str()
                    .is_some_and(|s| s.contains("downloaded file")),
            "per-score-set data license scope missing"
        );
        ensure!(
            manifest
                .assays
                .iter()
                .filter(|a| a.role == "validation_only")
                .count()
                == 5
                && manifest
                    .assays
                    .iter()
                    .filter(|a| a.role == "previously_exposed_reference")
                    .count()
                    == 1
                && manifest
                    .assays
                    .iter()
                    .any(|a| a.assay_id == "urn:mavedb:00000018-a-1"
                        && a.role == "previously_exposed_reference"),
            "validation/reference role counts differ from the fixed cohort"
        );
        let mut records = Vec::new();
        let mut source_rows = Vec::new();
        let mut accounting = Vec::new();
        let mut seen_assays = BTreeSet::new();
        for assay in &manifest.assays {
            ensure!(
                seen_assays.insert(&assay.assay_id)
                    && matches!(
                        assay.role.as_str(),
                        "primary" | "validation_only" | "previously_exposed_reference"
                    ),
                "duplicate assay or invalid role"
            );
            ensure!(
                assay.target_length == assay.target_sequence.len()
                    && (150..=1000).contains(&assay.target_length)
                    && assay.target_sequence.bytes().all(|b| b"ACGT".contains(&b)),
                "invalid cohort target DNA"
            );
            ensure!(
                assay.expected_source_rows > 0 && assay.expected_source_rows <= 100_000,
                "invalid source row budget"
            );
            for name in [
                &assay.metadata_file,
                &assay.reference_file,
                &assay.scores_file,
            ] {
                ensure!(
                    seen_sources.contains(name.as_str()),
                    "consumed source has no provenance record"
                );
            }
            let metadata: Value = serde_json::from_slice(get(&assay.metadata_file)?)?;
            validate_metadata(assay, &metadata, &license)?;
            let reference: Value = serde_json::from_slice(get(&assay.reference_file)?)?;
            validate_mapping(assay, &reference, &metadata)?;
            let mut reader = csv::Reader::from_reader(get(&assay.scores_file)?);
            ensure!(
                reader.headers()?.iter().eq(HEADERS),
                "unsupported regulatory CSV headers"
            );
            let mut ids = BTreeSet::new();
            let mut descriptions = BTreeSet::new();
            let mut keys = BTreeSet::new();
            let mut positions = BTreeSet::new();
            let mut representations = BTreeMap::<String, usize>::new();
            let mut statuses = BTreeMap::<String, usize>::new();
            let mut count = 0usize;
            for (index, row) in reader.records().enumerate() {
                count += 1;
                ensure!(
                    count <= assay.expected_source_rows,
                    "extra source rows beyond declared metadata"
                );
                let row = row?;
                let id = &row[0];
                let hgvs = &row[1];
                ensure!(
                    id.starts_with(&format!("{}#", assay.assay_id))
                        && ids.insert(id.to_owned())
                        && descriptions.insert(hgvs.to_owned()),
                    "duplicate or invalid source identity/description"
                );
                let (representation, snv) = parse_hgvs(hgvs, assay.target_sequence.as_bytes())?;
                let score = number(&row[4], "score");
                let p = number(&row[5], "probability");
                let tags = number(&row[6], "count");
                let dna = number(&row[7], "count");
                let rna = number(&row[8], "count");
                *representations
                    .entry(representation.to_owned())
                    .or_default() += 1;
                *statuses.entry(score.status.clone()).or_default() += 1;
                if let Some((position, original, alternate)) = snv {
                    let variant = map_snv(&assay.mapping, position, original, alternate)?;
                    ensure!(keys.insert(variant.key()), "duplicate normalized SNV");
                    positions.insert(position);
                    records.push(RegulatoryRecord {
                        assay_id: assay.assay_id.clone(),
                        element_id: assay.element_id.clone(),
                        id: id.into(),
                        role: assay.role.clone(),
                        target_position_1based: position,
                        reference: (original as char).to_string(),
                        alternate: (alternate as char).to_string(),
                        target_sequence: assay.target_sequence.clone(),
                        score: score.value,
                        score_status: score.status.clone(),
                        score_raw: score.raw.clone(),
                        unique_tags: tags.value,
                        dna_sequences: dna.value,
                        rna_sequences: rna.value,
                        cell_context: assay.cell_context.clone(),
                        assay_context: assay.assay_context.clone(),
                        source_hgvs_nt: hgvs.into(),
                        variant_key: variant.key(),
                        genomic_position_0based: variant.position().offset(),
                        genomic_reference: variant.reference().to_string(),
                        genomic_alternate: variant.alternate().to_string(),
                    });
                }
                source_rows.push(RegulatorySourceRow{assay_id:assay.assay_id.clone(),id:id.into(),source_row_1based:index+1,source_hgvs_nt:hgvs.into(),representation:representation.into(),score,p_value:p,unique_tags:tags,dna_sequences:dna,rna_sequences:rna,exclusion_reason_from_snv_representation:snv.is_none().then(||"Source representation is retained but not treated as an unambiguous SNV or neutral control".into())});
            }
            ensure!(
                count == assay.expected_source_rows,
                "source rows differ from metadata count"
            );
            accounting.push(json!({"assay_id":assay.assay_id,"element_id":assay.element_id,"role":assay.role,"all_source_rows":count,"representation_counts":representations,"score_status_counts":statuses,"reference_validated_snv_rows":keys.len(),"excluded_representation_rows":count-keys.len(),"unique_target_positions_with_snvs":positions.len()}));
        }
        ensure!(
            accounting == manifest.row_accounting,
            "normalized row accounting differs from frozen importer summary"
        );
        let primary: Vec<_> = manifest
            .assays
            .iter()
            .filter(|a| a.role == "primary")
            .collect();
        ensure!(
            primary.len() == 9,
            "primary cohort must contain nine elements"
        );
        for (i, a) in primary.iter().enumerate() {
            for b in &primary[i + 1..] {
                ensure!(
                    a.mapping.chromosome != b.mapping.chromosome
                        || a.mapping.start.max(b.mapping.start) >= a.mapping.end.min(b.mapping.end),
                    "overlapping primary loci"
                );
                ensure!(
                    a.target_sequence != b.target_sequence
                        && a.target_sequence != reverse(b.target_sequence.as_bytes())?,
                    "duplicate primary target sequence"
                );
            }
        }
        let report = json!({"schema_version":1,"cohort_id":manifest.cohort_id,"study_doi":manifest.study_doi,"score_unit":UNIT,"independent_study_count":1,"assay_count":15,"primary_elements":9,"validation_only_assays":5,"previously_exposed_reference_assays":1,"all_source_rows":source_rows.len(),"reference_validated_snv_records":records.len(),"source_equality_rows":source_rows.iter().filter(|r|r.representation=="source_equality").count(),"row_accounting":accounting,"assays":manifest.assays,"input_metadata_manifest_sha256":manifest.metadata_manifest_sha256,"input_frozen_protocol_sha256":manifest.frozen_protocol_sha256,"limitations":manifest.limitations,"ingestion_policy":"All reference-valid SNVs retained, including missing scores/counts. Quality, folds, validation and incomplete-position rules belong to the frozen analysis protocol; equality/complex source rows remain explicitly accounted."});
        Ok(Self {
            manifest,
            manifest_sha256: digest(&manifest_bytes),
            records,
            source_rows,
            report,
        })
    }
    pub fn records(&self) -> &[RegulatoryRecord] {
        &self.records
    }
    pub fn assays(&self) -> &[RegulatoryAssay] {
        &self.manifest.assays
    }
    pub fn source_rows(&self) -> &[RegulatorySourceRow] {
        &self.source_rows
    }
    pub fn report(&self) -> &Value {
        &self.report
    }
    pub fn manifest_sha256(&self) -> &str {
        &self.manifest_sha256
    }
}

fn validate_metadata(a: &RegulatoryAssay, m: &Value, license: &Value) -> Result<()> {
    ensure!(
        m["urn"] == a.assay_id
            && m["title"] == a.title
            && m["numVariants"].as_u64() == Some(a.expected_source_rows as u64),
        "cohort source identity/count differs"
    );
    ensure!(
        m["license"]["shortName"] == "CC0"
            && m["license"]["version"] == "1.0"
            && m["license"]["link"] == "https://creativecommons.org/publicdomain/zero/1.0/"
            && m["license"] == license["data_license"],
        "source license is not the verified CC0 grant"
    );
    let t = m["targetGenes"]
        .as_array()
        .context("target metadata missing")?;
    ensure!(
        t.len() == 1
            && t[0]["category"] == "regulatory"
            && t[0]["targetSequence"]["taxonomy"]["code"] == 9606
            && t[0]["targetSequence"]["sequenceType"] == "dna"
            && t[0]["targetSequence"]["sequence"] == a.target_sequence,
        "source target identity differs"
    );
    ensure!(
        m["title"]
            .as_str()
            .is_some_and(|s| s.to_ascii_lowercase().contains("promoter")),
        "not a promoter source"
    );
    let method = m["methodText"].as_str().context("method missing")?;
    ensure!(
        method.contains("Log2 variant expression effect") && method.contains("linear model"),
        "source endpoint mismatch"
    );
    ensure!(
        m["shortDescription"] == a.assay_context
            && m["shortDescription"]
                .as_str()
                .is_some_and(|s| s.contains(&a.cell_context)),
        "source context mismatch"
    );
    ensure!(
        m["experiment"]["primaryPublicationIdentifiers"]
            .as_array()
            .is_some_and(|v| v.iter().any(|p| p["identifier"] == "31395865")),
        "unexpected primary study"
    );
    Ok(())
}
fn validate_mapping(a: &RegulatoryAssay, r: &Value, m: &Value) -> Result<()> {
    ensure!(
        a.mapping.assembly == "GRCh38"
            && r["genome"] == "hg38"
            && r["chrom"] == format!("chr{}", a.mapping.chromosome),
        "reference assembly/chromosome mismatch"
    );
    let seq = r["dna"]
        .as_str()
        .context("reference DNA missing")?
        .to_ascii_uppercase();
    let start = r["start"].as_u64().context("reference start missing")?;
    let end = r["end"].as_u64().context("reference end missing")?;
    ensure!(
        start < end && end - start == seq.len() as u64,
        "reference length mismatch"
    );
    let rev = reverse(a.target_sequence.as_bytes())?;
    let mut matches = Vec::new();
    if seq.len() >= a.target_length {
        for (index, window) in seq.as_bytes().windows(a.target_length).enumerate() {
            for (t, strand) in [(a.target_sequence.as_bytes(), "+"), (rev.as_bytes(), "-")] {
                if window == t {
                    matches.push((start + index as u64, strand));
                }
            }
        }
    }
    ensure!(
        matches.len() == 1
            && matches[0].0 == a.mapping.start
            && matches[0].1 == a.mapping.target_reference_alignment_strand
            && a.mapping.end == a.mapping.start + a.target_length as u64,
        "source target does not reconcile uniquely with mapped reference"
    );
    ensure!(
        m["extraMetadata"]["reference"] == "hg38"
            && m["extraMetadata"]["chr"] == a.mapping.chromosome
            && m["extraMetadata"]["start"].as_u64() == Some(a.mapping.start + 1)
            && m["extraMetadata"]["end"].as_u64() == Some(a.mapping.end)
            && !a.mapping.genome_wide_uniqueness_established,
        "coordinate convention mismatch"
    );
    Ok(())
}
type TargetSnv = (usize, u8, u8);
type ParsedHgvs = (&'static str, Option<TargetSnv>);
fn parse_hgvs(s: &str, target: &[u8]) -> Result<ParsedHgvs> {
    let Some(body) = s.strip_prefix("n.") else {
        return Ok(("unsupported_coordinate_system", None));
    };
    let n = body.bytes().take_while(u8::is_ascii_digit).count();
    if n > 0 && !body.starts_with('0') {
        let p: usize = body[..n].parse().context("target coordinate overflow")?;
        ensure!(
            (1..=target.len()).contains(&p),
            "source coordinate outside target"
        );
        let suffix = &body[n..];
        if suffix == "=" {
            return Ok(("source_equality", None));
        }
        let b = suffix.as_bytes();
        if b.len() == 3 && b[1] == b'>' && b"ACGT".contains(&b[0]) && b"ACGT".contains(&b[2]) {
            ensure!(
                target[p - 1] == b[0] && b[0] != b[2],
                "source SNV reference/equality mismatch"
            );
            return Ok(("single_nucleotide_substitution", Some((p, b[0], b[2]))));
        }
    }
    let category = if body.contains("delins")
        || body.contains(';')
        || body.contains('[')
        || body.contains('_')
    {
        "unsupported_complex_allele"
    } else if body.contains("del") {
        "unsupported_deletion"
    } else if body.contains("ins") || body.contains("dup") {
        "unsupported_insertion_or_duplication"
    } else {
        "unsupported_or_malformed_description"
    };
    Ok((category, None))
}
fn map_snv(m: &Mapping, p: usize, r: u8, a: u8) -> Result<Variant<GRCh38>> {
    let (offset, r, a) = match m.target_reference_alignment_strand.as_str() {
        "+" => (m.start + p as u64 - 1, r, a),
        "-" => (m.end - p as u64, complement(r)?, complement(a)?),
        _ => bail!("unsupported target strand"),
    };
    Variant::snv(
        GenomicPosition::new(&m.chromosome, offset)?,
        &(r as char).to_string(),
        &(a as char).to_string(),
    )
}
fn complement(b: u8) -> Result<u8> {
    match b {
        b'A' => Ok(b'T'),
        b'C' => Ok(b'G'),
        b'G' => Ok(b'C'),
        b'T' => Ok(b'A'),
        _ => bail!("ambiguous DNA"),
    }
}
fn reverse(s: &[u8]) -> Result<String> {
    s.iter()
        .rev()
        .map(|&b| complement(b).map(char::from))
        .collect::<Result<String>>()
}
fn number(raw: &str, domain: &str) -> SourceNumber {
    let t = raw.trim();
    let (value, status) = if t.is_empty()
        || ["NA", "N/A", "NULL", "."]
            .iter()
            .any(|x| t.eq_ignore_ascii_case(x))
    {
        (None, "missing")
    } else {
        match t.parse::<f64>() {
            Err(_) => (None, "non_numeric"),
            Ok(n) if !n.is_finite() => (None, "non_finite"),
            Ok(n) => {
                if (domain == "probability" && !(0.0..=1.0).contains(&n))
                    || (domain == "count"
                        && (n < 0.0 || n.fract() != 0.0 || n > (1u64 << 53) as f64))
                {
                    (None, "outside_domain")
                } else {
                    (Some(n), "finite")
                }
            }
        }
    };
    SourceNumber {
        raw: raw.into(),
        value,
        status: status.into(),
    }
}
fn digest(b: &[u8]) -> String {
    format!("{:x}", Sha256::digest(b))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn strict_target_coordinates_and_reference_orientation() {
        assert_eq!(
            parse_hgvs("n.2C>T", b"ACGT").unwrap().1,
            Some((2, b'C', b'T'))
        );
        assert_eq!(parse_hgvs("n.2=", b"ACGT").unwrap().0, "source_equality");
        assert!(parse_hgvs("n.2A>T", b"ACGT").is_err());
        assert!(parse_hgvs("n.5A>T", b"ACGT").is_err());
        assert!(parse_hgvs("n.02C>T", b"ACGT").unwrap().1.is_none());
        assert!(parse_hgvs("n.[1A>C;2C>T]", b"ACGT").unwrap().1.is_none());
        let m = Mapping {
            assembly: "GRCh38".into(),
            chromosome: "11".into(),
            start: 100,
            end: 104,
            target_reference_alignment_strand: "-".into(),
            source_coordinate_convention: String::new(),
            genome_wide_uniqueness_established: false,
        };
        assert_eq!(
            map_snv(&m, 2, b'C', b'A').unwrap().key(),
            "GRCh38:11:102:G>T"
        );
    }
    #[test]
    fn missing_and_invalid_counts_are_not_silently_imputed() {
        assert_eq!(number("NaN", "score").status, "non_finite");
        assert_eq!(number("NA", "score").status, "missing");
        assert_eq!(number("10.5", "count").status, "outside_domain");
        assert_eq!(number("-0.5", "score").value, Some(-0.5));
    }
    #[test]
    fn frozen_cohort_maps_all_rows_and_preserves_roles() {
        let a = RegulatoryCohort::open(
            &Path::new(env!("CARGO_MANIFEST_DIR")).join("../../datasets/regulatory-cohort"),
        )
        .unwrap();
        assert_eq!(a.assays().len(), 15);
        assert_eq!(a.assays().iter().filter(|a| a.role == "primary").count(), 9);
        assert!(a.records().len() > 10000);
        assert!(a.source_rows().len() > a.records().len());
        assert_eq!(a.manifest_sha256().len(), 64);
    }
    fn copied_fixture() -> tempfile::TempDir {
        let temp = tempfile::tempdir().unwrap();
        let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../datasets/regulatory-cohort");
        for entry in fs::read_dir(source).unwrap() {
            let entry = entry.unwrap();
            if entry.file_type().unwrap().is_file() {
                fs::copy(entry.path(), temp.path().join(entry.file_name())).unwrap();
            }
        }
        temp
    }
    #[test]
    fn source_tampering_and_post_outcome_role_changes_are_rejected() {
        let temp = copied_fixture();
        fs::write(temp.path().join("00000015-a-1.scores.csv"), "tampered").unwrap();
        assert!(RegulatoryCohort::open(temp.path())
            .err()
            .unwrap()
            .to_string()
            .contains("checksum mismatch"));
        let temp = copied_fixture();
        let path = temp.path().join("manifest.json");
        let mut manifest: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        manifest["assays"][0]["role"] = json!("validation_only");
        fs::write(path, serde_json::to_vec(&manifest).unwrap()).unwrap();
        assert!(RegulatoryCohort::open(temp.path())
            .err()
            .unwrap()
            .to_string()
            .contains("roles differ"));
    }
    #[test]
    fn missing_score_keeps_source_row_and_snv_identity() {
        let temp = copied_fixture();
        let filename = "00000015-a-1.scores.csv";
        let mut reader = csv::Reader::from_path(temp.path().join(filename)).unwrap();
        let mut rows: Vec<Vec<String>> = reader
            .records()
            .map(|r| r.unwrap().iter().map(str::to_owned).collect())
            .collect();
        let index = rows.iter().position(|r| r[1].contains('>')).unwrap();
        let id = rows[index][0].clone();
        rows[index][4] = "NA".into();
        let mut writer = csv::Writer::from_writer(Vec::new());
        writer.write_record(HEADERS).unwrap();
        for row in rows {
            writer.write_record(row).unwrap();
        }
        let bytes = writer.into_inner().unwrap();
        fs::write(temp.path().join(filename), &bytes).unwrap();
        let path = temp.path().join("manifest.json");
        let mut manifest: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        manifest["files"][filename] = json!(digest(&bytes));
        for source in manifest["sources"].as_array_mut().unwrap() {
            if source["file"] == filename {
                source["sha256"] = json!(digest(&bytes));
                source["bytes"] = json!(bytes.len());
            }
        }
        manifest["row_accounting"][0]["score_status_counts"] = json!({"finite":983,"missing":1});
        fs::write(path, serde_json::to_vec(&manifest).unwrap()).unwrap();
        let cohort = RegulatoryCohort::open(temp.path()).unwrap();
        assert_eq!(cohort.records().len(), 15699);
        assert_eq!(cohort.source_rows().len(), 18468);
        let record = cohort.records().iter().find(|r| r.id == id).unwrap();
        assert!(record.score.is_none());
        assert_eq!(record.score_status, "missing");
    }
}
