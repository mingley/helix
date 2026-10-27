//! Checked, hash-pinned public data ingestion. Model inputs exclude held-out labels.
use anyhow::{bail, ensure, Context, Result};
use genome_types::{
    BiologicalContext, ClinicalLabel, Consequence, FeatureRecord, GRCh38, GenomicPosition, Variant,
};
use noodles_vcf as vcf;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::BufReader,
    path::{Path, PathBuf},
};

pub fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    pub schema_version: u32,
    pub assembly: String,
    pub gene: String,
    pub reference_accession: String,
    pub reference_start: u64,
    pub reference_end: u64,
    pub source_url: String,
    pub reference_url: String,
    pub query: String,
    pub retrieved_at: String,
    pub split: String,
    pub scope: String,
    pub files: BTreeMap<String, String>,
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct LabelRecord {
    pub id: String,
    pub label: ClinicalLabel,
    pub classification: String,
    pub review_status: String,
    pub source_url: String,
}
#[derive(Clone, Debug)]
pub struct TrainingExample {
    pub feature: FeatureRecord,
    pub observation: LabelRecord,
}

pub struct Dataset {
    pub manifest: Manifest,
    pub manifest_sha256: String,
    pub features: Vec<FeatureRecord>,
    pub training: Vec<TrainingExample>,
    heldout: Vec<LabelRecord>,
}
impl Dataset {
    pub fn open(directory: &Path) -> Result<Self> {
        let root = directory.canonicalize().context("dataset directory")?;
        let manifest_bytes = fs::read(root.join("manifest.json"))?;
        let manifest: Manifest = serde_json::from_slice(&manifest_bytes)?;
        ensure!(
            manifest.schema_version == 1 && manifest.assembly == "GRCh38" && manifest.gene == "HBB",
            "unsupported dataset schema, assembly, or gene"
        );
        ensure!(
            manifest.reference_accession == "NC_000011.10",
            "reference accession must identify GRCh38 chromosome 11"
        );
        ensure!(
            manifest.reference_start < manifest.reference_end
                && manifest.reference_end <= 135_086_622,
            "invalid reference interval"
        );
        let required = [
            "heldout.json",
            "reference.fasta",
            "source.json",
            "training.json",
            "variants.vcf",
        ];
        ensure!(
            manifest
                .files
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>()
                == required,
            "manifest must pin exactly the supported input files"
        );
        let mut contents = BTreeMap::new();
        for (name, expected) in &manifest.files {
            let path = checked_path(&root, name)?;
            let bytes = fs::read(path)?;
            ensure!(sha256(&bytes) == *expected, "checksum mismatch for {name}");
            contents.insert(name.as_str(), bytes);
        }
        let fasta = std::str::from_utf8(&contents["reference.fasta"])?;
        let header = fasta.lines().next().context("empty FASTA")?;
        let expected_header = format!(
            ">NC_000011.10:{}-{}",
            manifest.reference_start + 1,
            manifest.reference_end
        );
        ensure!(
            header.starts_with(&format!("{expected_header} ")) || header == expected_header,
            "FASTA header does not match manifest reference interval"
        );
        let reference: String = fasta
            .lines()
            .skip(1)
            .collect::<String>()
            .to_ascii_uppercase();
        ensure!(
            reference.len() as u64 == manifest.reference_end - manifest.reference_start,
            "FASTA length mismatch"
        );
        ensure!(
            reference.bytes().all(|b| b"ACGT".contains(&b)),
            "ambiguous reference sequence"
        );
        let raw_vcf = std::str::from_utf8(&contents["variants.vcf"])?;
        ensure!(
            raw_vcf.lines().any(|line| line == "##reference=GRCh38"),
            "VCF must declare GRCh38"
        );
        let mut reader = vcf::io::Reader::new(BufReader::new(raw_vcf.as_bytes()));
        let header = reader.read_header()?;
        let mut features = Vec::new();
        let mut ids = BTreeSet::new();
        let mut keys = BTreeSet::new();
        for record in reader.record_bufs(&header) {
            let record = record.context("invalid VCF record")?;
            ensure!(
                record.reference_sequence_name() == "11",
                "this slice supports HBB chromosome 11 only"
            );
            let position =
                usize::from(record.variant_start().context("missing VCF position")?) as u64;
            let offset = position
                .checked_sub(1)
                .context("VCF position must be positive")?;
            ensure!(
                (manifest.reference_start..manifest.reference_end).contains(&offset),
                "variant outside reference slice"
            );
            let alternate: &[String] = record.alternate_bases().as_ref();
            ensure!(
                alternate.len() == 1,
                "multiallelic VCF records require explicit normalization"
            );
            let variant = Variant::<GRCh38>::snv(
                GenomicPosition::new("11", offset)?,
                record.reference_bases(),
                &alternate[0],
            )?;
            ensure!(
                reference.as_bytes()[(offset - manifest.reference_start) as usize] as char
                    == variant.reference(),
                "reference allele mismatch at {}",
                variant.key()
            );
            ensure!(
                record.ids().as_ref().len() == 1,
                "each VCF record needs one versioned ClinVar ID"
            );
            let id = record.ids().as_ref().iter().next().unwrap().clone();
            ensure!(
                ids.insert(id.clone()) && keys.insert(variant.key()),
                "duplicate VCF ID or variant"
            );
            let consequence = match record.info().get("MC") {
                Some(Some(vcf::variant::record_buf::info::field::Value::String(value))) => {
                    serde_json::from_value::<Consequence>(serde_json::Value::String(value.clone()))?
                }
                _ => bail!("missing molecular consequence for {id}"),
            };
            features.push(FeatureRecord {
                id,
                variant,
                consequence,
                context: BiologicalContext::hbb_germline(),
            });
        }
        ensure!(!features.is_empty(), "empty VCF");
        let training_labels: Vec<LabelRecord> = serde_json::from_slice(&contents["training.json"])?;
        let heldout: Vec<LabelRecord> = serde_json::from_slice(&contents["heldout.json"])?;
        ensure!(
            !training_labels.is_empty() && !heldout.is_empty(),
            "both partitions must be nonempty"
        );
        let source: serde_json::Value = serde_json::from_slice(&contents["source.json"])?;
        let by_id: BTreeMap<_, _> = features.iter().map(|f| (f.id.as_str(), f)).collect();
        let mut assigned = BTreeSet::new();
        let mut loci = BTreeMap::new();
        for (name, labels) in [("training", &training_labels), ("heldout", &heldout)] {
            for label in labels {
                ensure!(
                    assigned.insert(label.id.as_str()),
                    "duplicate or overlapping train/heldout ID: {}",
                    label.id
                );
                let feature = by_id
                    .get(label.id.as_str())
                    .context("label has no VCF variant")?;
                let locus = feature.variant.position().offset();
                if let Some(previous) = loci.insert(locus, name) {
                    ensure!(previous == name, "locus leakage across partitions");
                }
                let group = Sha256::digest(format!("NC_000011.10:{locus}").as_bytes())[0] % 2;
                ensure!(
                    (group == 0) == (name == "training"),
                    "partition assignment violates frozen locus split"
                );
                validate_source(label, feature, &source)?;
            }
        }
        ensure!(
            assigned.len() == features.len(),
            "every variant must belong to exactly one partition"
        );
        let training = training_labels
            .into_iter()
            .map(|observation| TrainingExample {
                feature: by_id[observation.id.as_str()].clone(),
                observation,
            })
            .collect();
        Ok(Self {
            manifest,
            manifest_sha256: sha256(&manifest_bytes),
            features,
            training,
            heldout,
        })
    }
    /// Read only after prediction. Agents receive training data and label-free features.
    pub fn heldout(&self) -> &[LabelRecord] {
        &self.heldout
    }
    pub fn heldout_features(&self) -> Vec<FeatureRecord> {
        let training: BTreeSet<_> = self.training.iter().map(|t| &t.feature.id).collect();
        self.features
            .iter()
            .filter(|f| !training.contains(&f.id))
            .cloned()
            .collect()
    }
}

fn checked_path(root: &Path, name: &str) -> Result<PathBuf> {
    ensure!(
        Path::new(name).components().count() == 1,
        "dataset paths must be basenames"
    );
    let path = root.join(name).canonicalize()?;
    ensure!(path.starts_with(root), "dataset symlink escapes root");
    Ok(path)
}
fn validate_source(
    label: &LabelRecord,
    feature: &FeatureRecord,
    source: &serde_json::Value,
) -> Result<()> {
    let uid = label
        .source_url
        .strip_prefix("https://www.ncbi.nlm.nih.gov/clinvar/variation/")
        .and_then(|s| s.strip_suffix('/'))
        .context("invalid ClinVar provenance URL")?;
    ensure!(
        !uid.is_empty() && uid.bytes().all(|b| b.is_ascii_digit()),
        "invalid ClinVar accession URL"
    );
    let record = &source["result"][uid];
    ensure!(
        record["genes"]
            .as_array()
            .is_some_and(|genes| genes.iter().any(|g| g["symbol"] == "HBB")),
        "source does not identify HBB"
    );
    ensure!(
        record["accession_version"].as_str() == Some(label.id.as_str()),
        "source accession mismatch"
    );
    ensure!(
        record["germline_classification"]["description"].as_str()
            == Some(label.classification.as_str()),
        "source classification mismatch"
    );
    ensure!(
        record["germline_classification"]["review_status"].as_str()
            == Some(label.review_status.as_str()),
        "source review status mismatch"
    );
    let expected = match label.classification.as_str() {
        "Benign" | "Likely benign" => ClinicalLabel::Benign,
        "Pathogenic" | "Likely pathogenic" => ClinicalLabel::Pathogenic,
        _ => bail!("uncertain/conflicting classifications cannot be binary truth"),
    };
    ensure!(
        label.label == expected,
        "binary label disagrees with source"
    );
    let v = &feature.variant;
    let spdi = format!(
        "NC_000011.10:{}:{}:{}",
        v.position().offset(),
        v.reference(),
        v.alternate()
    );
    ensure!(
        record["variation_set"][0]["canonical_spdi"].as_str() == Some(spdi.as_str()),
        "source SPDI mismatch"
    );
    let names = record["molecular_consequence_list"]
        .as_array()
        .context("missing source molecular consequences")?
        .iter()
        .map(|name| name.as_str().context("invalid molecular consequence"))
        .collect::<Result<Vec<_>>>()?
        .join(" ")
        .to_ascii_lowercase();
    let expected_consequence = [
        ("splice", Consequence::Splice),
        ("stop", Consequence::Stop),
        ("missense", Consequence::Missense),
        ("synonymous", Consequence::Synonymous),
        ("intron", Consequence::Intron),
        ("utr", Consequence::Utr),
    ]
    .into_iter()
    .find(|(term, _)| names.contains(term))
    .map(|(_, group)| group)
    .unwrap_or(Consequence::Other);
    ensure!(
        feature.consequence == expected_consequence,
        "molecular consequence disagrees with source"
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../datasets/clinvar-hbb")
    }
    #[test]
    fn public_fixture_loads_and_has_disjoint_loci() {
        let ds = Dataset::open(&fixture()).unwrap();
        assert!(!ds.training.is_empty() && !ds.heldout().is_empty());
        assert_eq!(ds.features.len(), ds.training.len() + ds.heldout().len());
    }
    #[test]
    fn tampering_is_rejected_before_ingestion() {
        let temp = tempfile::tempdir().unwrap();
        for entry in fs::read_dir(fixture()).unwrap() {
            let entry = entry.unwrap();
            if entry.file_type().unwrap().is_file() {
                fs::copy(entry.path(), temp.path().join(entry.file_name())).unwrap();
            }
        }
        fs::write(temp.path().join("variants.vcf"), "tampered").unwrap();
        assert!(Dataset::open(temp.path())
            .err()
            .unwrap()
            .to_string()
            .contains("checksum mismatch"));
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
        let mut manifest: Manifest = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        manifest.files.insert(name.into(), sha256(bytes));
        fs::write(path, serde_json::to_vec(&manifest).unwrap()).unwrap();
    }
    #[test]
    fn a_new_hash_does_not_hide_wrong_alleles_or_source_annotations() {
        let temp = copied_fixture();
        let vcf = fs::read_to_string(temp.path().join("variants.vcf")).unwrap();
        let first = vcf.lines().find(|s| !s.starts_with('#')).unwrap();
        let mut columns: Vec<_> = first.split('\t').map(str::to_owned).collect();
        columns[3] = "G".into(); // First curated record is C>A.
        let mutated = vcf.replacen(first, &columns.join("\t"), 1);
        recommit(temp.path(), "variants.vcf", mutated.as_bytes());
        assert!(Dataset::open(temp.path())
            .err()
            .unwrap()
            .to_string()
            .contains("reference allele mismatch"));
        let mutated = vcf.replacen("MC=utr", "MC=other", 1);
        recommit(temp.path(), "variants.vcf", mutated.as_bytes());
        assert!(Dataset::open(temp.path())
            .err()
            .unwrap()
            .to_string()
            .contains("molecular consequence disagrees"));
    }
    #[test]
    fn overlap_and_labels_inconsistent_with_source_are_rejected() {
        let temp = copied_fixture();
        let mut train: Vec<LabelRecord> =
            serde_json::from_slice(&fs::read(temp.path().join("training.json")).unwrap()).unwrap();
        train[0].label = if train[0].label.positive() {
            ClinicalLabel::Benign
        } else {
            ClinicalLabel::Pathogenic
        };
        recommit(
            temp.path(),
            "training.json",
            &serde_json::to_vec(&train).unwrap(),
        );
        assert!(Dataset::open(temp.path())
            .err()
            .unwrap()
            .to_string()
            .contains("binary label disagrees"));
        let temp = copied_fixture();
        let train: Vec<LabelRecord> =
            serde_json::from_slice(&fs::read(temp.path().join("training.json")).unwrap()).unwrap();
        let mut holdout: Vec<LabelRecord> =
            serde_json::from_slice(&fs::read(temp.path().join("heldout.json")).unwrap()).unwrap();
        holdout.push(train[0].clone());
        recommit(
            temp.path(),
            "heldout.json",
            &serde_json::to_vec(&holdout).unwrap(),
        );
        assert!(Dataset::open(temp.path())
            .err()
            .unwrap()
            .to_string()
            .contains("overlapping"));
    }
}
