//! Frozen prior-biology PFM features, independent of assay outcomes.
use anyhow::{ensure, Context, Result};
use serde::Serialize;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs,
    path::{Component, Path},
};

#[derive(Clone, Debug)]
pub struct Motif {
    pub matrix_id: String,
    pub name: String,
    pub width: usize,
    log_odds: Vec<[f64; 4]>,
}

#[derive(Clone, Debug, Serialize)]
pub struct MotifHit {
    pub start_0based: usize,
    pub end_0based: usize,
    pub strand: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct MotifEffect {
    pub matrix_id: String,
    pub reference_max_log2_odds: f64,
    pub alternate_max_log2_odds: f64,
    pub delta_max_log2_odds: f64,
    pub reference_argmax: MotifHit,
    pub alternate_argmax: MotifHit,
    pub candidate_window_orientation_count: usize,
    pub support_start_0based: usize,
    pub support_end_0based: usize,
}

impl Motif {
    pub fn from_counts(matrix_id: &str, name: &str, counts: &[Vec<f64>; 4]) -> Result<Self> {
        let width = counts[0].len();
        ensure!(
            !matrix_id.is_empty() && width > 0 && width <= 32,
            "invalid motif identity/width"
        );
        ensure!(
            counts
                .iter()
                .all(|row| row.len() == width && row.iter().all(|v| v.is_finite() && *v >= 0.0)),
            "invalid motif count matrix"
        );
        let mut log_odds = Vec::new();
        for column in 0..width {
            let total = counts.iter().map(|row| row[column]).sum::<f64>();
            ensure!(total.is_finite(), "motif column count overflow");
            let mut weights = [0.0; 4];
            for base in 0..4 {
                weights[base] = ((counts[base][column] + 0.25) / (total + 1.0) / 0.25).log2();
            }
            log_odds.push(weights);
        }
        Ok(Self {
            matrix_id: matrix_id.into(),
            name: name.into(),
            width,
            log_odds,
        })
    }

    /// All complete overlapping windows in both orientations, same candidate
    /// set for both alleles. Maxima are independent; no binding threshold.
    /// None explicitly denotes a target too short for a complete motif window.
    pub fn effect(
        &self,
        sequence: &str,
        position_1based: usize,
        alternate: &str,
    ) -> Result<Option<MotifEffect>> {
        let sequence = sequence.as_bytes();
        ensure!(
            !sequence.is_empty() && sequence.iter().all(|b| b"ACGT".contains(b)),
            "motif input requires unambiguous uppercase DNA"
        );
        ensure!(
            (1..=sequence.len()).contains(&position_1based),
            "motif target position out of bounds"
        );
        let alternate = alternate.as_bytes();
        ensure!(
            alternate.len() == 1 && b"ACGT".contains(&alternate[0]),
            "motif alternate must be one DNA base"
        );
        if sequence.len() < self.width {
            return Ok(None);
        }
        let position = position_1based - 1;
        let start_min = (position + 1).saturating_sub(self.width);
        let start_max = position.min(sequence.len() - self.width);
        let mut reference_best: Option<(f64, MotifHit)> = None;
        let mut alternate_best: Option<(f64, MotifHit)> = None;
        let mut candidates = 0;
        for start in start_min..=start_max {
            for strand in ["+", "-"] {
                let mut original = 0.0;
                let mut changed = 0.0;
                for column in 0..self.width {
                    let original_base = sequence[start + column];
                    let changed_base = if start + column == position {
                        alternate[0]
                    } else {
                        original_base
                    };
                    let (column, original_base, changed_base) = if strand == "+" {
                        (column, original_base, changed_base)
                    } else {
                        (
                            self.width - 1 - column,
                            complement(original_base),
                            complement(changed_base),
                        )
                    };
                    original += self.log_odds[column][base_index(original_base)?];
                    changed += self.log_odds[column][base_index(changed_base)?];
                }
                let hit = MotifHit {
                    start_0based: start,
                    end_0based: start + self.width,
                    strand: strand.into(),
                };
                // Strict comparison preserves the first start/orientation on exact ties.
                if reference_best
                    .as_ref()
                    .is_none_or(|(score, _)| original > *score)
                {
                    reference_best = Some((original, hit.clone()));
                }
                if alternate_best
                    .as_ref()
                    .is_none_or(|(score, _)| changed > *score)
                {
                    alternate_best = Some((changed, hit));
                }
                candidates += 1;
            }
        }
        let (reference_max_log2_odds, reference_argmax) =
            reference_best.context("no reference motif candidates")?;
        let (alternate_max_log2_odds, alternate_argmax) =
            alternate_best.context("no alternate motif candidates")?;
        Ok(Some(MotifEffect {
            matrix_id: self.matrix_id.clone(),
            reference_max_log2_odds,
            alternate_max_log2_odds,
            delta_max_log2_odds: alternate_max_log2_odds - reference_max_log2_odds,
            reference_argmax,
            alternate_argmax,
            candidate_window_orientation_count: candidates,
            support_start_0based: start_min,
            support_end_0based: start_max + self.width,
        }))
    }
}

pub struct MotifSet {
    pub motifs: Vec<Motif>,
    pub metadata: Value,
    manifest_sha256: String,
}
impl MotifSet {
    pub fn open(directory: &Path) -> Result<Self> {
        let root = directory.canonicalize()?;
        let manifest_bytes = fs::read(root.join("manifest.json"))?;
        let manifest: Value = serde_json::from_slice(&manifest_bytes)?;
        ensure!(
            manifest["schema_version"] == 1 && manifest["motifs_file"] == "motifs.json",
            "unsupported motif manifest"
        );
        let files = manifest["files"]
            .as_object()
            .context("motif file commitments missing")?;
        let mut contents = BTreeMap::new();
        for (name, hash) in files {
            let mut parts = Path::new(name).components();
            ensure!(
                matches!(parts.next(), Some(Component::Normal(_))) && parts.next().is_none(),
                "motif files must be basenames"
            );
            let path = root.join(name).canonicalize()?;
            ensure!(path.starts_with(&root), "motif file escapes bundle");
            let bytes = fs::read(path)?;
            ensure!(
                hash.as_str() == Some(digest(&bytes).as_str()),
                "motif checksum mismatch: {name}"
            );
            contents.insert(name.as_str(), bytes);
        }
        ensure!(
            manifest["terms"]["license"] == "CC BY 4.0",
            "unsupported motif source license"
        );
        let metadata: Value = serde_json::from_slice(
            contents
                .get("motifs.json")
                .context("missing motif definitions")?,
        )?;
        ensure!(
            metadata["schema_version"] == 1
                && metadata["scoring"]["id"] == "jaspar-overlap-max-delta-v1"
                && metadata["scoring"]["pseudocount_total_per_column"] == 1.0
                && metadata["scoring"]["maximum_motif_width"] == 11
                && metadata["scoring"]["strict_disjoint_support_center_guard"] == 20,
            "unsupported motif scoring definition"
        );
        for base in ["A", "C", "G", "T"] {
            ensure!(
                metadata["scoring"]["background"][base] == 0.25,
                "unsupported motif background"
            );
        }
        let definitions = metadata["motifs"].as_array().context("missing motifs")?;
        let expected = ["MA0035.4", "MA0493.3", "MA0060.3", "MA0108.3"];
        ensure!(
            definitions.len() == expected.len(),
            "unexpected prior motif count"
        );
        let mut motifs = Vec::new();
        for (definition, expected_id) in definitions.iter().zip(expected) {
            ensure!(
                definition["matrix_id"] == expected_id,
                "prior motif identity/order changed"
            );
            let file = definition["source_file"]
                .as_str()
                .context("motif source file missing")?;
            let raw: Value =
                serde_json::from_slice(contents.get(file).context("unpinned motif source")?)?;
            ensure!(
                raw["matrix_id"] == definition["matrix_id"]
                    && raw["name"] == definition["name"]
                    && raw["pfm"] == definition["pfm"],
                "projected motif differs from raw JASPAR source"
            );
            let rows: Vec<Vec<f64>> = ["A", "C", "G", "T"]
                .into_iter()
                .map(|base| {
                    definition["pfm"][base]
                        .as_array()
                        .context("PFM row missing")?
                        .iter()
                        .map(|v| v.as_f64().context("invalid PFM count"))
                        .collect()
                })
                .collect::<Result<_>>()?;
            let counts: [Vec<f64>; 4] = rows
                .try_into()
                .map_err(|_| anyhow::anyhow!("PFM must have four rows"))?;
            let motif = Motif::from_counts(
                expected_id,
                definition["name"].as_str().context("motif name missing")?,
                &counts,
            )?;
            ensure!(
                definition["width"].as_u64() == Some(motif.width as u64),
                "motif width mismatch"
            );
            motifs.push(motif);
        }
        Ok(Self {
            motifs,
            metadata,
            manifest_sha256: digest(&manifest_bytes),
        })
    }
    pub fn manifest_sha256(&self) -> &str {
        &self.manifest_sha256
    }
    pub fn report(&self) -> Value {
        json!({"manifest_sha256":self.manifest_sha256,"scoring":self.metadata["scoring"],
            "motifs":self.metadata["motifs"],"interpretation":"Prior motif log-odds similarity features are not measured binding, occupancy, energy, expression or clinical effects"})
    }
}

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn base_index(base: u8) -> Result<usize> {
    b"ACGT"
        .iter()
        .position(|&b| b == base)
        .context("invalid DNA base")
}
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => unreachable!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn toy() -> Motif {
        Motif::from_counts(
            "toy",
            "AC",
            &[
                vec![10.0, 0.0],
                vec![0.0, 10.0],
                vec![0.0, 0.0],
                vec![0.0, 0.0],
            ],
        )
        .unwrap()
    }
    #[test]
    fn pseudocounts_and_allele_delta_match_hand_calculation() {
        let effect = toy().effect("AC", 1, "T").unwrap().unwrap();
        assert!((effect.reference_max_log2_odds - 2.0 * (41.0f64 / 11.0).log2()).abs() < 1e-12);
        assert!((effect.delta_max_log2_odds + 41.0f64.log2()).abs() < 1e-12);
        assert_eq!(effect.candidate_window_orientation_count, 2);
        assert_eq!(effect.reference_argmax.start_0based, 0);
    }
    #[test]
    fn scores_respect_both_strands_and_identical_window_candidates() {
        let motif = toy();
        let a = motif.effect("AACGTT", 2, "T").unwrap().unwrap();
        let b = motif.effect("AACGTT", 5, "A").unwrap().unwrap(); // reverse-complement palindrome; transformed allele
        assert!((a.reference_max_log2_odds - b.reference_max_log2_odds).abs() < 1e-12);
        assert!((a.alternate_max_log2_odds - b.alternate_max_log2_odds).abs() < 1e-12);
        let unchanged = motif.effect("AACGTT", 2, "A").unwrap().unwrap();
        assert_eq!(unchanged.delta_max_log2_odds, 0.0);
        assert!(motif.effect("A", 1, "C").unwrap().is_none());
        assert!(motif.effect("AN", 1, "C").is_err());
        assert!(motif.effect("AC", 0, "C").is_err());
    }
    #[test]
    fn frozen_prior_matrices_are_source_checked_without_assay_data() {
        let set = MotifSet::open(
            &Path::new(env!("CARGO_MANIFEST_DIR")).join("../../datasets/regulatory-motifs"),
        )
        .unwrap();
        assert_eq!(set.motifs.len(), 4);
        assert_eq!(set.motifs.iter().map(|m| m.width).max(), Some(11));
        assert_eq!(set.metadata["motifs"][3]["species"], json!([]));
        assert!(set.metadata["motifs"][3]["species_caveat"].is_string());
    }
}
