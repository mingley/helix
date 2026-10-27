//! Paired sequence/context metadata, separately supplied outcomes, and group splits.
use crate::features::{canonical_reverse_complement, sequence_features, SequenceFeatures};
use anyhow::{ensure, Context, Result};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug, Serialize)]
pub struct SequenceRecord {
    pub id: String,
    /// A protocol-defined leakage group, shared by every occurrence of this DNA.
    pub group: String,
    pub sequence: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct MatchedPanel {
    sequences: Vec<SequenceRecord>,
    contexts: Vec<String>,
    sequence_features: Vec<SequenceFeatures>,
    metadata_sha256: String,
}

impl MatchedPanel {
    /// Metadata only. Every sequence is intended to have one outcome per context.
    /// The separate OutcomeMatrix constructor enforces that complete pairing.
    pub fn new(mut sequences: Vec<SequenceRecord>, mut contexts: Vec<String>) -> Result<Self> {
        ensure!(!sequences.is_empty(), "sequence panel is empty");
        ensure!(!contexts.is_empty(), "context list is empty");
        contexts.sort();
        ensure!(
            contexts.windows(2).all(|pair| pair[0] != pair[1]),
            "duplicate context"
        );
        ensure!(
            contexts.iter().all(|context| valid_label(context)),
            "empty or padded context label"
        );
        sequences.sort_by(|left, right| left.id.cmp(&right.id));
        ensure!(
            sequences.windows(2).all(|pair| pair[0].id != pair[1].id),
            "duplicate sequence ID"
        );
        let mut canonical_groups: BTreeMap<String, (&str, &str)> = BTreeMap::new();
        let mut features = Vec::with_capacity(sequences.len());
        for sequence in &sequences {
            ensure!(
                valid_label(&sequence.id) && valid_label(&sequence.group),
                "empty or padded sequence/group label"
            );
            let canonical = canonical_reverse_complement(&sequence.sequence)
                .with_context(|| format!("invalid sequence {}", sequence.id))?;
            if let Some((group, previous_id)) = canonical_groups.get(&canonical) {
                ensure!(
                    *group == sequence.group,
                    "identical/RC-identical DNA crosses groups: {} and {}",
                    previous_id,
                    sequence.id
                );
            } else {
                canonical_groups.insert(canonical, (&sequence.group, &sequence.id));
            }
            features.push(sequence_features(&sequence.sequence)?);
        }
        let mut hasher = Sha256::new();
        hash_field(&mut hasher, "matched-panel-metadata-v1");
        hasher.update((contexts.len() as u64).to_le_bytes());
        for context in &contexts {
            hash_field(&mut hasher, context);
        }
        hash_field(&mut hasher, "sequence-records");
        hasher.update((sequences.len() as u64).to_le_bytes());
        for sequence in &sequences {
            for value in [&sequence.id, &sequence.group, &sequence.sequence] {
                hash_field(&mut hasher, value);
            }
        }
        Ok(Self {
            sequences,
            contexts,
            sequence_features: features,
            metadata_sha256: format!("{:x}", hasher.finalize()),
        })
    }

    pub fn sequences(&self) -> &[SequenceRecord] {
        &self.sequences
    }
    pub fn contexts(&self) -> &[String] {
        &self.contexts
    }
    pub fn features(&self) -> &[SequenceFeatures] {
        &self.sequence_features
    }
    pub fn metadata_sha256(&self) -> &str {
        &self.metadata_sha256
    }
    pub fn groups(&self) -> Vec<String> {
        self.sequences
            .iter()
            .map(|sequence| sequence.group.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }
    pub fn sequence_index(&self, id: &str) -> Result<usize> {
        self.sequences
            .binary_search_by(|sequence| sequence.id.as_str().cmp(id))
            .map_err(|_| anyhow::anyhow!("unknown sequence ID: {id}"))
    }
    pub fn context_index(&self, context: &str) -> Result<usize> {
        self.contexts
            .binary_search_by(|value| value.as_str().cmp(context))
            .map_err(|_| anyhow::anyhow!("unknown context: {context}"))
    }
}

fn valid_label(value: &str) -> bool {
    !value.is_empty() && value.trim() == value
}
fn hash_field(hasher: &mut Sha256, value: &str) {
    hasher.update((value.len() as u64).to_le_bytes());
    hasher.update(value.as_bytes());
}

#[derive(Clone, Debug, Serialize)]
pub struct MeasuredOutcome {
    pub sequence_id: String,
    pub context: String,
    pub value: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct OutcomeMatrix {
    panel_metadata_sha256: String,
    /// Sequence rows and context columns follow the panel's sorted identifiers.
    values: Vec<Vec<f64>>,
    outcomes_sha256: String,
}

impl OutcomeMatrix {
    pub fn from_records(panel: &MatchedPanel, records: &[MeasuredOutcome]) -> Result<Self> {
        let mut values = vec![vec![None; panel.contexts.len()]; panel.sequences.len()];
        for record in records {
            ensure!(
                record.value.is_finite(),
                "nonfinite outcome for {} / {}",
                record.sequence_id,
                record.context
            );
            let row = panel.sequence_index(&record.sequence_id)?;
            let column = panel.context_index(&record.context)?;
            ensure!(
                values[row][column].replace(record.value).is_none(),
                "duplicate outcome for {} / {}",
                record.sequence_id,
                record.context
            );
        }
        let values: Vec<Vec<f64>> = values
            .into_iter()
            .enumerate()
            .map(|(i, row)| {
                row.into_iter()
                    .enumerate()
                    .map(|(j, value)| {
                        value.with_context(|| {
                            format!(
                                "missing paired outcome for {} / {}",
                                panel.sequences[i].id, panel.contexts[j]
                            )
                        })
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<_>>()?;
        let mut hasher = Sha256::new();
        hash_field(&mut hasher, "matched-panel-outcomes-f64-v1");
        hash_field(&mut hasher, panel.metadata_sha256());
        for row in &values {
            for value in row {
                hasher.update(value.to_bits().to_le_bytes());
            }
        }
        Ok(Self {
            panel_metadata_sha256: panel.metadata_sha256.clone(),
            values,
            outcomes_sha256: format!("{:x}", hasher.finalize()),
        })
    }

    pub fn values(&self) -> &[Vec<f64>] {
        &self.values
    }
    pub fn outcomes_sha256(&self) -> &str {
        &self.outcomes_sha256
    }
    pub(crate) fn verify_panel(&self, panel: &MatchedPanel) -> Result<()> {
        ensure!(
            self.panel_metadata_sha256 == panel.metadata_sha256,
            "outcomes belong to different panel metadata"
        );
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct GroupSplit {
    panel_metadata_sha256: String,
    training_sequence_indices: Vec<usize>,
    test_sequence_indices: Vec<usize>,
    pub training_sequence_ids: Vec<String>,
    pub test_sequence_ids: Vec<String>,
    pub training_groups: Vec<String>,
    pub test_groups: Vec<String>,
    /// Optional outcome-free physical-neighborhood exclusions, removed in every context.
    pub excluded_guard_sequence_ids: Vec<String>,
    pub contexts: Vec<String>,
    pub shared_groups: usize,
    pub shared_exact_or_reverse_complement_dna: usize,
}

impl GroupSplit {
    /// Hold out whole declared groups, with identical sequence membership in every context.
    /// Similarity/physical-distance policy belongs in the outcome-free group assignment.
    pub fn new(panel: &MatchedPanel, test_groups: &[String]) -> Result<Self> {
        Self::with_guard(panel, test_groups, &[])
    }

    /// Guard membership must be derived from metadata by the caller's frozen rule.
    /// Test groups remain intact; training sequences can be additionally excluded.
    pub fn with_guard(
        panel: &MatchedPanel,
        test_groups: &[String],
        guarded_sequence_ids: &[String],
    ) -> Result<Self> {
        let test: BTreeSet<_> = test_groups.iter().cloned().collect();
        ensure!(test.len() == test_groups.len(), "duplicate test group");
        let all: BTreeSet<_> = panel.groups().into_iter().collect();
        ensure!(
            !test.is_empty() && test.len() < all.len(),
            "split needs nonempty training and test groups"
        );
        ensure!(test.is_subset(&all), "unknown test group");
        let guarded: BTreeSet<_> = guarded_sequence_ids.iter().cloned().collect();
        ensure!(
            guarded.len() == guarded_sequence_ids.len(),
            "duplicate guard sequence ID"
        );
        for id in &guarded {
            let index = panel.sequence_index(id)?;
            ensure!(
                !test.contains(&panel.sequences[index].group),
                "test sequence cannot be a training guard exclusion: {id}"
            );
        }
        let mut training_sequence_indices = Vec::new();
        let mut test_sequence_indices = Vec::new();
        for (i, sequence) in panel.sequences.iter().enumerate() {
            if test.contains(&sequence.group) {
                test_sequence_indices.push(i);
            } else if !guarded.contains(&sequence.id) {
                training_sequence_indices.push(i);
            }
        }
        ensure!(
            !training_sequence_indices.is_empty(),
            "guard removes all training sequences"
        );
        let training_groups = training_sequence_indices
            .iter()
            .map(|index| panel.sequences[*index].group.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let training_sequence_ids = training_sequence_indices
            .iter()
            .map(|i| panel.sequences[*i].id.clone())
            .collect();
        let test_sequence_ids = test_sequence_indices
            .iter()
            .map(|i| panel.sequences[*i].id.clone())
            .collect();
        Ok(Self {
            panel_metadata_sha256: panel.metadata_sha256.clone(),
            training_sequence_indices,
            test_sequence_indices,
            training_sequence_ids,
            test_sequence_ids,
            training_groups,
            test_groups: test.into_iter().collect(),
            excluded_guard_sequence_ids: guarded.into_iter().collect(),
            contexts: panel.contexts.clone(),
            shared_groups: 0,
            shared_exact_or_reverse_complement_dna: 0,
        })
    }

    /// Every sequence group must occur in exactly one test fold.
    pub fn partition(panel: &MatchedPanel, fold_test_groups: &[Vec<String>]) -> Result<Vec<Self>> {
        ensure!(
            fold_test_groups.len() >= 2,
            "cross-validation needs at least two folds"
        );
        let mut seen = BTreeSet::new();
        for groups in fold_test_groups {
            for group in groups {
                ensure!(
                    seen.insert(group.clone()),
                    "group occurs in multiple test folds: {group}"
                );
            }
        }
        ensure!(
            seen == panel.groups().into_iter().collect(),
            "test-fold groups do not partition the panel"
        );
        fold_test_groups
            .iter()
            .map(|groups| Self::new(panel, groups))
            .collect()
    }

    pub fn training_indices(&self) -> &[usize] {
        &self.training_sequence_indices
    }
    pub fn test_indices(&self) -> &[usize] {
        &self.test_sequence_indices
    }
    pub(crate) fn verify_panel(&self, panel: &MatchedPanel) -> Result<()> {
        ensure!(
            self.panel_metadata_sha256 == panel.metadata_sha256,
            "split belongs to different panel metadata"
        );
        Ok(())
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    pub(crate) fn panel() -> MatchedPanel {
        MatchedPanel::new(
            vec![
                SequenceRecord {
                    id: "a".into(),
                    group: "g1".into(),
                    sequence: "AAAAAA".into(),
                },
                SequenceRecord {
                    id: "b".into(),
                    group: "g2".into(),
                    sequence: "CCCCCC".into(),
                },
                SequenceRecord {
                    id: "c".into(),
                    group: "g3".into(),
                    sequence: "ACACAC".into(),
                },
            ],
            vec!["ctx2".into(), "ctx1".into()],
        )
        .unwrap()
    }
    pub(crate) fn outcomes(panel: &MatchedPanel) -> Vec<MeasuredOutcome> {
        panel
            .sequences()
            .iter()
            .enumerate()
            .flat_map(|(i, sequence)| {
                panel
                    .contexts()
                    .iter()
                    .enumerate()
                    .map(move |(j, context)| MeasuredOutcome {
                        sequence_id: sequence.id.clone(),
                        context: context.clone(),
                        value: i as f64 + 10.0 * j as f64,
                    })
            })
            .collect()
    }

    #[test]
    fn complete_pairing_and_source_identity_are_enforced() {
        let panel = panel();
        let rows = outcomes(&panel);
        let matrix = OutcomeMatrix::from_records(&panel, &rows).unwrap();
        assert_eq!(
            matrix.values(),
            &[vec![0.0, 10.0], vec![1.0, 11.0], vec![2.0, 12.0]]
        );
        assert!(OutcomeMatrix::from_records(&panel, &rows[..rows.len() - 1]).is_err());
        let mut duplicates = rows.clone();
        duplicates.push(rows[0].clone());
        assert!(OutcomeMatrix::from_records(&panel, &duplicates).is_err());
        let mut nonfinite = rows;
        nonfinite[0].value = f64::NAN;
        assert!(OutcomeMatrix::from_records(&panel, &nonfinite).is_err());
    }

    #[test]
    fn exact_and_rc_identity_cannot_cross_groups() {
        let mut rows = panel().sequences().to_vec();
        rows.push(SequenceRecord {
            id: "duplicate".into(),
            group: "different".into(),
            sequence: "TTTTTT".into(),
        });
        assert!(MatchedPanel::new(rows.clone(), vec!["ctx".into()]).is_err());
        rows.last_mut().unwrap().group = "g1".into();
        assert!(MatchedPanel::new(rows, vec!["ctx".into()]).is_ok());
    }

    #[test]
    fn folds_partition_whole_groups_without_outcomes() {
        let panel = panel();
        let folds = GroupSplit::partition(
            &panel,
            &[vec!["g1".into()], vec!["g2".into()], vec!["g3".into()]],
        )
        .unwrap();
        assert!(folds
            .iter()
            .all(|split| split.training_indices().len() == 2 && split.test_indices().len() == 1));
        assert_eq!(folds[0].test_sequence_ids, vec!["a"]);
        assert_eq!(folds[0].contexts, panel.contexts());
        assert!(GroupSplit::partition(
            &panel,
            &[
                vec!["g1".into()],
                vec!["g1".into(), "g2".into(), "g3".into()]
            ]
        )
        .is_err());
        assert!(GroupSplit::partition(&panel, &[vec!["g1".into()], vec!["g2".into()]]).is_err());
        assert!(GroupSplit::new(&panel, &["unknown".into()]).is_err());
    }

    #[test]
    fn metadata_guard_excludes_a_sequence_in_every_context() {
        let panel = panel();
        let split = GroupSplit::with_guard(&panel, &["g3".into()], &["b".into()]).unwrap();
        assert_eq!(split.training_sequence_ids, vec!["a"]);
        assert_eq!(split.test_sequence_ids, vec!["c"]);
        assert_eq!(split.excluded_guard_sequence_ids, vec!["b"]);
        assert_eq!(split.training_groups, vec!["g1"]);
        assert!(GroupSplit::with_guard(&panel, &["g3".into()], &["c".into()]).is_err());
        assert!(GroupSplit::with_guard(&panel, &["g3".into()], &["a".into(), "b".into()]).is_err());
    }
}
