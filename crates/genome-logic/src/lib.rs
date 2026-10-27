//! Context-sensitive contradictions and finite propositional implication chains.
use anyhow::{ensure, Result};
use genome_knowledge::EvidenceStore;
use genome_types::Probability;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Serialize)]
pub struct Contradiction {
    pub left_evidence: String,
    pub right_evidence: String,
    pub variant_key: String,
}
pub fn contradictions(store: &EvidenceStore) -> Vec<Contradiction> {
    let records: Vec<_> = store.records().iter().collect();
    let mut result = Vec::new();
    for (i, (left_id, left)) in records.iter().enumerate() {
        for (right_id, right) in &records[i + 1..] {
            let a = &left.proposition;
            let b = &right.proposition;
            if a.variant_key == b.variant_key
                && a.context == b.context
                && a.classification != b.classification
            {
                result.push(Contradiction {
                    left_evidence: (*left_id).clone(),
                    right_evidence: (*right_id).clone(),
                    variant_key: a.variant_key.clone(),
                });
            }
        }
    }
    result
}

#[derive(Clone, Debug, Serialize)]
pub struct Rule {
    pub premise: String,
    pub conclusion: String,
    pub assumption: String,
}
#[derive(Debug, Serialize)]
pub struct Derivation {
    assumptions: Vec<String>,
    rules: Vec<Rule>,
    conclusion: String,
    proof_expression: String,
}
impl Derivation {
    pub fn derive(assumptions: Vec<String>, rules: Vec<Rule>, conclusion: String) -> Result<Self> {
        ensure!(
            !assumptions.is_empty() && assumptions.len() <= 128 && rules.len() <= 128,
            "proof budget exceeded or no premises"
        );
        let mut known: BTreeMap<String, String> = BTreeMap::new();
        for (index, atom) in assumptions.iter().enumerate() {
            ensure!(!atom.trim().is_empty(), "empty premise");
            known
                .entry(atom.clone())
                .or_insert_with(|| format!("h{index}"));
        }
        for _ in 0..=rules.len() {
            let mut changed = false;
            for (index, rule) in rules.iter().enumerate() {
                ensure!(
                    !rule.assumption.trim().is_empty(),
                    "implications need explicit scientific assumptions"
                );
                if !known.contains_key(&rule.conclusion) {
                    if let Some(expression) = known.get(&rule.premise).cloned() {
                        known.insert(rule.conclusion.clone(), format!("(r{index} {expression})"));
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
        }
        let proof_expression = known
            .get(&conclusion)
            .cloned()
            .ok_or_else(|| anyhow::anyhow!("conclusion is not entailed by supplied assumptions"))?;
        Ok(Self {
            assumptions,
            rules,
            conclusion,
            proof_expression,
        })
    }
    pub fn lean(&self) -> String {
        let atoms: BTreeSet<_> = self
            .assumptions
            .iter()
            .chain(self.rules.iter().flat_map(|r| [&r.premise, &r.conclusion]))
            .chain([&self.conclusion])
            .collect();
        let names: BTreeMap<_, _> = atoms
            .into_iter()
            .enumerate()
            .map(|(i, atom)| (atom, format!("p{i}")))
            .collect();
        let mut text = String::from("theorem conditional_conclusion ");
        for name in names.values() {
            text.push_str(&format!("({name} : Prop) "));
        }
        for (i, atom) in self.assumptions.iter().enumerate() {
            text.push_str(&format!("(h{i} : {}) ", names[atom]));
        }
        for (i, rule) in self.rules.iter().enumerate() {
            text.push_str(&format!(
                "(r{i} : {} → {}) ",
                names[&rule.premise], names[&rule.conclusion]
            ));
        }
        text.push_str(&format!(
            ": {} := {}\n",
            names[&self.conclusion], self.proof_expression
        ));
        text
    }
}
pub fn probability_obligations(probabilities: &[Probability]) -> String {
    let mut text = String::new();
    for (i, p) in probabilities.iter().enumerate() {
        text.push_str(&format!(
            "theorem probability_{i} : ({} : Nat) ≤ {} ∧ (0 : Nat) < {} := by decide\n",
            p.numerator(),
            p.denominator(),
            p.denominator()
        ));
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;
    use genome_knowledge::{Evidence, EvidenceKind, Proposition};
    use genome_types::{BiologicalContext, ClinicalLabel};
    #[test]
    fn cannot_invent_conclusions_or_escape_via_atom_text() {
        let rule = Rule {
            premise: "a".into(),
            conclusion: "b".into(),
            assumption: "conditional relation".into(),
        };
        assert!(Derivation::derive(vec!["c".into()], vec![rule.clone()], "b".into()).is_err());
        let d = Derivation::derive(vec!["a".into()], vec![rule], "b".into()).unwrap();
        assert!(d.lean().contains("(r0 h0)"));
        let d = Derivation::derive(
            vec!["evil\naxiom injected : False".into()],
            vec![],
            "evil\naxiom injected : False".into(),
        )
        .unwrap();
        assert!(!d.lean().contains("axiom"));
    }
    #[test]
    fn opposite_labels_only_conflict_in_identical_context() {
        let mut store = EvidenceStore::default();
        for (label, gene) in [
            (ClinicalLabel::Benign, "HBB"),
            (ClinicalLabel::Pathogenic, "HBB"),
            (ClinicalLabel::Pathogenic, "CFTR"),
        ] {
            let mut context = BiologicalContext::hbb_germline();
            context.gene = gene.into();
            store
                .insert(Evidence {
                    proposition: Proposition {
                        variant_key: "v".into(),
                        context,
                        classification: label,
                    },
                    kind: EvidenceKind::CuratedClinical,
                    source_url: "https://example.org".into(),
                    source_sha256: "a".repeat(64),
                    accession: format!("{label:?}-{gene}"),
                    study_id: None,
                })
                .unwrap();
        }
        assert_eq!(contradictions(&store).len(), 1);
    }
}
