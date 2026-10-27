//! Evidence-backed claims and a deterministic, tamper-evident revision ledger.
use anyhow::{ensure, Context, Result};
use genome_types::{BiologicalContext, ClinicalLabel};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    marker::PhantomData,
};

fn digest(value: &impl Serialize) -> Result<String> {
    Ok(format!("{:x}", Sha256::digest(serde_json::to_vec(value)?)))
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Proposition {
    pub variant_key: String,
    pub context: BiologicalContext,
    pub classification: ClinicalLabel,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EvidenceKind {
    CuratedClinical,
    Measured,
    Predicted,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Evidence {
    pub proposition: Proposition,
    pub kind: EvidenceKind,
    pub source_url: String,
    pub source_sha256: String,
    pub accession: String,
    /// Independent experiment/cohort identity; database accession is not independence.
    pub study_id: Option<String>,
}
#[derive(Default, Debug, Serialize)]
pub struct EvidenceStore {
    records: BTreeMap<String, Evidence>,
}
impl EvidenceStore {
    pub fn insert(&mut self, evidence: Evidence) -> Result<String> {
        ensure!(
            evidence.source_url.starts_with("https://"),
            "evidence needs an HTTPS source"
        );
        ensure!(
            evidence.source_sha256.len() == 64
                && evidence
                    .source_sha256
                    .bytes()
                    .all(|b| b.is_ascii_hexdigit()),
            "evidence needs a SHA-256 source commitment"
        );
        ensure!(
            !evidence.accession.trim().is_empty(),
            "evidence accession is required"
        );
        let id = digest(&evidence)?;
        self.records.insert(id.clone(), evidence);
        Ok(id)
    }
    pub fn get(&self, id: &str) -> Option<&Evidence> {
        self.records.get(id)
    }
    pub fn records(&self) -> &BTreeMap<String, Evidence> {
        &self.records
    }
}

pub struct Hypothesis;
pub struct Supported;
pub struct Replicated;
pub struct Claim<State> {
    proposition: Proposition,
    evidence: Vec<String>,
    state: PhantomData<State>,
}
impl Claim<Hypothesis> {
    pub fn new(proposition: Proposition) -> Self {
        Self {
            proposition,
            evidence: Vec::new(),
            state: PhantomData,
        }
    }
    pub fn support(self, store: &EvidenceStore, evidence: Vec<String>) -> Result<Claim<Supported>> {
        ensure!(!evidence.is_empty(), "support requires evidence");
        let mut seen = BTreeSet::new();
        for id in &evidence {
            ensure!(
                seen.insert(id),
                "duplicate evidence is not additional support"
            );
            let record = store.get(id).context("unknown evidence reference")?;
            ensure!(
                record.kind != EvidenceKind::Predicted,
                "a prediction cannot support itself"
            );
            ensure!(
                record.proposition == self.proposition,
                "evidence has a different variant, context, or conclusion"
            );
        }
        Ok(Claim {
            proposition: self.proposition,
            evidence,
            state: PhantomData,
        })
    }
}
impl Claim<Supported> {
    pub fn replicate(self, store: &EvidenceStore) -> Result<Claim<Replicated>> {
        let studies: BTreeSet<_> = self
            .evidence
            .iter()
            .filter_map(|id| store.get(id))
            .filter(|e| e.kind == EvidenceKind::Measured)
            .filter_map(|e| e.study_id.as_deref())
            .filter(|s| !s.trim().is_empty())
            .collect();
        ensure!(
            studies.len() >= 2,
            "replication needs measurements from at least two independent studies"
        );
        Ok(Claim {
            proposition: self.proposition,
            evidence: self.evidence,
            state: PhantomData,
        })
    }
}
impl<S> Claim<S> {
    pub fn proposition(&self) -> &Proposition {
        &self.proposition
    }
    pub fn evidence(&self) -> &[String] {
        &self.evidence
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Revision {
    pub sequence: usize,
    pub previous: Option<String>,
    pub kind: String,
    pub payload: serde_json::Value,
    pub hash: String,
}
#[derive(Default, Debug, Serialize)]
pub struct Ledger {
    revisions: Vec<Revision>,
}
#[derive(Serialize)]
struct Commitment<'a> {
    sequence: usize,
    previous: &'a Option<String>,
    kind: &'a str,
    payload: &'a serde_json::Value,
}
impl Ledger {
    pub fn append(&mut self, kind: &str, payload: impl Serialize) -> Result<()> {
        let sequence = self.revisions.len();
        let previous = self.revisions.last().map(|r| r.hash.clone());
        let payload = serde_json::to_value(payload)?;
        let hash = digest(&Commitment {
            sequence,
            previous: &previous,
            kind,
            payload: &payload,
        })?;
        self.revisions.push(Revision {
            sequence,
            previous,
            kind: kind.into(),
            payload,
            hash,
        });
        Ok(())
    }
    pub fn revisions(&self) -> &[Revision] {
        &self.revisions
    }
    pub fn verify(revisions: &[Revision]) -> Result<()> {
        ensure!(!revisions.is_empty(), "empty ledger");
        let mut previous = None;
        for (sequence, revision) in revisions.iter().enumerate() {
            ensure!(
                revision.sequence == sequence && revision.previous == previous,
                "broken revision chain"
            );
            let hash = digest(&Commitment {
                sequence,
                previous: &previous,
                kind: &revision.kind,
                payload: &revision.payload,
            })?;
            ensure!(
                hash == revision.hash,
                "revision payload/hash mismatch at revision {sequence} ({})",
                revision.kind
            );
            previous = Some(hash);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn proposition(label: ClinicalLabel) -> Proposition {
        Proposition {
            variant_key: "GRCh38:11:1:A>T".into(),
            context: BiologicalContext::hbb_germline(),
            classification: label,
        }
    }
    #[test]
    fn predictions_and_wrong_context_cannot_be_promoted() {
        let p = proposition(ClinicalLabel::Pathogenic);
        let mut store = EvidenceStore::default();
        let mut evidence = Evidence {
            proposition: p.clone(),
            kind: EvidenceKind::Predicted,
            source_url: "https://example.org".into(),
            source_sha256: "a".repeat(64),
            accession: "test".into(),
            study_id: None,
        };
        let id = store.insert(evidence.clone()).unwrap();
        assert!(Claim::new(p.clone()).support(&store, vec![id]).is_err());
        evidence.kind = EvidenceKind::CuratedClinical;
        evidence.proposition.context.gene = "CFTR".into();
        let id = store.insert(evidence).unwrap();
        assert!(Claim::new(p).support(&store, vec![id]).is_err());
    }
    #[test]
    fn curated_database_records_do_not_establish_replication() {
        let p = proposition(ClinicalLabel::Benign);
        let mut store = EvidenceStore::default();
        let ids = ["a", "b"].map(|study| {
            store
                .insert(Evidence {
                    proposition: p.clone(),
                    kind: EvidenceKind::CuratedClinical,
                    source_url: "https://example.org".into(),
                    source_sha256: "a".repeat(64),
                    accession: study.into(),
                    study_id: Some(study.into()),
                })
                .unwrap()
        });
        let claim = Claim::new(p).support(&store, ids.into()).unwrap();
        assert!(claim.replicate(&store).is_err());
    }
    #[test]
    fn ledger_detects_modification_and_reordering() {
        let mut ledger = Ledger::default();
        ledger.append("observe", "one").unwrap();
        ledger.append("update", "two").unwrap();
        Ledger::verify(ledger.revisions()).unwrap();
        let mut tampered = ledger.revisions().to_vec();
        tampered[0].payload = serde_json::json!("different");
        assert!(Ledger::verify(&tampered).is_err());
        tampered.swap(0, 1);
        assert!(Ledger::verify(&tampered).is_err());
    }
    #[test]
    fn numeric_evaluation_commitments_survive_disk_round_trip() {
        let mut ledger = Ledger::default();
        let metrics: Vec<f64> = (1..1000).map(|n| (n as f64 / 1003.0).ln()).collect();
        ledger.append("evaluation", metrics).unwrap();
        let saved = serde_json::to_vec_pretty(ledger.revisions()).unwrap();
        let restored: Vec<Revision> = serde_json::from_slice(&saved).unwrap();
        Ledger::verify(&restored).unwrap();
    }
}
