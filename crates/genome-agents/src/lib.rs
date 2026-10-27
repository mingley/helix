//! Bounded offline research agent with an adapter boundary for future LLM proposals.
use anyhow::{ensure, Result};
use genome_data::TrainingExample;
use genome_inference::{ConsequenceModel, Prediction, Predictor};
use genome_knowledge::{Claim, Evidence, EvidenceKind, EvidenceStore, Ledger, Proposition};
use genome_logic::{contradictions, Contradiction, Derivation, Rule};
use genome_types::{Consequence, FeatureRecord};
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct ResearchQuestion {
    pub consequence: Consequence,
    pub question: String,
    pub proposed_test: String,
    pub supporting_evidence: Vec<String>,
}
pub struct ResearchInput<'a> {
    pub training: &'a [TrainingExample],
    pub targets: &'a [FeatureRecord],
    pub source_sha256: &'a str,
    pub max_steps: usize,
}
pub struct ResearchResult {
    pub predictions: Vec<Prediction>,
    pub questions: Vec<ResearchQuestion>,
    pub evidence: EvidenceStore,
    pub contradictions: Vec<Contradiction>,
    pub derivation: Derivation,
    pub model: ConsequenceModel,
    pub ledger: Ledger,
}
pub trait ResearchAgent {
    fn run(&self, input: ResearchInput<'_>) -> Result<ResearchResult>;
}
/// Deterministic proposals are the reproducible baseline, not an LLM or a claim
/// of autonomous biological discovery. Agents never get held-out annotations.
pub struct OfflineAgent;
impl ResearchAgent for OfflineAgent {
    fn run(&self, input: ResearchInput<'_>) -> Result<ResearchResult> {
        let required = input
            .training
            .len()
            .checked_add(input.targets.len())
            .ok_or_else(|| anyhow::anyhow!("agent budget overflow"))?;
        ensure!(
            required <= input.max_steps && input.max_steps <= 100_000,
            "agent step budget exhausted or exceeds hard limit"
        );
        let mut evidence = EvidenceStore::default();
        let mut ledger = Ledger::default();
        ledger.append(
            "discover",
            serde_json::json!({
                "objective": "Predict HBB binary ClinVar classification from consequence groups",
                "scope": "Public computational demonstration; no clinical recommendation",
                "training_count": input.training.len(), "target_count": input.targets.len()
            }),
        )?;
        let mut by_consequence = std::collections::BTreeMap::<Consequence, Vec<String>>::new();
        for example in input.training {
            let proposition = Proposition {
                variant_key: example.feature.variant.key(),
                context: example.feature.context.clone(),
                classification: example.observation.label,
            };
            let id = evidence.insert(Evidence {
                proposition: proposition.clone(),
                kind: EvidenceKind::CuratedClinical,
                source_url: example.observation.source_url.clone(),
                source_sha256: input.source_sha256.into(),
                accession: example.observation.id.clone(),
                study_id: None,
            })?;
            let supported = Claim::new(proposition).support(&evidence, vec![id.clone()])?;
            ledger.append(
                "supported_curated_assertion",
                serde_json::json!({
                    "proposition": supported.proposition(), "evidence": supported.evidence(),
                    "replicated": false, "empirically_measured": false
                }),
            )?;
            by_consequence
                .entry(example.feature.consequence)
                .or_default()
                .push(id);
        }
        let contradictions = contradictions(&evidence);
        ledger.append("investigate_contradictions", &contradictions)?;
        let model = ConsequenceModel::fit(input.training)?;
        let mut questions = Vec::new();
        for (consequence, counts) in &model.groups {
            if counts.total < 3 || (counts.pathogenic > 0 && counts.pathogenic < counts.total) {
                questions.push(ResearchQuestion {
                    consequence: *consequence,
                    question: format!("Does {consequence:?} separate HBB clinical classifications across independent studies?"),
                    proposed_test: "Freeze the model; compare calibration against a new independently curated cohort and a functional assay. Mixed labels are model gaps, not automatically contradictory measurements.".into(),
                    supporting_evidence: by_consequence[consequence].clone(),
                });
            }
        }
        ledger.append("formalize_testable_questions", &questions)?;
        let derivation = Derivation::derive(
            vec!["Training assertions and molecular consequence annotations are reliable for this context".into()],
            vec![
                Rule {
                    premise: "Training assertions and molecular consequence annotations are reliable for this context".into(),
                    conclusion: "Consequence-group frequencies estimate target classification probabilities".into(),
                    assumption: "Target records are exchangeable with training records within consequence groups; Beta(1,1) smoothing is appropriate".into(),
                },
                Rule {
                    premise: "Consequence-group frequencies estimate target classification probabilities".into(),
                    conclusion: "Thresholding at probability >= 0.5 yields the stated binary predictions under the assumed decision loss".into(),
                    assumption: "Equal false-positive and false-negative costs; clinical utility is not established".into(),
                }
            ],
            "Thresholding at probability >= 0.5 yields the stated binary predictions under the assumed decision loss".into(),
        )?;
        ledger.append("conditional_reasoning", &derivation)?;
        let predictions: Vec<_> = input
            .targets
            .iter()
            .map(|f| model.predict(f))
            .collect::<Result<_>>()?;
        ledger.append("freeze_predictions_before_evaluation", &predictions)?;
        Ok(ResearchResult {
            predictions,
            questions,
            evidence,
            contradictions,
            derivation,
            model,
            ledger,
        })
    }
}
