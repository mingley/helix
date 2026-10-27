//! Complete-panel models with a common 33-feature, training-only scaler.
use crate::{
    data::{GroupSplit, MatchedPanel, OutcomeMatrix},
    features::{feature_names, SequenceFeatures},
};
use anyhow::{ensure, Result};
use genome_regulatory::math::{RidgeDesign, StandardizedRidge};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ModelKind {
    GlobalMean,
    ContextMean,
    SharedDna,
    AdditiveContextDna,
    IndependentContextDna,
}

impl ModelKind {
    pub const ALL: [Self; 5] = [
        Self::GlobalMean,
        Self::ContextMean,
        Self::SharedDna,
        Self::AdditiveContextDna,
        Self::IndependentContextDna,
    ];
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub struct ModelPredictions {
    pub global_mean: f64,
    pub context_mean: f64,
    pub shared_dna: f64,
    pub additive_context_dna: f64,
    pub independent_context_dna: f64,
}

impl ModelPredictions {
    pub fn get(&self, model: ModelKind) -> f64 {
        match model {
            ModelKind::GlobalMean => self.global_mean,
            ModelKind::ContextMean => self.context_mean,
            ModelKind::SharedDna => self.shared_dna,
            ModelKind::AdditiveContextDna => self.additive_context_dna,
            ModelKind::IndependentContextDna => self.independent_context_dna,
        }
    }
    pub fn values(&self) -> [f64; 5] {
        ModelKind::ALL.map(|model| self.get(model))
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ContextFit {
    pub panel_metadata_sha256: String,
    /// Depends exclusively on training rows; the full outcome digest is reported by evaluation.
    pub training_outcomes_sha256: String,
    pub training_sequence_ids: Vec<String>,
    pub contexts: Vec<String>,
    pub training_sequences: usize,
    pub training_measurements: usize,
    pub alpha: f64,
    pub feature_names: Vec<String>,
    pub global_mean: f64,
    pub training_context_means: BTreeMap<String, f64>,
    pub shared_dna_fit: StandardizedRidge,
    pub independent_context_fits: BTreeMap<String, StandardizedRidge>,
    pub objective: String,
}

impl ContextFit {
    /// Shared DNA fits one context-mean response per sequence using sum SSE +
    /// alpha||beta||². Each independent context uses the same n training sequences
    /// and alpha. The equivalent shared pooled objective has penalty C*alpha.
    ///
    /// Complete pairing implies the additive model shares S's slopes exactly;
    /// its free context offsets are training_context_mean - global_mean.
    pub fn fit(
        panel: &MatchedPanel,
        split: &GroupSplit,
        outcomes: &OutcomeMatrix,
        alpha: f64,
    ) -> Result<Self> {
        split.verify_panel(panel)?;
        outcomes.verify_panel(panel)?;
        ensure!(
            alpha.is_finite() && alpha > 0.0,
            "alpha must be finite and positive"
        );
        let indices = split.training_indices();
        let features: Vec<_> = indices
            .iter()
            .map(|index| panel.features()[*index].vector())
            .collect();
        // The same outcome-independent factorization/scaler is reused by every fit.
        let design = RidgeDesign::prepare(&features, alpha)?;
        let shared_outcomes: Vec<f64> = indices
            .iter()
            .map(|index| finite_mean(&outcomes.values()[*index]))
            .collect::<Result<_>>()?;
        let shared_dna_fit = design.fit(&shared_outcomes)?;
        let global_mean = shared_dna_fit.intercept;
        let mut training_context_means = BTreeMap::new();
        let mut independent_context_fits = BTreeMap::new();
        for (j, context) in panel.contexts().iter().enumerate() {
            let values: Vec<_> = indices
                .iter()
                .map(|index| outcomes.values()[*index][j])
                .collect();
            training_context_means.insert(context.clone(), finite_mean(&values)?);
            independent_context_fits.insert(context.clone(), design.fit(&values)?);
        }
        let mut digest = Sha256::new();
        digest.update(b"context-training-outcomes-f64-v1");
        digest.update(panel.metadata_sha256().as_bytes());
        for index in indices {
            digest.update((*index as u64).to_le_bytes());
            for value in &outcomes.values()[*index] {
                digest.update(value.to_bits().to_le_bytes());
            }
        }
        Ok(Self {
            panel_metadata_sha256: panel.metadata_sha256().to_owned(),
            training_outcomes_sha256: format!("{:x}",digest.finalize()),
            training_sequence_ids: split.training_sequence_ids.clone(), contexts: panel.contexts().to_vec(),
            training_sequences: indices.len(), training_measurements: indices.len()*panel.contexts().len(),
            alpha, feature_names: feature_names(), global_mean, training_context_means, shared_dna_fit,
            independent_context_fits,
            objective: "S: sum over unique training sequences of squared error against the arithmetic mean across contexts + alpha*||beta||²; I: same n sequences and alpha independently per context; A: S plus unpenalized training context offsets; population-SD training-only scaling shared across S/A/I; free intercepts".into(),
        })
    }

    pub fn predict(&self, features: &SequenceFeatures, context: &str) -> Result<ModelPredictions> {
        let context_mean = *self
            .training_context_means
            .get(context)
            .ok_or_else(|| anyhow::anyhow!("context absent from training panel: {context}"))?;
        let row = features.vector();
        ensure!(
            row.len() == self.shared_dna_fit.feature_count
                && row.iter().all(|value| value.is_finite()),
            "invalid context DNA feature row"
        );
        let shared_sequence_effect: f64 = row
            .iter()
            .enumerate()
            .map(|(j, x)| {
                (x - self.shared_dna_fit.training_feature_means[j])
                    / self.shared_dna_fit.training_feature_scales[j]
                    * self.shared_dna_fit.standardized_coefficients[j]
            })
            .sum();
        let shared_dna = self.global_mean + shared_sequence_effect;
        let independent_context_dna = self.independent_context_fits[context].predict(&row)?;
        let predictions = ModelPredictions {
            global_mean: self.global_mean,
            context_mean,
            shared_dna,
            additive_context_dna: context_mean + shared_sequence_effect,
            independent_context_dna,
        };
        ensure!(
            predictions.values().iter().all(|x| x.is_finite()),
            "context prediction overflow"
        );
        Ok(predictions)
    }

    pub fn predict_sequence(
        &self,
        panel: &MatchedPanel,
        id: &str,
        context: &str,
    ) -> Result<ModelPredictions> {
        ensure!(
            self.panel_metadata_sha256 == panel.metadata_sha256(),
            "prediction panel metadata mismatch"
        );
        self.predict(&panel.features()[panel.sequence_index(id)?], context)
    }
}

pub(crate) fn finite_mean(values: &[f64]) -> Result<f64> {
    ensure!(
        !values.is_empty() && values.iter().all(|value| value.is_finite()),
        "mean needs nonempty finite values"
    );
    let mean = if values.iter().all(|value| *value == values[0]) {
        values[0]
    } else {
        values.iter().map(|value| value / values.len() as f64).sum()
    };
    ensure!(mean.is_finite(), "mean overflow");
    Ok(mean)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{
        tests::{outcomes, panel},
        MeasuredOutcome,
    };

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-12, "{a} != {b}");
    }

    #[test]
    fn shared_additive_and_independent_slopes_obey_balanced_panel_algebra() {
        let panel = panel();
        let matrix = OutcomeMatrix::from_records(&panel, &outcomes(&panel)).unwrap();
        let split = GroupSplit::new(&panel, &["g3".into()]).unwrap();
        let fit = ContextFit::fit(&panel, &split, &matrix, 10.0).unwrap();
        close(fit.global_mean, 5.5);
        close(fit.training_context_means["ctx1"], 0.5);
        close(fit.training_context_means["ctx2"], 10.5);
        for sequence in panel.sequences() {
            let one = fit.predict_sequence(&panel, &sequence.id, "ctx1").unwrap();
            let two = fit.predict_sequence(&panel, &sequence.id, "ctx2").unwrap();
            close(one.additive_context_dna - one.shared_dna, -5.0);
            close(two.additive_context_dna - two.shared_dna, 5.0);
            close(one.independent_context_dna, one.additive_context_dna);
            close(two.independent_context_dna, two.additive_context_dna);
            close(
                (one.independent_context_dna + two.independent_context_dna) / 2.0,
                one.shared_dna,
            );
        }
        for j in 0..33 {
            close(
                fit.shared_dna_fit.standardized_coefficients[j],
                fit.independent_context_fits
                    .values()
                    .map(|context| context.standardized_coefficients[j] / 2.0)
                    .sum(),
            );
        }
        for context in fit.independent_context_fits.values() {
            assert_eq!(
                context.training_feature_means,
                fit.shared_dna_fit.training_feature_means
            );
            assert_eq!(
                context.training_feature_scales,
                fit.shared_dna_fit.training_feature_scales
            );
        }
    }

    #[test]
    fn changing_held_out_outcomes_cannot_change_fits_or_predictions() {
        let panel = panel();
        let split = GroupSplit::new(&panel, &["g3".into()]).unwrap();
        let rows = outcomes(&panel);
        let matrix = OutcomeMatrix::from_records(&panel, &rows).unwrap();
        let fit = ContextFit::fit(&panel, &split, &matrix, 10.0).unwrap();
        let changed: Vec<_> = rows
            .into_iter()
            .map(|mut row| {
                if row.sequence_id == "c" {
                    row.value = 1e9;
                }
                row
            })
            .collect();
        let changed_matrix = OutcomeMatrix::from_records(&panel, &changed).unwrap();
        let changed_fit = ContextFit::fit(&panel, &split, &changed_matrix, 10.0).unwrap();
        assert_ne!(matrix.outcomes_sha256(), changed_matrix.outcomes_sha256());
        assert_eq!(
            fit.training_outcomes_sha256,
            changed_fit.training_outcomes_sha256
        );
        assert_eq!(
            fit.shared_dna_fit.training_feature_means,
            changed_fit.shared_dna_fit.training_feature_means
        );
        for context in panel.contexts() {
            assert_eq!(
                fit.predict_sequence(&panel, "c", context).unwrap(),
                changed_fit.predict_sequence(&panel, "c", context).unwrap()
            );
        }
        // Only AAA and CCC training rows contribute to the GC scaler.
        close(fit.shared_dna_fit.training_feature_means[32], 0.5);
        close(fit.shared_dna_fit.training_feature_scales[32], 0.5);
    }

    #[test]
    fn constant_outcomes_and_unknown_context_are_explicit() {
        let panel = panel();
        let rows: Vec<MeasuredOutcome> = outcomes(&panel)
            .into_iter()
            .map(|mut row| {
                row.value = 0.1;
                row
            })
            .collect();
        let matrix = OutcomeMatrix::from_records(&panel, &rows).unwrap();
        let split = GroupSplit::new(&panel, &["g3".into()]).unwrap();
        let fit = ContextFit::fit(&panel, &split, &matrix, 10.0).unwrap();
        assert_eq!(
            fit.predict_sequence(&panel, "c", "ctx1").unwrap().values(),
            [0.1; 5]
        );
        assert!(fit.predict_sequence(&panel, "c", "unknown").is_err());
        assert!(ContextFit::fit(&panel, &split, &matrix, 0.0).is_err());
    }

    #[test]
    fn one_feature_closed_form_uses_sum_sse_ridge_penalty() {
        // This verifies the reused solver against beta=X'y/(X'X+alpha).
        let fit = StandardizedRidge::fit(&[vec![0.0], vec![2.0]], &[0.0, 2.0], 2.0).unwrap();
        close(fit.training_feature_means[0], 1.0);
        close(fit.training_feature_scales[0], 1.0);
        close(fit.standardized_coefficients[0], 0.5);
        close(fit.predict(&[0.0]).unwrap(), 0.5);
        close(fit.predict(&[2.0]).unwrap(), 1.5);
    }
}
