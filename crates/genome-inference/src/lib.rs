//! Auditable, deterministic model adapters. No held-out labels enter fitting/prediction.
use anyhow::{ensure, Result};
use genome_data::TrainingExample;
use genome_types::{Consequence, FeatureRecord, Probability};
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Clone, Debug, Serialize)]
pub struct Prediction {
    pub id: String,
    pub variant_key: String,
    pub probability_pathogenic: Probability,
    pub training_support: usize,
    pub extrapolated: bool,
    pub model: String,
}
pub trait Predictor {
    fn predict(&self, feature: &FeatureRecord) -> Result<Prediction>;
}
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct Counts {
    pub pathogenic: u32,
    pub total: u32,
}
impl Counts {
    fn add(&mut self, positive: bool) -> Result<()> {
        self.total = self
            .total
            .checked_add(1)
            .ok_or_else(|| anyhow::anyhow!("too many training observations"))?;
        if positive {
            self.pathogenic += 1;
        }
        Ok(())
    }
    fn probability(self) -> Result<Probability> {
        Probability::new(
            self.pathogenic
                .checked_add(1)
                .ok_or_else(|| anyhow::anyhow!("count overflow"))?,
            self.total
                .checked_add(2)
                .ok_or_else(|| anyhow::anyhow!("count overflow"))?,
        )
    }
}
#[derive(Debug, Serialize)]
pub struct ConsequenceModel {
    pub groups: BTreeMap<Consequence, Counts>,
    pub overall: Counts,
}
impl ConsequenceModel {
    pub fn fit(training: &[TrainingExample]) -> Result<Self> {
        ensure!(!training.is_empty(), "cannot train on an empty dataset");
        let mut groups = BTreeMap::<Consequence, Counts>::new();
        let mut overall = Counts::default();
        for example in training {
            groups
                .entry(example.feature.consequence)
                .or_default()
                .add(example.observation.label.positive())?;
            overall.add(example.observation.label.positive())?;
        }
        Ok(Self { groups, overall })
    }
}
impl Predictor for ConsequenceModel {
    fn predict(&self, feature: &FeatureRecord) -> Result<Prediction> {
        let counts = self.groups.get(&feature.consequence);
        let estimate = counts.copied().unwrap_or(self.overall);
        Ok(Prediction {
            id: feature.id.clone(),
            variant_key: feature.variant.key(),
            probability_pathogenic: estimate.probability()?,
            training_support: estimate.total as usize,
            extrapolated: counts.is_none(),
            model: "consequence-beta-binomial-v1".into(),
        })
    }
}
pub struct PrevalenceBaseline {
    counts: Counts,
}
impl PrevalenceBaseline {
    pub fn fit(training: &[TrainingExample]) -> Result<Self> {
        Ok(Self {
            counts: ConsequenceModel::fit(training)?.overall,
        })
    }
}
impl Predictor for PrevalenceBaseline {
    fn predict(&self, feature: &FeatureRecord) -> Result<Prediction> {
        Ok(Prediction {
            id: feature.id.clone(),
            variant_key: feature.variant.key(),
            probability_pathogenic: self.counts.probability()?,
            training_support: self.counts.total as usize,
            extrapolated: false,
            model: "training-prevalence-beta-binomial-v1".into(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn smoothing_is_not_certainty() {
        assert_eq!(
            Counts {
                pathogenic: 0,
                total: 3
            }
            .probability()
            .unwrap()
            .value(),
            0.2
        );
        assert_eq!(
            Counts {
                pathogenic: 3,
                total: 3
            }
            .probability()
            .unwrap()
            .value(),
            0.8
        );
    }
}
