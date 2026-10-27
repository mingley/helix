//! Paired group resampling with caller-specified choices and retained common draws.
use crate::{
    evaluation::GroupErrors,
    model::{finite_mean, ModelKind},
};
use anyhow::{ensure, Result};
use genome_research_stats::SplitMix64;
use serde::Serialize;
use std::collections::BTreeSet;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize)]
pub struct ModelComparison {
    pub baseline: ModelKind,
    pub candidate: ModelKind,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum LossWeighting {
    /// Each sampled component contributes its sequence count to the denominator.
    EqualSequences,
    /// Each sampled component mean contributes equally, regardless of size.
    EqualGroups,
}

#[derive(Clone, Debug, Serialize)]
pub struct BootstrapComparison {
    pub comparison: ModelComparison,
    /// Positive means the candidate's held-out error is lower.
    pub observed_mse_improvement: f64,
    pub observed_mae_improvement: f64,
    pub mse_interval: [f64; 2],
    pub mae_interval: [f64; 2],
}

#[derive(Clone, Debug, Serialize)]
pub struct BootstrapModelSummary {
    pub model: ModelKind,
    pub observed_mse: f64,
    pub observed_mae: f64,
    pub mse_interval: [f64; 2],
    pub mae_interval: [f64; 2],
}

#[derive(Clone, Debug, Serialize)]
pub struct GroupBootstrap {
    pub group_order: Vec<String>,
    /// Every model statistic array follows this fixed order.
    pub model_order: [ModelKind; 5],
    pub comparison_order: Vec<ModelComparison>,
    pub replicates: usize,
    pub seed: u64,
    pub interval_probabilities: [f64; 2],
    pub weighting: LossWeighting,
    /// One index table shared across all comparisons and both error measures.
    pub resampling_indices: Vec<Vec<usize>>,
    pub sampled_mse_improvements: Vec<Vec<f64>>,
    pub sampled_mae_improvements: Vec<Vec<f64>>,
    pub observed_model_mse: [f64; 5],
    pub observed_model_mae: [f64; 5],
    pub sampled_model_mse: Vec<[f64; 5]>,
    pub sampled_model_mae: Vec<[f64; 5]>,
    pub models: Vec<BootstrapModelSummary>,
    pub comparisons: Vec<BootstrapComparison>,
    pub quantile_method: &'static str,
    pub interpretation_limit: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct BootstrapIntervals {
    pub interval_probabilities: [f64; 2],
    pub comparisons: Vec<BootstrapComparison>,
    pub models: Vec<BootstrapModelSummary>,
}

impl GroupBootstrap {
    /// Additional prespecified percentile levels using the same saved statistics
    /// and draws, without resampling or fitting again.
    pub fn intervals_at(&self, probabilities: [f64; 2]) -> Result<BootstrapIntervals> {
        validate_probabilities(probabilities)?;
        let comparisons = self
            .comparisons
            .iter()
            .enumerate()
            .map(|(j, comparison)| BootstrapComparison {
                comparison: comparison.comparison,
                observed_mse_improvement: comparison.observed_mse_improvement,
                observed_mae_improvement: comparison.observed_mae_improvement,
                mse_interval: interval(&self.sampled_mse_improvements, j, probabilities),
                mae_interval: interval(&self.sampled_mae_improvements, j, probabilities),
            })
            .collect();
        Ok(BootstrapIntervals {
            interval_probabilities: probabilities,
            comparisons,
            models: model_summaries(
                self.observed_model_mse,
                self.observed_model_mae,
                &self.sampled_model_mse,
                &self.sampled_model_mae,
                probabilities,
            ),
        })
    }
}

/// Sample G groups with replacement, using explicitly chosen loss weights.
/// This does not refit models or resample individual contexts, sequences or tags.
/// No study defaults are selected: every design choice is an explicit argument.
pub fn paired_group_bootstrap(
    groups: &[GroupErrors],
    comparisons: &[ModelComparison],
    replicates: usize,
    seed: u64,
    interval_probabilities: [f64; 2],
    weighting: LossWeighting,
) -> Result<GroupBootstrap> {
    ensure!(
        !groups.is_empty() && !comparisons.is_empty(),
        "bootstrap needs groups and comparisons"
    );
    ensure!(
        (1..=100_000).contains(&replicates),
        "invalid bootstrap replicate count"
    );
    ensure!(
        groups
            .len()
            .checked_mul(replicates)
            .is_some_and(|count| count <= 50_000_000),
        "retained bootstrap draw table exceeds the 50-million-index resource budget"
    );
    validate_probabilities(interval_probabilities)?;
    let unique: BTreeSet<_> = comparisons.iter().copied().collect();
    ensure!(
        unique.len() == comparisons.len()
            && comparisons
                .iter()
                .all(|comparison| comparison.baseline != comparison.candidate),
        "duplicate or self model comparison"
    );
    let mut groups: Vec<_> = groups.iter().collect();
    groups.sort_by(|a, b| a.group.cmp(&b.group));
    ensure!(
        groups.windows(2).all(|pair| pair[0].group != pair[1].group),
        "duplicate bootstrap group"
    );
    for group in &groups {
        ensure!(
            group.sequences > 0
                && group
                    .model_losses
                    .mse
                    .iter()
                    .chain(&group.model_losses.mae)
                    .all(|value| value.is_finite() && *value >= 0.0),
            "invalid held-out group losses"
        );
    }
    let mut mse = vec![vec![0.0; comparisons.len()]; groups.len()];
    let mut mae = mse.clone();
    let weights: Vec<f64> = groups
        .iter()
        .map(|group| match weighting {
            LossWeighting::EqualSequences => group.sequences as f64,
            LossWeighting::EqualGroups => 1.0,
        })
        .collect();
    for (i, group) in groups.iter().enumerate() {
        for (j, comparison) in comparisons.iter().enumerate() {
            let baseline = model_index(comparison.baseline);
            let candidate = model_index(comparison.candidate);
            mse[i][j] = group.model_losses.mse[baseline] - group.model_losses.mse[candidate];
            mae[i][j] = group.model_losses.mae[baseline] - group.model_losses.mae[candidate];
        }
    }
    let mut rng = SplitMix64::new(seed);
    let mut resampling_indices = Vec::with_capacity(replicates);
    let mut sampled_mse_improvements = Vec::with_capacity(replicates);
    let mut sampled_mae_improvements = Vec::with_capacity(replicates);
    let mut observed_model_mse = [0.0; 5];
    let mut observed_model_mae = [0.0; 5];
    for j in 0..5 {
        observed_model_mse[j] = weighted_mean(
            &groups
                .iter()
                .map(|group| group.model_losses.mse[j])
                .collect::<Vec<_>>(),
            &weights,
        )?;
        observed_model_mae[j] = weighted_mean(
            &groups
                .iter()
                .map(|group| group.model_losses.mae[j])
                .collect::<Vec<_>>(),
            &weights,
        )?;
    }
    let mut sampled_model_mse = Vec::with_capacity(replicates);
    let mut sampled_model_mae = Vec::with_capacity(replicates);
    for _ in 0..replicates {
        let indices = (0..groups.len())
            .map(|_| rng.index(groups.len()))
            .collect::<Result<Vec<_>>>()?;
        let mut sampled_mse = Vec::with_capacity(comparisons.len());
        let mut sampled_mae = Vec::with_capacity(comparisons.len());
        let sampled_weights: Vec<_> = indices.iter().map(|i| weights[*i]).collect();
        for j in 0..comparisons.len() {
            sampled_mse.push(weighted_mean(
                &indices.iter().map(|i| mse[*i][j]).collect::<Vec<_>>(),
                &sampled_weights,
            )?);
            sampled_mae.push(weighted_mean(
                &indices.iter().map(|i| mae[*i][j]).collect::<Vec<_>>(),
                &sampled_weights,
            )?);
        }
        // Reuse these exact indices and weights; this reporting work consumes no RNG.
        let mut model_mse = [0.0; 5];
        let mut model_mae = [0.0; 5];
        for j in 0..5 {
            model_mse[j] = weighted_mean(
                &indices
                    .iter()
                    .map(|i| groups[*i].model_losses.mse[j])
                    .collect::<Vec<_>>(),
                &sampled_weights,
            )?;
            model_mae[j] = weighted_mean(
                &indices
                    .iter()
                    .map(|i| groups[*i].model_losses.mae[j])
                    .collect::<Vec<_>>(),
                &sampled_weights,
            )?;
        }
        resampling_indices.push(indices);
        sampled_mse_improvements.push(sampled_mse);
        sampled_mae_improvements.push(sampled_mae);
        sampled_model_mse.push(model_mse);
        sampled_model_mae.push(model_mae);
    }
    let summaries = comparisons
        .iter()
        .enumerate()
        .map(|(j, comparison)| {
            Ok(BootstrapComparison {
                comparison: *comparison,
                observed_mse_improvement: weighted_mean(
                    &mse.iter().map(|row| row[j]).collect::<Vec<_>>(),
                    &weights,
                )?,
                observed_mae_improvement: weighted_mean(
                    &mae.iter().map(|row| row[j]).collect::<Vec<_>>(),
                    &weights,
                )?,
                mse_interval: interval(&sampled_mse_improvements, j, interval_probabilities),
                mae_interval: interval(&sampled_mae_improvements, j, interval_probabilities),
            })
        })
        .collect::<Result<_>>()?;
    let models = model_summaries(
        observed_model_mse,
        observed_model_mae,
        &sampled_model_mse,
        &sampled_model_mae,
        interval_probabilities,
    );
    Ok(GroupBootstrap{group_order:groups.iter().map(|group|group.group.clone()).collect(),model_order:ModelKind::ALL,comparison_order:comparisons.to_vec(),
        replicates,seed,interval_probabilities,weighting,resampling_indices,sampled_mse_improvements,sampled_mae_improvements,
        observed_model_mse,observed_model_mae,sampled_model_mse,sampled_model_mae,models,comparisons:summaries,
        quantile_method:"Type7 linear interpolation at the explicitly supplied probabilities",
        interpretation_limit:"Conditional paired resampling of fixed held-out group losses, not independent-study replication, model-refit uncertainty, guaranteed population coverage or a simultaneous interval family.",})
}

fn model_summaries(
    observed_mse: [f64; 5],
    observed_mae: [f64; 5],
    sampled_mse: &[[f64; 5]],
    sampled_mae: &[[f64; 5]],
    probabilities: [f64; 2],
) -> Vec<BootstrapModelSummary> {
    ModelKind::ALL
        .iter()
        .enumerate()
        .map(|(j, model)| BootstrapModelSummary {
            model: *model,
            observed_mse: observed_mse[j],
            observed_mae: observed_mae[j],
            mse_interval: interval(sampled_mse, j, probabilities),
            mae_interval: interval(sampled_mae, j, probabilities),
        })
        .collect()
}

fn weighted_mean(values: &[f64], weights: &[f64]) -> Result<f64> {
    if weights.iter().all(|weight| *weight == weights[0]) {
        return finite_mean(values);
    }
    let denominator: f64 = weights.iter().sum();
    let mean: f64 = values
        .iter()
        .zip(weights)
        .map(|(value, weight)| value * (weight / denominator))
        .sum();
    ensure!(mean.is_finite(), "weighted loss overflow");
    Ok(mean)
}

fn validate_probabilities(probabilities: [f64; 2]) -> Result<()> {
    ensure!(
        probabilities
            .iter()
            .all(|p| p.is_finite() && (0.0..=1.0).contains(p))
            && probabilities[0] < probabilities[1],
        "invalid interval probabilities"
    );
    Ok(())
}

fn model_index(model: ModelKind) -> usize {
    ModelKind::ALL
        .iter()
        .position(|value| *value == model)
        .expect("exhaustive model list")
}
fn interval<T: AsRef<[f64]>>(
    statistics: &[T],
    comparison: usize,
    probabilities: [f64; 2],
) -> [f64; 2] {
    let mut sorted: Vec<_> = statistics
        .iter()
        .map(|row| row.as_ref()[comparison])
        .collect();
    sorted.sort_by(f64::total_cmp);
    probabilities.map(|probability| {
        let location = (sorted.len() - 1) as f64 * probability;
        let lower = location.floor() as usize;
        let upper = location.ceil() as usize;
        sorted[lower] + (location - lower as f64) * (sorted[upper] - sorted[lower])
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evaluation::ModelLosses;

    #[test]
    fn all_model_statistics_reuse_the_same_draw_stream_and_explicit_weights() {
        let groups = vec![
            GroupErrors {
                group: "a".into(),
                sequences: 1,
                measurements: 3,
                model_losses: ModelLosses {
                    mse: [1.0, 4.0, 9.0, 16.0, 25.0],
                    mae: [0.5, 1.0, 1.5, 2.0, 2.5],
                },
                by_context: vec![],
            },
            GroupErrors {
                group: "b".into(),
                sequences: 3,
                measurements: 9,
                model_losses: ModelLosses {
                    mse: [6.0, 2.0, 7.0, 1.0, 8.0],
                    mae: [1.0, 0.5, 2.0, 0.25, 1.25],
                },
                by_context: vec![],
            },
            GroupErrors {
                group: "c".into(),
                sequences: 2,
                measurements: 6,
                model_losses: ModelLosses {
                    mse: [3.0, 8.0, 2.0, 11.0, 4.0],
                    mae: [0.75, 1.75, 0.25, 2.5, 1.0],
                },
                by_context: vec![],
            },
        ];
        let comparisons = [ModelComparison {
            baseline: ModelKind::AdditiveContextDna,
            candidate: ModelKind::IndependentContextDna,
        }];
        for weighting in [LossWeighting::EqualSequences, LossWeighting::EqualGroups] {
            let report = paired_group_bootstrap(
                &groups,
                &comparisons,
                17,
                0x2026102702,
                [0.025, 0.975],
                weighting,
            )
            .unwrap();
            let mut independent_stream = SplitMix64::new(0x2026102702);
            let weights: Vec<f64> = groups
                .iter()
                .map(|group| match weighting {
                    LossWeighting::EqualSequences => group.sequences as f64,
                    LossWeighting::EqualGroups => 1.0,
                })
                .collect();
            let expected = |indices: &[usize], j: usize, mae: bool| {
                let numerator: f64 = indices
                    .iter()
                    .map(|i| {
                        weights[*i]
                            * if mae {
                                groups[*i].model_losses.mae[j]
                            } else {
                                groups[*i].model_losses.mse[j]
                            }
                    })
                    .sum();
                numerator / indices.iter().map(|i| weights[*i]).sum::<f64>()
            };
            assert_eq!(report.model_order, ModelKind::ALL);
            assert_eq!(report.models.len(), 5);
            for j in 0..5 {
                assert!(
                    (report.observed_model_mse[j] - expected(&[0, 1, 2], j, false)).abs() < 1e-12
                );
                assert!(
                    (report.observed_model_mae[j] - expected(&[0, 1, 2], j, true)).abs() < 1e-12
                );
            }
            for (replicate, indices) in report.resampling_indices.iter().enumerate() {
                assert_eq!(
                    *indices,
                    (0..3)
                        .map(|_| independent_stream.index(3).unwrap())
                        .collect::<Vec<_>>()
                );
                for j in 0..5 {
                    assert!(
                        (report.sampled_model_mse[replicate][j] - expected(indices, j, false))
                            .abs()
                            < 1e-12
                    );
                    assert!(
                        (report.sampled_model_mae[replicate][j] - expected(indices, j, true)).abs()
                            < 1e-12
                    );
                }
                assert!(
                    (report.sampled_mse_improvements[replicate][0]
                        - (report.sampled_model_mse[replicate][3]
                            - report.sampled_model_mse[replicate][4]))
                        .abs()
                        < 1e-12
                );
                assert!(
                    (report.sampled_mae_improvements[replicate][0]
                        - (report.sampled_model_mae[replicate][3]
                            - report.sampled_model_mae[replicate][4]))
                        .abs()
                        < 1e-12
                );
            }
            let ranges = report.intervals_at([0.0, 1.0]).unwrap();
            assert_eq!(ranges.models.len(), 5);
            for (j, model) in ranges.models.iter().enumerate() {
                assert_eq!(model.observed_mse, report.observed_model_mse[j]);
                assert_eq!(model.observed_mae, report.observed_model_mae[j]);
                assert_eq!(
                    model.mse_interval,
                    [
                        report
                            .sampled_model_mse
                            .iter()
                            .map(|row| row[j])
                            .fold(f64::INFINITY, f64::min),
                        report
                            .sampled_model_mse
                            .iter()
                            .map(|row| row[j])
                            .fold(f64::NEG_INFINITY, f64::max)
                    ]
                );
            }
        }
    }

    #[test]
    fn common_draws_preserve_paired_opposite_comparisons_and_group_order() {
        let group = |id: &str, value: f64| GroupErrors {
            group: id.into(),
            sequences: 1,
            measurements: 2,
            model_losses: ModelLosses {
                mse: [value, 0.0, 0.0, 0.0, 0.0],
                mae: [value, 0.0, 0.0, 0.0, 0.0],
            },
            by_context: vec![],
        };
        let groups = vec![group("b", 3.0), group("a", 1.0)];
        let comparisons = [
            ModelComparison {
                baseline: ModelKind::GlobalMean,
                candidate: ModelKind::ContextMean,
            },
            ModelComparison {
                baseline: ModelKind::ContextMean,
                candidate: ModelKind::GlobalMean,
            },
        ];
        let report = paired_group_bootstrap(
            &groups,
            &comparisons,
            20,
            123,
            [0.025, 0.975],
            LossWeighting::EqualGroups,
        )
        .unwrap();
        assert_eq!(report.group_order, vec!["a", "b"]);
        assert_eq!(report.comparisons[0].observed_mse_improvement, 2.0);
        for ((indices, mse), mae) in report
            .resampling_indices
            .iter()
            .zip(&report.sampled_mse_improvements)
            .zip(&report.sampled_mae_improvements)
        {
            assert_eq!(mse[0], -mse[1]);
            assert_eq!(mse, mae);
            assert_eq!(
                mse[0],
                indices
                    .iter()
                    .map(|i| if *i == 0 { 0.5 } else { 1.5 })
                    .sum::<f64>()
            );
        }
        let repeat = paired_group_bootstrap(
            &groups,
            &comparisons,
            20,
            123,
            [0.025, 0.975],
            LossWeighting::EqualGroups,
        )
        .unwrap();
        assert_eq!(report.resampling_indices, repeat.resampling_indices);
        let wider = report.intervals_at([0.0, 1.0]).unwrap();
        assert!(wider.comparisons[0].mse_interval[0] <= report.comparisons[0].mse_interval[0]);
        assert!(wider.comparisons[0].mse_interval[1] >= report.comparisons[0].mse_interval[1]);
        let mut unequal = groups;
        unequal[0].sequences = 3;
        let weighted = paired_group_bootstrap(
            &unequal,
            &comparisons,
            20,
            123,
            [0.025, 0.975],
            LossWeighting::EqualSequences,
        )
        .unwrap();
        assert_eq!(weighted.resampling_indices, report.resampling_indices);
        assert_eq!(weighted.comparisons[0].observed_mse_improvement, 2.5);
        for (indices, losses) in weighted
            .resampling_indices
            .iter()
            .zip(&weighted.sampled_mse_improvements)
        {
            let numerator: f64 = indices
                .iter()
                .map(|i| if *i == 0 { 1.0 } else { 9.0 })
                .sum();
            let denominator: f64 = indices
                .iter()
                .map(|i| if *i == 0 { 1.0 } else { 3.0 })
                .sum();
            assert!((losses[0] - numerator / denominator).abs() < 1e-12);
        }
    }
}
