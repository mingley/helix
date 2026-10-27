//! Group-held-out predictions and errors; every model uses identical target rows.
use crate::{
    data::{GroupSplit, MatchedPanel, OutcomeMatrix},
    model::{finite_mean, ContextFit, ModelKind, ModelPredictions},
};
use anyhow::{ensure, Result};
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Clone, Debug, Serialize)]
pub struct ContextPrediction {
    pub fold_index_0based: usize,
    pub sequence_id: String,
    pub group: String,
    pub context: String,
    pub observed: f64,
    pub predictions: ModelPredictions,
}

/// All arrays follow ContextEvaluation.model_order.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct ModelLosses {
    pub mse: [f64; 5],
    pub mae: [f64; 5],
}

#[derive(Clone, Debug, Serialize)]
pub struct ContextLosses {
    pub context: String,
    pub model_losses: ModelLosses,
}

#[derive(Clone, Debug, Serialize)]
pub struct SequenceErrors {
    pub sequence_id: String,
    pub group: String,
    /// Equal weight to the paired contexts within this sequence.
    pub model_losses: ModelLosses,
    pub by_context: Vec<ContextLosses>,
}

#[derive(Clone, Debug, Serialize)]
pub struct GroupErrors {
    pub group: String,
    pub sequences: usize,
    pub measurements: usize,
    /// Equal sequence weight within this group; equal context weight within a sequence.
    pub model_losses: ModelLosses,
    pub by_context: Vec<ContextLosses>,
}

#[derive(Clone, Debug, Serialize)]
pub struct ModelSummary {
    pub model: ModelKind,
    pub macro_group_mse: f64,
    pub macro_group_mae: f64,
    pub sequence_weighted_mse: f64,
    pub sequence_weighted_mae: f64,
    pub macro_group_by_context: Vec<ContextSummary>,
}

#[derive(Clone, Debug, Serialize)]
pub struct ContextSummary {
    pub context: String,
    pub mse: f64,
    pub mae: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct ContrastPrediction {
    pub fold_index_0based: usize,
    pub sequence_id: String,
    pub group: String,
    pub context_left: String,
    pub context_right: String,
    /// Right context minus left context, on the original endpoint scale.
    pub observed_difference: f64,
    pub predicted_difference: ModelPredictions,
}

#[derive(Clone, Debug, Serialize)]
pub struct ContrastGroupErrors {
    pub group: String,
    pub context_left: String,
    pub context_right: String,
    pub sequences: usize,
    pub model_losses: ModelLosses,
}

#[derive(Clone, Debug, Serialize)]
pub struct FoldFit {
    pub fold_index_0based: usize,
    pub split: GroupSplit,
    pub fit: ContextFit,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct BalancedDesignChecks {
    /// A shares the fitted S coefficient object, rather than fitting it again.
    pub max_abs_additive_offset_prediction_residual: f64,
    pub max_abs_mean_independent_beta_minus_shared_beta: f64,
    pub max_abs_mean_independent_prediction_minus_shared_prediction: f64,
    /// MSE(A)-MSE(I) = (C-1)/(2C) * mean_pairwise[contrastMSE(A)-contrastMSE(I)].
    pub max_abs_per_sequence_mse_contrast_identity_residual: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct ContextEvaluation {
    pub schema_version: u32,
    pub panel_metadata_sha256: String,
    pub outcomes_sha256: String,
    pub alpha: f64,
    pub model_order: [ModelKind; 5],
    pub contexts: Vec<String>,
    pub sequence_count: usize,
    pub group_count: usize,
    pub measurement_count: usize,
    pub folds: Vec<FoldFit>,
    pub predictions: Vec<ContextPrediction>,
    pub per_sequence_errors: Vec<SequenceErrors>,
    pub per_group_errors: Vec<GroupErrors>,
    pub models: Vec<ModelSummary>,
    pub context_contrast_predictions: Vec<ContrastPrediction>,
    pub context_contrast_group_errors: Vec<ContrastGroupErrors>,
    pub balanced_design_checks: BalancedDesignChecks,
    pub interpretation_limits: Vec<String>,
}

/// No fold or scientific default is inferred here: caller supplies group folds,
/// endpoint values, and the preselected alpha. Evaluation rejects incomplete OOF
/// coverage. Every feature statistic is fitted exclusively on that fold's train set.
pub fn evaluate(
    panel: &MatchedPanel,
    folds: &[GroupSplit],
    outcomes: &OutcomeMatrix,
    alpha: f64,
) -> Result<ContextEvaluation> {
    ensure!(folds.len() >= 2, "evaluation needs at least two folds");
    outcomes.verify_panel(panel)?;
    let mut counts = vec![0usize; panel.sequences().len()];
    for split in folds {
        split.verify_panel(panel)?;
        for index in split.test_indices() {
            counts[*index] += 1;
        }
    }
    ensure!(
        counts.iter().all(|count| *count == 1),
        "OOF test folds must cover every sequence exactly once"
    );
    let mut fits = Vec::with_capacity(folds.len());
    let mut predictions = Vec::with_capacity(panel.sequences().len() * panel.contexts().len());
    let mut per_sequence_errors = Vec::with_capacity(panel.sequences().len());
    let mut group_errors: BTreeMap<String, GroupAccumulator> = BTreeMap::new();
    let mut contrast_errors: BTreeMap<(String, usize, usize), Accumulator> = BTreeMap::new();
    let mut context_contrast_predictions = Vec::new();
    let mut checks = BalancedDesignChecks::default();
    for (fold_index, split) in folds.iter().enumerate() {
        let fit = ContextFit::fit(panel, split, outcomes, alpha)?;
        for j in 0..fit.shared_dna_fit.feature_count {
            let independent = fit
                .independent_context_fits
                .values()
                .map(|context| context.standardized_coefficients[j])
                .collect::<Vec<_>>();
            checks.max_abs_mean_independent_beta_minus_shared_beta =
                checks.max_abs_mean_independent_beta_minus_shared_beta.max(
                    (finite_mean(&independent)? - fit.shared_dna_fit.standardized_coefficients[j])
                        .abs(),
                );
        }
        for index in split.test_indices() {
            let sequence = &panel.sequences()[*index];
            let mut sequence_errors = Accumulator::default();
            let mut by_context = Vec::with_capacity(panel.contexts().len());
            let mut sequence_predictions = Vec::with_capacity(panel.contexts().len());
            for (j, context) in panel.contexts().iter().enumerate() {
                let prediction = fit.predict(&panel.features()[*index], context)?;
                let observed = outcomes.values()[*index][j];
                let mut context_errors = Accumulator::default();
                context_errors.push(observed, prediction)?;
                sequence_errors.push(observed, prediction)?;
                by_context.push(ContextLosses {
                    context: context.clone(),
                    model_losses: context_errors.finish()?,
                });
                let residual = prediction.additive_context_dna
                    - prediction.shared_dna
                    - (fit.training_context_means[context] - fit.global_mean);
                checks.max_abs_additive_offset_prediction_residual = checks
                    .max_abs_additive_offset_prediction_residual
                    .max(residual.abs());
                sequence_predictions.push(prediction);
                predictions.push(ContextPrediction {
                    fold_index_0based: fold_index,
                    sequence_id: sequence.id.clone(),
                    group: sequence.group.clone(),
                    context: context.clone(),
                    observed,
                    predictions: prediction,
                });
            }
            let mean_independent = finite_mean(
                &sequence_predictions
                    .iter()
                    .map(|prediction| prediction.independent_context_dna)
                    .collect::<Vec<_>>(),
            )?;
            checks.max_abs_mean_independent_prediction_minus_shared_prediction = checks
                .max_abs_mean_independent_prediction_minus_shared_prediction
                .max((mean_independent - sequence_predictions[0].shared_dna).abs());
            let losses = sequence_errors.finish()?;
            let group = group_errors
                .entry(sequence.group.clone())
                .or_insert_with(|| GroupAccumulator::new(panel.contexts().len()));
            group.sequences += 1;
            group.overall.merge(&sequence_errors);
            for (j, context) in by_context.iter().enumerate() {
                group.contexts[j].add_losses(context.model_losses);
            }
            let mut contrast_loss_difference = Vec::new();
            for left in 0..panel.contexts().len() {
                for right in left + 1..panel.contexts().len() {
                    let observed =
                        outcomes.values()[*index][right] - outcomes.values()[*index][left];
                    let prediction = subtract_predictions(
                        sequence_predictions[right],
                        sequence_predictions[left],
                    );
                    contrast_errors
                        .entry((sequence.group.clone(), left, right))
                        .or_default()
                        .push(observed, prediction)?;
                    let mut errors = Accumulator::default();
                    errors.push(observed, prediction)?;
                    let contrast = errors.finish()?;
                    contrast_loss_difference.push(contrast.mse[3] - contrast.mse[4]);
                    context_contrast_predictions.push(ContrastPrediction {
                        fold_index_0based: fold_index,
                        sequence_id: sequence.id.clone(),
                        group: sequence.group.clone(),
                        context_left: panel.contexts()[left].clone(),
                        context_right: panel.contexts()[right].clone(),
                        observed_difference: observed,
                        predicted_difference: prediction,
                    });
                }
            }
            if !contrast_loss_difference.is_empty() {
                let factor =
                    (panel.contexts().len() - 1) as f64 / (2 * panel.contexts().len()) as f64;
                checks.max_abs_per_sequence_mse_contrast_identity_residual = checks
                    .max_abs_per_sequence_mse_contrast_identity_residual
                    .max(
                        (losses.mse[3]
                            - losses.mse[4]
                            - factor * finite_mean(&contrast_loss_difference)?)
                        .abs(),
                    );
            }
            per_sequence_errors.push(SequenceErrors {
                sequence_id: sequence.id.clone(),
                group: sequence.group.clone(),
                model_losses: losses,
                by_context,
            });
        }
        fits.push(FoldFit {
            fold_index_0based: fold_index,
            split: split.clone(),
            fit,
        });
    }
    let per_group_errors: Vec<_> = group_errors
        .into_iter()
        .map(|(group, errors)| {
            Ok(GroupErrors {
                group,
                sequences: errors.sequences,
                measurements: errors.overall.count,
                model_losses: errors.overall.finish()?,
                by_context: errors
                    .contexts
                    .iter()
                    .enumerate()
                    .map(|(j, accumulator)| {
                        Ok(ContextLosses {
                            context: panel.contexts()[j].clone(),
                            model_losses: accumulator.finish()?,
                        })
                    })
                    .collect::<Result<_>>()?,
            })
        })
        .collect::<Result<_>>()?;
    let models = ModelKind::ALL
        .iter()
        .enumerate()
        .map(|(j, model)| {
            Ok(ModelSummary {
                model: *model,
                macro_group_mse: finite_mean(
                    &per_group_errors
                        .iter()
                        .map(|group| group.model_losses.mse[j])
                        .collect::<Vec<_>>(),
                )?,
                macro_group_mae: finite_mean(
                    &per_group_errors
                        .iter()
                        .map(|group| group.model_losses.mae[j])
                        .collect::<Vec<_>>(),
                )?,
                sequence_weighted_mse: finite_mean(
                    &per_sequence_errors
                        .iter()
                        .map(|sequence| sequence.model_losses.mse[j])
                        .collect::<Vec<_>>(),
                )?,
                sequence_weighted_mae: finite_mean(
                    &per_sequence_errors
                        .iter()
                        .map(|sequence| sequence.model_losses.mae[j])
                        .collect::<Vec<_>>(),
                )?,
                macro_group_by_context: panel
                    .contexts()
                    .iter()
                    .enumerate()
                    .map(|(context_index, context)| {
                        Ok(ContextSummary {
                            context: context.clone(),
                            mse: finite_mean(
                                &per_group_errors
                                    .iter()
                                    .map(|group| {
                                        group.by_context[context_index].model_losses.mse[j]
                                    })
                                    .collect::<Vec<_>>(),
                            )?,
                            mae: finite_mean(
                                &per_group_errors
                                    .iter()
                                    .map(|group| {
                                        group.by_context[context_index].model_losses.mae[j]
                                    })
                                    .collect::<Vec<_>>(),
                            )?,
                        })
                    })
                    .collect::<Result<_>>()?,
            })
        })
        .collect::<Result<_>>()?;
    let context_contrast_group_errors = contrast_errors
        .into_iter()
        .map(|((group, left, right), errors)| {
            Ok(ContrastGroupErrors {
                group,
                context_left: panel.contexts()[left].clone(),
                context_right: panel.contexts()[right].clone(),
                sequences: errors.count,
                model_losses: errors.finish()?,
            })
        })
        .collect::<Result<_>>()?;
    // Stable output order independent of the supplied fold list order within each fold.
    predictions.sort_by(|a, b| (&a.sequence_id, &a.context).cmp(&(&b.sequence_id, &b.context)));
    per_sequence_errors.sort_by(|a, b| a.sequence_id.cmp(&b.sequence_id));
    context_contrast_predictions.sort_by(|a, b| {
        (&a.sequence_id, &a.context_left, &a.context_right).cmp(&(
            &b.sequence_id,
            &b.context_left,
            &b.context_right,
        ))
    });
    Ok(ContextEvaluation{schema_version:1,panel_metadata_sha256:panel.metadata_sha256().into(),outcomes_sha256:outcomes.outcomes_sha256().into(),
        alpha,model_order:ModelKind::ALL,contexts:panel.contexts().to_vec(),sequence_count:panel.sequences().len(),group_count:per_group_errors.len(),
        measurement_count:predictions.len(),folds:fits,predictions,per_sequence_errors,per_group_errors,models,
        context_contrast_predictions,context_contrast_group_errors,balanced_design_checks:checks,
        interpretation_limits:vec![
            "Measured endpoint and source/cell/assay identities must be supplied and pinned by the caller; this generic model does not supply provenance or a study protocol.".into(),
            "Held-out groups are caller-defined. Exact/RC-identical DNA cannot cross groups; near-duplicate, source-parent and physical-neighborhood policy require an independently frozen grouping procedure.".into(),
            "Primary group weighting gives each held-out group equal weight regardless of sequence count; sequence-weighted losses are separately reported.".into(),
            "Context-specific slopes encode predictive differences on this paired reporter panel; a gain does not identify causal factors, endogenous activity or clinical effects.".into(),
            "A source may already have been explored; this API makes no claim of untouched outcomes, independent laboratory replication or discovery.".into(),
        ]})
}

fn subtract_predictions(right: ModelPredictions, left: ModelPredictions) -> ModelPredictions {
    ModelPredictions {
        global_mean: right.global_mean - left.global_mean,
        context_mean: right.context_mean - left.context_mean,
        shared_dna: right.shared_dna - left.shared_dna,
        additive_context_dna: right.additive_context_dna - left.additive_context_dna,
        independent_context_dna: right.independent_context_dna - left.independent_context_dna,
    }
}

#[derive(Default)]
struct Accumulator {
    count: usize,
    squared: [f64; 5],
    absolute: [f64; 5],
}
impl Accumulator {
    fn push(&mut self, observed: f64, predicted: ModelPredictions) -> Result<()> {
        for (j, prediction) in predicted.values().iter().enumerate() {
            let residual = observed - prediction;
            ensure!(
                residual.is_finite() && residual.powi(2).is_finite(),
                "context error overflow"
            );
            self.squared[j] += residual.powi(2);
            self.absolute[j] += residual.abs();
        }
        self.count += 1;
        Ok(())
    }
    fn merge(&mut self, other: &Self) {
        self.count += other.count;
        for j in 0..5 {
            self.squared[j] += other.squared[j];
            self.absolute[j] += other.absolute[j];
        }
    }
    fn add_losses(&mut self, losses: ModelLosses) {
        self.count += 1;
        for j in 0..5 {
            self.squared[j] += losses.mse[j];
            self.absolute[j] += losses.mae[j];
        }
    }
    fn finish(&self) -> Result<ModelLosses> {
        ensure!(self.count > 0, "empty context error accumulator");
        ensure!(
            self.squared
                .iter()
                .chain(&self.absolute)
                .all(|x| x.is_finite()),
            "context error sum overflow"
        );
        Ok(ModelLosses {
            mse: self.squared.map(|value| value / self.count as f64),
            mae: self.absolute.map(|value| value / self.count as f64),
        })
    }
}
struct GroupAccumulator {
    sequences: usize,
    overall: Accumulator,
    contexts: Vec<Accumulator>,
}
impl GroupAccumulator {
    fn new(contexts: usize) -> Self {
        Self {
            sequences: 0,
            overall: Accumulator::default(),
            contexts: (0..contexts).map(|_| Accumulator::default()).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{
        tests::{outcomes, panel},
        MeasuredOutcome, SequenceRecord,
    };

    #[test]
    fn all_models_target_the_same_held_out_rows_and_preserve_algebra() {
        let panel = panel();
        let mut rows = outcomes(&panel);
        for row in &mut rows {
            if row.context == "ctx2" {
                row.value = 13.0 - 2.0 * row.value;
            }
        }
        let matrix = OutcomeMatrix::from_records(&panel, &rows).unwrap();
        let folds = GroupSplit::partition(
            &panel,
            &[vec!["g1".into()], vec!["g2".into()], vec!["g3".into()]],
        )
        .unwrap();
        let report = evaluate(&panel, &folds, &matrix, 10.0).unwrap();
        assert_eq!(report.predictions.len(), 6);
        assert_eq!(report.per_sequence_errors.len(), 3);
        assert_eq!(report.context_contrast_predictions.len(), 3);
        assert_eq!(report.models.len(), 5);
        assert!(
            report
                .balanced_design_checks
                .max_abs_mean_independent_beta_minus_shared_beta
                < 1e-12
        );
        assert!(
            report
                .balanced_design_checks
                .max_abs_mean_independent_prediction_minus_shared_prediction
                < 1e-12
        );
        assert!(
            report
                .balanced_design_checks
                .max_abs_per_sequence_mse_contrast_identity_residual
                < 1e-12
        );
        assert!(evaluate(&panel, &folds[..2], &matrix, 10.0).is_err());
    }

    #[test]
    fn equal_group_weight_differs_from_equal_sequence_weight() {
        let panel = MatchedPanel::new(
            vec![
                SequenceRecord {
                    id: "a".into(),
                    group: "large".into(),
                    sequence: "AAAAAA".into(),
                },
                SequenceRecord {
                    id: "b".into(),
                    group: "large".into(),
                    sequence: "CCCCCC".into(),
                },
                SequenceRecord {
                    id: "c".into(),
                    group: "small".into(),
                    sequence: "ACACAC".into(),
                },
            ],
            vec!["ctx".into()],
        )
        .unwrap();
        let rows = vec![
            MeasuredOutcome {
                sequence_id: "a".into(),
                context: "ctx".into(),
                value: 0.0,
            },
            MeasuredOutcome {
                sequence_id: "b".into(),
                context: "ctx".into(),
                value: 4.0,
            },
            MeasuredOutcome {
                sequence_id: "c".into(),
                context: "ctx".into(),
                value: 10.0,
            },
        ];
        let matrix = OutcomeMatrix::from_records(&panel, &rows).unwrap();
        let folds =
            GroupSplit::partition(&panel, &[vec!["large".into()], vec!["small".into()]]).unwrap();
        let report = evaluate(&panel, &folds, &matrix, 10.0).unwrap();
        // Large held-out group errors:100,36 =>68; small group:(10-2)^2=64.
        let global = &report.models[0];
        assert_eq!(global.macro_group_mse, 66.0);
        assert!((global.sequence_weighted_mse - 200.0 / 3.0).abs() < 1e-12);
        assert!(report.context_contrast_predictions.is_empty());
    }

    #[test]
    fn three_context_interactions_obey_pairwise_contrast_mse_identity() {
        let panel = MatchedPanel::new(
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
                SequenceRecord {
                    id: "d".into(),
                    group: "g4".into(),
                    sequence: "AGAGAG".into(),
                },
            ],
            vec!["x".into(), "y".into(), "z".into()],
        )
        .unwrap();
        let rows = panel
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
                        value: (i * i + 3 * j) as f64 - (i * j * j) as f64,
                    })
            })
            .collect::<Vec<_>>();
        let outcomes = OutcomeMatrix::from_records(&panel, &rows).unwrap();
        let folds = GroupSplit::partition(
            &panel,
            &[
                vec!["g1".into()],
                vec!["g2".into()],
                vec!["g3".into()],
                vec!["g4".into()],
            ],
        )
        .unwrap();
        let report = evaluate(&panel, &folds, &outcomes, 10.0).unwrap();
        assert_eq!(report.predictions.len(), 12);
        assert_eq!(report.context_contrast_predictions.len(), 12);
        assert!(
            report
                .balanced_design_checks
                .max_abs_mean_independent_beta_minus_shared_beta
                < 1e-12
        );
        assert!(
            report
                .balanced_design_checks
                .max_abs_mean_independent_prediction_minus_shared_prediction
                < 1e-12
        );
        assert!(
            report
                .balanced_design_checks
                .max_abs_per_sequence_mse_contrast_identity_residual
                < 1e-12
        );
        assert!(report
            .predictions
            .iter()
            .any(|row| (row.predictions.additive_context_dna
                - row.predictions.independent_context_dna)
                .abs()
                > 1e-6));
    }
}
