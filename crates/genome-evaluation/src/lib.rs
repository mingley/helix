//! Evaluation receives labels only after predictions have been frozen.
use anyhow::{ensure, Context, Result};
use genome_data::LabelRecord;
use genome_inference::Prediction;
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Default, Debug, Serialize)]
pub struct Confusion {
    pub true_positive: usize,
    pub true_negative: usize,
    pub false_positive: usize,
    pub false_negative: usize,
}
#[derive(Debug, Serialize)]
pub struct CalibrationBin {
    pub lower: f64,
    pub upper: f64,
    pub count: usize,
    pub mean_prediction: Option<f64>,
    pub observed_frequency: Option<f64>,
}
#[derive(Debug, Serialize)]
pub struct Disagreement {
    pub id: String,
    pub probability_pathogenic: f64,
    pub observed_label: String,
}
#[derive(Debug, Serialize)]
pub struct Evaluation {
    pub count: usize,
    pub positive_count: usize,
    pub accuracy: f64,
    pub balanced_accuracy: Option<f64>,
    pub brier_score: f64,
    pub log_loss: f64,
    pub expected_calibration_error: f64,
    pub confusion: Confusion,
    pub calibration: Vec<CalibrationBin>,
    pub disagreements: Vec<Disagreement>,
}
pub fn evaluate(predictions: &[Prediction], labels: &[LabelRecord]) -> Result<Evaluation> {
    ensure!(
        !labels.is_empty(),
        "zero examples cannot validate a workflow"
    );
    ensure!(
        predictions.len() == labels.len(),
        "prediction/target count mismatch"
    );
    let by_id: BTreeMap<_, _> = predictions.iter().map(|p| (p.id.as_str(), p)).collect();
    ensure!(by_id.len() == predictions.len(), "duplicate prediction ID");
    let mut seen = std::collections::BTreeSet::new();
    let mut confusion = Confusion::default();
    let (mut brier, mut log_loss) = (0.0, 0.0);
    let mut bins = [(0usize, 0.0, 0.0); 5];
    let mut disagreements = Vec::new();
    for label in labels {
        ensure!(seen.insert(&label.id), "duplicate evaluation label");
        let prediction = by_id
            .get(label.id.as_str())
            .context("missing prediction for target")?;
        let p = prediction.probability_pathogenic.value();
        let actual = label.label.positive();
        let positive = prediction.probability_pathogenic.label().positive();
        match (positive, actual) {
            (true, true) => confusion.true_positive += 1,
            (false, false) => confusion.true_negative += 1,
            (true, false) => confusion.false_positive += 1,
            (false, true) => confusion.false_negative += 1,
        }
        let y = f64::from(actual);
        brier += (p - y).powi(2);
        // Explicit numerical clipping for log-loss, never a change to stored predictions.
        let q = p.clamp(1e-15, 1.0 - 1e-15);
        log_loss -= y * q.ln() + (1.0 - y) * (1.0 - q).ln();
        let index = ((p * 5.0).floor() as usize).min(4);
        bins[index].0 += 1;
        bins[index].1 += p;
        bins[index].2 += y;
        if positive != actual {
            disagreements.push(Disagreement {
                id: label.id.clone(),
                probability_pathogenic: p,
                observed_label: format!("{:?}", label.label),
            });
        }
    }
    let n = labels.len() as f64;
    let positive_count = confusion.true_positive + confusion.false_negative;
    let negative_count = confusion.true_negative + confusion.false_positive;
    let balanced_accuracy = if positive_count > 0 && negative_count > 0 {
        Some(
            0.5 * (confusion.true_positive as f64 / positive_count as f64
                + confusion.true_negative as f64 / negative_count as f64),
        )
    } else {
        None
    };
    let mut ece = 0.0;
    let calibration = bins
        .into_iter()
        .enumerate()
        .map(|(i, (count, sum, actual))| {
            if count > 0 {
                ece += (sum - actual).abs() / n;
            }
            CalibrationBin {
                lower: i as f64 / 5.0,
                upper: (i + 1) as f64 / 5.0,
                count,
                mean_prediction: (count > 0).then(|| sum / count as f64),
                observed_frequency: (count > 0).then(|| actual / count as f64),
            }
        })
        .collect();
    Ok(Evaluation {
        count: labels.len(),
        positive_count,
        accuracy: (confusion.true_positive + confusion.true_negative) as f64 / n,
        balanced_accuracy,
        brier_score: brier / n,
        log_loss: log_loss / n,
        expected_calibration_error: ece,
        confusion,
        calibration,
        disagreements,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use genome_types::{ClinicalLabel, Probability};
    #[test]
    fn metrics_match_hand_calculation_and_reject_missing_targets() {
        let labels: Vec<_> = [ClinicalLabel::Pathogenic, ClinicalLabel::Benign]
            .into_iter()
            .enumerate()
            .map(|(i, label)| LabelRecord {
                id: i.to_string(),
                label,
                classification: String::new(),
                review_status: String::new(),
                source_url: String::new(),
            })
            .collect();
        let predictions: Vec<_> = [
            Probability::new(4, 5).unwrap(),
            Probability::new(2, 5).unwrap(),
        ]
        .into_iter()
        .enumerate()
        .map(|(i, p)| Prediction {
            id: i.to_string(),
            variant_key: i.to_string(),
            probability_pathogenic: p,
            training_support: 1,
            extrapolated: false,
            model: "test".into(),
        })
        .collect();
        let result = evaluate(&predictions, &labels).unwrap();
        assert_eq!(result.accuracy, 1.0);
        assert!((result.brier_score - 0.1).abs() < 1e-12);
        assert!((result.expected_calibration_error - 0.3).abs() < 1e-12);
        assert_eq!(result.confusion.true_positive, 1);
        assert!(evaluate(&predictions[..1], &labels).is_err());
        assert!(evaluate(&[], &[]).is_err());
    }
}
