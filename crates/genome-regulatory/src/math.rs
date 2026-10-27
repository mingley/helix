//! Small deterministic ridge solver and descriptive regression metrics.
use anyhow::{ensure, Result};
use serde::Serialize;

const MAX_FEATURES: usize = 256;

#[derive(Clone, Debug, Serialize)]
pub struct StandardizedRidge {
    pub training_rows: usize,
    pub feature_count: usize,
    pub alpha: f64,
    pub intercept: f64,
    pub training_feature_means: Vec<f64>,
    pub training_feature_scales: Vec<f64>,
    pub standardized_coefficients: Vec<f64>,
}

impl StandardizedRidge {
    /// Minimize sum squared error + alpha * ||beta||², after training-only
    /// population-SD standardization. The response mean is an unpenalized
    /// intercept. A constant feature receives scale one and coefficient zero.
    pub fn fit(features: &[Vec<f64>], outcomes: &[f64], alpha: f64) -> Result<Self> {
        RidgeDesign::prepare(features, alpha)?.fit(outcomes)
    }

    pub fn predict(&self, row: &[f64]) -> Result<f64> {
        ensure!(
            row.len() == self.feature_count && row.iter().all(|x| x.is_finite()),
            "invalid ridge prediction row"
        );
        let value = self.intercept
            + row
                .iter()
                .enumerate()
                .map(|(j, x)| {
                    (x - self.training_feature_means[j]) / self.training_feature_scales[j]
                        * self.standardized_coefficients[j]
                })
                .sum::<f64>();
        ensure!(value.is_finite(), "ridge prediction overflow");
        Ok(value)
    }
}

/// Outcome-independent design preparation reused for real and null fits.
pub struct RidgeDesign {
    training_rows: usize,
    feature_count: usize,
    alpha: f64,
    means: Vec<f64>,
    scales: Vec<f64>,
    standardized_rows: Vec<Vec<f64>>,
    lower: Vec<Vec<f64>>,
}
impl RidgeDesign {
    pub fn prepare(features: &[Vec<f64>], alpha: f64) -> Result<Self> {
        ensure!(
            !features.is_empty(),
            "ridge requires nonempty training rows"
        );
        ensure!(
            alpha.is_finite() && alpha > 0.0,
            "ridge alpha must be finite and positive"
        );
        let width = features[0].len();
        ensure!(width <= MAX_FEATURES, "ridge feature budget exceeded");
        ensure!(
            features
                .iter()
                .all(|row| row.len() == width && row.iter().all(|x| x.is_finite())),
            "nonfinite or inconsistent ridge features"
        );
        let n = features.len() as f64;
        let mut means = vec![0.0; width];
        for row in features {
            for (j, x) in row.iter().enumerate() {
                means[j] += x / n;
            }
        }
        let mut scales = vec![0.0; width];
        for row in features {
            for (j, x) in row.iter().enumerate() {
                scales[j] += (x - means[j]).powi(2) / n;
            }
        }
        for scale in &mut scales {
            *scale = scale.sqrt();
            if *scale <= 1e-12 {
                *scale = 1.0;
            }
        }
        ensure!(
            means.iter().chain(&scales).all(|x| x.is_finite()),
            "ridge centering/scale overflow"
        );
        let standardized_rows: Vec<Vec<f64>> = features
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .map(|(j, x)| (x - means[j]) / scales[j])
                    .collect()
            })
            .collect();
        let mut gram = vec![vec![0.0; width]; width];
        for row in &standardized_rows {
            for j in 0..width {
                for k in 0..=j {
                    gram[j][k] += row[j] * row[k];
                }
            }
        }
        for j in 0..width {
            gram[j][j] += alpha;
            for k in 0..j {
                gram[k][j] = gram[j][k];
            }
        }
        let mut lower = vec![vec![0.0; width]; width];
        for i in 0..width {
            for j in 0..=i {
                let product: f64 = (0..j).map(|k| lower[i][k] * lower[j][k]).sum();
                if i == j {
                    let diagonal = gram[i][i] - product;
                    ensure!(
                        diagonal.is_finite() && diagonal > 0.0,
                        "ridge matrix is not numerically positive definite"
                    );
                    lower[i][j] = diagonal.sqrt();
                } else {
                    lower[i][j] = (gram[i][j] - product) / lower[j][j];
                }
            }
        }
        Ok(Self {
            training_rows: features.len(),
            feature_count: width,
            alpha,
            means,
            scales,
            standardized_rows,
            lower,
        })
    }
    pub fn fit(&self, outcomes: &[f64]) -> Result<StandardizedRidge> {
        ensure!(
            outcomes.len() == self.training_rows && outcomes.iter().all(|y| y.is_finite()),
            "nonfinite or inconsistent ridge outcomes"
        );
        let intercept = outcomes
            .iter()
            .map(|y| y / self.training_rows as f64)
            .sum::<f64>();
        ensure!(intercept.is_finite(), "ridge outcome centering overflow");
        let mut rhs = vec![0.0; self.feature_count];
        for (row, y) in self.standardized_rows.iter().zip(outcomes) {
            for (j, x) in row.iter().enumerate() {
                rhs[j] += x * (y - intercept);
            }
        }
        let mut intermediate = vec![0.0; self.feature_count];
        for i in 0..self.feature_count {
            let product: f64 = (0..i).map(|j| self.lower[i][j] * intermediate[j]).sum();
            intermediate[i] = (rhs[i] - product) / self.lower[i][i];
        }
        let mut beta = vec![0.0; self.feature_count];
        for i in (0..self.feature_count).rev() {
            let product: f64 = ((i + 1)..self.feature_count)
                .map(|j| self.lower[j][i] * beta[j])
                .sum();
            beta[i] = (intermediate[i] - product) / self.lower[i][i];
        }
        ensure!(beta.iter().all(|x| x.is_finite()), "ridge solve overflow");
        Ok(StandardizedRidge {
            training_rows: self.training_rows,
            feature_count: self.feature_count,
            alpha: self.alpha,
            intercept,
            training_feature_means: self.means.clone(),
            training_feature_scales: self.scales.clone(),
            standardized_coefficients: beta,
        })
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct RegressionMetrics {
    pub count: usize,
    pub mse: f64,
    pub mae: f64,
    /// Undefined for constant observed values; negative values are retained.
    pub r_squared: Option<f64>,
    /// Undefined if either argument is constant.
    pub pearson: Option<f64>,
    /// Average ranks for exact ties; undefined for a constant rank vector.
    pub spearman: Option<f64>,
}

pub fn regression_metrics(observed: &[f64], predicted: &[f64]) -> Result<RegressionMetrics> {
    ensure!(
        !observed.is_empty() && observed.len() == predicted.len(),
        "metrics require paired nonempty observations"
    );
    ensure!(
        observed.iter().chain(predicted).all(|x| x.is_finite()),
        "nonfinite metric input"
    );
    let n = observed.len() as f64;
    let mean = observed.iter().map(|y| y / n).sum::<f64>();
    let sse = observed
        .iter()
        .zip(predicted)
        .map(|(y, p)| (y - p).powi(2))
        .sum::<f64>();
    let mae = observed
        .iter()
        .zip(predicted)
        .map(|(y, p)| (y - p).abs() / n)
        .sum::<f64>();
    let sst = observed.iter().map(|y| (y - mean).powi(2)).sum::<f64>();
    ensure!(
        sse.is_finite() && mae.is_finite() && sst.is_finite(),
        "regression metric overflow"
    );
    let ranks_y = average_ranks(observed)?;
    let ranks_p = average_ranks(predicted)?;
    Ok(RegressionMetrics {
        count: observed.len(),
        mse: sse / n,
        mae,
        r_squared: (sst > 0.0 && !observed.iter().all(|y| *y == observed[0]))
            .then(|| 1.0 - sse / sst),
        pearson: pearson(observed, predicted),
        spearman: pearson(&ranks_y, &ranks_p),
    })
}

/// The observed rank vector can come from exact source-decimal position means.
/// Predicted values retain their actual f64 ordering; no epsilon or rounding.
pub fn regression_metrics_with_observed_ranks(
    observed: &[f64],
    predicted: &[f64],
    observed_ranks: &[f64],
) -> Result<RegressionMetrics> {
    ensure!(
        observed_ranks.len() == observed.len()
            && observed_ranks.iter().all(|rank| rank.is_finite()),
        "invalid exact observed ranks"
    );
    let mut metrics = regression_metrics(observed, predicted)?;
    metrics.spearman = pearson(observed_ranks, &average_ranks(predicted)?);
    Ok(metrics)
}

fn pearson(left: &[f64], right: &[f64]) -> Option<f64> {
    if left.iter().all(|value| *value == left[0]) || right.iter().all(|value| *value == right[0]) {
        return None;
    }
    let n = left.len() as f64;
    let a = left.iter().map(|x| x / n).sum::<f64>();
    let b = right.iter().map(|x| x / n).sum::<f64>();
    let covariance = left
        .iter()
        .zip(right)
        .map(|(x, y)| (x - a) * (y - b))
        .sum::<f64>();
    let aa = left.iter().map(|x| (x - a).powi(2)).sum::<f64>();
    let bb = right.iter().map(|x| (x - b).powi(2)).sum::<f64>();
    if aa > 0.0 && bb > 0.0 {
        let value = covariance / aa.sqrt() / bb.sqrt();
        value.is_finite().then(|| value.clamp(-1.0, 1.0))
    } else {
        None
    }
}

fn average_ranks(values: &[f64]) -> Result<Vec<f64>> {
    ensure!(
        values.iter().all(|x| x.is_finite()),
        "cannot rank nonfinite observations"
    );
    let mut indexes: Vec<_> = (0..values.len()).collect();
    indexes.sort_by(|&a, &b| values[a].total_cmp(&values[b]));
    let mut ranks = vec![0.0; values.len()];
    let mut start = 0;
    while start < indexes.len() {
        let mut end = start + 1;
        while end < indexes.len() && values[indexes[start]] == values[indexes[end]] {
            end += 1;
        }
        let rank = (start + end + 1) as f64 / 2.0;
        for &index in &indexes[start..end] {
            ranks[index] = rank;
        }
        start = end;
    }
    Ok(ranks)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ridge_matches_hand_solution_with_unpenalized_intercept() {
        let fit = StandardizedRidge::fit(&[vec![-1.0], vec![1.0]], &[1.0, 3.0], 2.0).unwrap();
        assert_eq!(fit.intercept, 2.0);
        assert!((fit.standardized_coefficients[0] - 0.5).abs() < 1e-12);
        assert!((fit.predict(&[-1.0]).unwrap() - 1.5).abs() < 1e-12);
        assert!((fit.predict(&[1.0]).unwrap() - 2.5).abs() < 1e-12);
    }
    #[test]
    fn standardization_is_training_only_and_constant_features_stay_zero() {
        let fit =
            StandardizedRidge::fit(&[vec![-2.0, 3.0], vec![2.0, 3.0]], &[1.0, 3.0], 2.0).unwrap();
        assert_eq!(fit.training_feature_means, vec![0.0, 3.0]);
        assert_eq!(fit.training_feature_scales, vec![2.0, 1.0]);
        assert_eq!(fit.standardized_coefficients[1], 0.0);
        assert!((fit.predict(&[2.0, 1_000_000.0]).unwrap() - 2.5).abs() < 1e-12);
        let shifted =
            StandardizedRidge::fit(&[vec![-20.0], vec![20.0]], &[11.0, 13.0], 2.0).unwrap();
        assert!((shifted.predict(&[20.0]).unwrap() - 12.5).abs() < 1e-12);
    }
    #[test]
    fn invalid_inputs_and_negative_r_squared_are_not_hidden() {
        assert!(StandardizedRidge::fit(&[], &[], 1.0).is_err());
        assert!(StandardizedRidge::fit(&[vec![f64::NAN]], &[1.0], 1.0).is_err());
        assert!(StandardizedRidge::fit(&[vec![1.0]], &[1.0], 0.0).is_err());
        let result = regression_metrics(&[1.0, 2.0], &[10.0, 10.0]).unwrap();
        assert!(result.r_squared.unwrap() < 0.0);
        assert!(result.pearson.is_none());
        assert!(result.spearman.is_none());
    }
    #[test]
    fn metrics_and_tied_ranks_have_hand_checked_values() {
        let result = regression_metrics(&[1.0, 2.0, 3.0], &[1.0, 2.0, 4.0]).unwrap();
        assert!((result.mse - 1.0 / 3.0).abs() < 1e-12);
        assert!((result.mae - 1.0 / 3.0).abs() < 1e-12);
        assert!((result.r_squared.unwrap() - 0.5).abs() < 1e-12);
        assert_eq!(
            average_ranks(&[4.0, 1.0, 1.0, 3.0]).unwrap(),
            vec![4.0, 1.5, 1.5, 3.0]
        );
        assert!((result.spearman.unwrap() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn decimal_constant_vectors_are_undefined_despite_centering_roundoff() {
        let result = regression_metrics(&[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6], &[0.1; 7]).unwrap();
        assert!(result.pearson.is_none());
        assert!(result.spearman.is_none());
        let result = regression_metrics(&[0.1; 7], &[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]).unwrap();
        assert!(result.r_squared.is_none());
        assert!(result.pearson.is_none());
        assert!(result.spearman.is_none());
    }
}
