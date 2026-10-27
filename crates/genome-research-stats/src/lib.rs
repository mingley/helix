//! Descriptive uncertainty and paired comparisons for small research datasets.
//!
//! Calculations do not establish random sampling, independent cohorts, empirical
//! validity or discovery. Callers must define the sampling unit and comparison.
use anyhow::{ensure, Result};
use serde::Serialize;

mod regulatory;
pub use regulatory::*;

const Z_975: f64 = 1.959_963_984_540_054;
const MAX_EXACT_F64_INTEGER: u64 = 1 << 53;
const MAX_DISCORDANT_PAIRS: usize = 1_000_000;

#[derive(Clone, Debug, Serialize)]
pub struct ProportionInterval {
    pub successes: usize,
    pub total: usize,
    pub estimate: f64,
    pub lower: f64,
    pub upper: f64,
    pub confidence_level: f64,
    pub method: &'static str,
    pub assumption: &'static str,
}

/// Wilson score interval at nominal 95% confidence. Independent, identically
/// distributed Bernoulli sampling is assumed, not inferred from counts. On a
/// selected fixture this is descriptive sensitivity to sample size, not a
/// population confidence guarantee.
pub fn proportion_interval(successes: usize, total: usize) -> Result<ProportionInterval> {
    ensure!(total > 0, "a proportion requires at least one observation");
    ensure!(successes <= total, "successes exceed total observations");
    ensure!(
        total as u128 <= u128::from(MAX_EXACT_F64_INTEGER),
        "sample count exceeds exact floating-point integer representation"
    );
    let n = total as f64;
    let estimate = successes as f64 / n;
    let z2 = Z_975 * Z_975;
    let denominator = 1.0 + z2 / n;
    let center = (estimate + z2 / (2.0 * n)) / denominator;
    let half_width =
        Z_975 * (estimate * (1.0 - estimate) / n + z2 / (4.0 * n * n)).sqrt() / denominator;
    Ok(ProportionInterval {
        successes,
        total,
        estimate,
        lower: (center - half_width).max(0.0),
        upper: (center + half_width).min(1.0),
        confidence_level: 0.95,
        method: "Wilson score interval; nominal 95%",
        assumption: "Independent Bernoulli observations with one common success probability; not established by a selected annotation fixture.",
    })
}

#[derive(Clone, Debug, Serialize)]
pub struct PairedComparison {
    pub total: usize,
    pub both_correct: usize,
    pub left_only_correct: usize,
    pub right_only_correct: usize,
    pub both_incorrect: usize,
    pub discordant_pairs: usize,
    pub left_accuracy: ProportionInterval,
    pub right_accuracy: ProportionInterval,
    /// Left minus right; not a confidence interval.
    pub accuracy_difference: f64,
    pub two_sided_exact_p_value: f64,
    /// Natural logarithm, retained when the displayed probability underflows.
    pub log_two_sided_exact_p_value: f64,
    pub p_value_underflowed: bool,
    pub test: &'static str,
    pub null_hypothesis: &'static str,
    pub assumptions: &'static str,
    pub interpretation_limit: &'static str,
}

/// Exact conditional McNemar comparison for paired correctness vectors.
///
/// Under the null, left-only successes among discordant pairs follow
/// Binomial(discordant_pairs, 0.5). Two-sided p is twice the smaller inclusive
/// tail, capped at one. This uses the exact finite-sample distribution, evaluated
/// in floating point with log-domain summation, rather than chi-squared inference.
/// Callers must enforce identical ordering, target identities and ground truth.
pub fn paired_comparison(
    left_correct: &[bool],
    right_correct: &[bool],
) -> Result<PairedComparison> {
    ensure!(
        left_correct.len() == right_correct.len(),
        "paired correctness vectors must have equal lengths"
    );
    ensure!(
        !left_correct.is_empty(),
        "a paired comparison requires at least one target"
    );
    let (mut both_correct, mut left_only, mut right_only, mut both_incorrect) = (0, 0, 0, 0);
    for (&left, &right) in left_correct.iter().zip(right_correct) {
        match (left, right) {
            (true, true) => both_correct += 1,
            (true, false) => left_only += 1,
            (false, true) => right_only += 1,
            (false, false) => both_incorrect += 1,
        }
    }
    let discordant_pairs = left_only + right_only;
    ensure!(
        discordant_pairs <= MAX_DISCORDANT_PAIRS,
        "discordant-pair count exceeds exact-test computation budget"
    );
    // With discordant cells differing by at most one, the doubled inclusive
    // lower tail is at least one; return the exact capped result directly.
    let log_p_value = if left_only.abs_diff(right_only) <= 1 {
        0.0
    } else {
        let smaller = left_only.min(right_only);
        let mut log_term = -(discordant_pairs as f64) * std::f64::consts::LN_2;
        let mut log_sum = f64::NEG_INFINITY;
        for k in 0..=smaller {
            log_sum = log_add(log_sum, log_term);
            if k < smaller {
                log_term += ((discordant_pairs - k) as f64).ln() - ((k + 1) as f64).ln();
            }
        }
        (log_sum + std::f64::consts::LN_2).min(0.0)
    };
    let p_value = log_p_value.exp().min(1.0);
    let total = left_correct.len();
    Ok(PairedComparison {
        total,
        both_correct,
        left_only_correct: left_only,
        right_only_correct: right_only,
        both_incorrect,
        discordant_pairs,
        left_accuracy: proportion_interval(both_correct + left_only, total)?,
        right_accuracy: proportion_interval(both_correct + right_only, total)?,
        accuracy_difference: (left_only as f64 - right_only as f64) / total as f64,
        two_sided_exact_p_value: p_value,
        log_two_sided_exact_p_value: log_p_value,
        p_value_underflowed: p_value == 0.0,
        test: "Two-sided exact conditional McNemar/binomial test; twice the smaller inclusive tail, capped at one",
        null_hypothesis: "Among discordant pairs, left-only and right-only correctness are equally likely",
        assumptions: "Paired targets and ground truth are identical; experimental units are independent; comparison was specified without selecting on these outcomes.",
        interpretation_limit: "Descriptive on a selected, previously inspected fixture; no population-level significance or independent biological discovery is established.",
    })
}

fn log_add(left: f64, right: f64) -> f64 {
    if left == f64::NEG_INFINITY {
        return right;
    }
    let larger = left.max(right);
    larger + (left.min(right) - larger).exp().ln_1p()
}

#[derive(Clone, Debug, Serialize)]
pub struct AgreementSummary {
    pub total_targets: usize,
    pub compared_targets: usize,
    pub agreements: usize,
    pub disagreements: usize,
    pub not_compared: usize,
    pub coverage: ProportionInterval,
    /// None if no targets were comparable, rather than zero agreement.
    pub agreement_among_compared: Option<ProportionInterval>,
    pub interpretation_limit: &'static str,
}

/// Separate coverage from agreement. Unsupported or context-incompatible records
/// are not silently counted as successes/failures; callers must disclose reasons.
pub fn agreement_summary(
    total_targets: usize,
    compared_targets: usize,
    agreements: usize,
) -> Result<AgreementSummary> {
    ensure!(
        compared_targets <= total_targets,
        "compared targets exceed total targets"
    );
    ensure!(
        agreements <= compared_targets,
        "agreements exceed compared targets"
    );
    let coverage = proportion_interval(compared_targets, total_targets)?;
    let agreement_among_compared = if compared_targets == 0 {
        None
    } else {
        Some(proportion_interval(agreements, compared_targets)?)
    };
    Ok(AgreementSummary {
        total_targets,
        compared_targets,
        agreements,
        disagreements: compared_targets - agreements,
        not_compared: total_targets - compared_targets,
        coverage,
        agreement_among_compared,
        interpretation_limit: "Agreement with an annotation is consistency, not independent functional validation. Exclusions must be disclosed; score intervals assume independent Bernoulli units.",
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1e-10,
            "{actual} differs from {expected}"
        );
    }

    #[test]
    fn wilson_matches_hand_checked_reference_values() {
        let middle = proportion_interval(5, 10).unwrap();
        close(middle.lower, 0.236_593_090_512_564);
        close(middle.upper, 0.763_406_909_487_436);
        let zero = proportion_interval(0, 10).unwrap();
        close(zero.lower, 0.0);
        close(zero.upper, 0.277_532_799_862_889);
        let all = proportion_interval(10, 10).unwrap();
        close(all.lower, 0.722_467_200_137_111);
        close(all.upper, 1.0);
        let one = proportion_interval(1, 1).unwrap();
        close(one.lower, 0.206_549_314_377_237);
        close(one.upper, 1.0);
    }

    #[test]
    fn intervals_are_bounded_and_contain_the_estimate() {
        for total in 1..100 {
            for successes in 0..=total {
                let interval = proportion_interval(successes, total).unwrap();
                assert!(0.0 <= interval.lower && interval.upper <= 1.0);
                assert!(interval.lower <= interval.estimate + 1e-15);
                assert!(interval.upper + 1e-15 >= interval.estimate);
            }
        }
    }

    #[test]
    fn empty_and_inconsistent_counts_are_errors() {
        assert!(proportion_interval(0, 0).is_err());
        assert!(proportion_interval(2, 1).is_err());
        assert!(paired_comparison(&[], &[]).is_err());
        assert!(paired_comparison(&[true], &[]).is_err());
        assert!(agreement_summary(0, 0, 0).is_err());
        assert!(agreement_summary(3, 4, 0).is_err());
        assert!(agreement_summary(3, 2, 3).is_err());
    }

    #[test]
    fn exact_mcnemar_uses_the_small_discordant_distribution() {
        // Six left-only: 2 * (1 / 2^6).
        let result = paired_comparison(&[true; 6], &[false; 6]).unwrap();
        close(result.two_sided_exact_p_value, 0.03125);
        assert_eq!(result.left_only_correct, 6);
        assert_eq!(result.right_only_correct, 0);
        close(result.accuracy_difference, 1.0);
        // Five left-only, one right-only: 2 * (1 + 6) / 2^6.
        let result = paired_comparison(
            &[true, true, true, true, true, false],
            &[false, false, false, false, false, true],
        )
        .unwrap();
        close(result.two_sided_exact_p_value, 0.21875);
        let four_to_one = paired_comparison(
            &[true, true, true, true, false],
            &[false, false, false, false, true],
        )
        .unwrap();
        close(four_to_one.two_sided_exact_p_value, 0.375);
        let balanced = paired_comparison(&[true, false], &[false, true]).unwrap();
        close(balanced.two_sided_exact_p_value, 1.0);
        let identical = paired_comparison(&[true, false], &[true, false]).unwrap();
        close(identical.two_sided_exact_p_value, 1.0);
        assert_eq!(identical.discordant_pairs, 0);
    }

    #[test]
    fn reversing_models_preserves_p_value_and_reverses_difference() {
        let left = [true, true, true, true, false, true, false, false];
        let right = [true, false, false, false, true, true, false, false];
        let a = paired_comparison(&left, &right).unwrap();
        let b = paired_comparison(&right, &left).unwrap();
        close(a.two_sided_exact_p_value, b.two_sided_exact_p_value);
        close(a.accuracy_difference, -b.accuracy_difference);
        assert_eq!(a.both_correct, 2);
        assert_eq!(a.both_incorrect, 2);
    }

    #[test]
    fn log_probability_survives_underflow_without_fabricating_certainty() {
        let result = paired_comparison(&vec![true; 10_000], &vec![false; 10_000]).unwrap();
        assert!(result.p_value_underflowed);
        assert_eq!(result.two_sided_exact_p_value, 0.0);
        close(
            result.log_two_sided_exact_p_value,
            -9_999.0 * std::f64::consts::LN_2,
        );
        assert!(result.log_two_sided_exact_p_value.is_finite());
    }

    #[test]
    fn missing_coverage_does_not_masquerade_as_zero_agreement() {
        let report = agreement_summary(10, 0, 0).unwrap();
        assert!(report.agreement_among_compared.is_none());
        assert_eq!(report.not_compared, 10);
        let report = agreement_summary(10, 8, 7).unwrap();
        close(report.coverage.estimate, 0.8);
        close(report.agreement_among_compared.unwrap().estimate, 0.875);
        assert_eq!(report.disagreements, 1);
    }
}
