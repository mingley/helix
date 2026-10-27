//! Fixed-protocol conditional resampling and multiplicity diagnostics.
use anyhow::{ensure, Result};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

pub const COMPARISON_IDS: [&str; 3] = [
    "delta3mer_vs_chemistry",
    "motif_vs_chemistry",
    "full_vs_coordinate",
];

/// The exact wrapping-unsigned stream fixed by regulatory-prediction-v1.
/// Modulo indexing follows the protocol; it is not a cryptographic sampler.
#[derive(Clone, Debug)]
pub struct SplitMix64 {
    state: u64,
}
impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    pub fn index(&mut self, upper_exclusive: usize) -> Result<usize> {
        ensure!(upper_exclusive > 0, "cannot sample an empty index set");
        let upper = u64::try_from(upper_exclusive)?;
        Ok((self.next_u64() % upper) as usize)
    }
}

/// Each entry is a position-level mean of paired squared-error differences,
/// ordered as COMPARISON_IDS. Average alleles before constructing this value.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct PositionLosses {
    pub position_1based: usize,
    pub loss_differences: [f64; 3],
}

#[derive(Clone, Debug, Serialize)]
pub struct BootstrapComparison {
    pub comparison_id: &'static str,
    pub observed_mean_loss_improvement: f64,
    pub percentile_95_lower: f64,
    pub percentile_95_upper: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct BootstrapSummary {
    pub unit_kind: &'static str,
    pub unit_count: usize,
    pub original_position_or_element_count: usize,
    pub replicates: usize,
    pub seed: u64,
    pub comparisons: Vec<BootstrapComparison>,
    /// Identical resampling draws are used for all three comparisons.
    pub resampling_indices: Vec<Vec<usize>>,
    /// Each row uses the exact comparison order in COMPARISON_IDS.
    pub sampled_statistics: Vec<[f64; 3]>,
    pub quantile_method: &'static str,
    pub interpretation_limit: &'static str,
}

/// Physical blocks are anchored at position1, not at the first observed site.
/// Empty physical blocks are omitted. A shorter/partly observed block remains
/// shorter; each draw concatenates all positions in the selected blocks, so
/// statistics weight positions equally rather than weighting block means equally.
pub fn paired_block_bootstrap(
    positions: &[PositionLosses],
    block_bp: usize,
    replicates: usize,
    seed: u64,
) -> Result<BootstrapSummary> {
    ensure!(block_bp > 0, "physical block width must be positive");
    let mut seen = BTreeSet::new();
    let mut blocks = BTreeMap::<usize, BTreeMap<usize, [f64; 3]>>::new();
    for position in positions {
        ensure!(
            position.position_1based > 0 && seen.insert(position.position_1based),
            "position IDs must be positive and unique"
        );
        validate_values(&position.loss_differences)?;
        let block = (position.position_1based - 1) / block_bp;
        blocks
            .entry(block)
            .or_default()
            .insert(position.position_1based, position.loss_differences);
    }
    let units: Vec<Vec<_>> = blocks
        .into_values()
        .map(|block| block.into_values().collect())
        .collect();
    bootstrap(
        &units,
        replicates,
        seed,
        "Physical position blocks; all included positions in each sampled block",
    )
}

/// Each array contains the three paired per-assay MSE differences, with assays
/// in lexicographic full-URN order. The protocol's primary analysis supplies nine
/// elements. This generic helper also supports small hand-computable fixtures.
pub fn paired_element_bootstrap(
    element_loss_differences: &[[f64; 3]],
    replicates: usize,
    seed: u64,
) -> Result<BootstrapSummary> {
    let units: Vec<_> = element_loss_differences
        .iter()
        .map(|values| vec![*values])
        .collect();
    bootstrap(
        &units,
        replicates,
        seed,
        "Primary promoter elements; each selected element has equal weight",
    )
}

fn validate_values(values: &[f64; 3]) -> Result<()> {
    ensure!(
        values.iter().all(|v| v.is_finite()),
        "all paired loss differences must be finite"
    );
    Ok(())
}

fn bootstrap(
    units: &[Vec<[f64; 3]>],
    replicates: usize,
    seed: u64,
    unit_kind: &'static str,
) -> Result<BootstrapSummary> {
    ensure!(!units.is_empty(), "bootstrap requires observed units");
    ensure!(
        (1..=100_000).contains(&replicates),
        "bootstrap replicate count is empty or exceeds computation budget"
    );
    let mut unit_sums = Vec::with_capacity(units.len());
    let mut observed_sum = [0.0; 3];
    let mut observed_count = 0usize;
    for unit in units {
        ensure!(!unit.is_empty(), "bootstrap units cannot be empty");
        let mut sum = [0.0; 3];
        for values in unit {
            validate_values(values)?;
            for j in 0..3 {
                sum[j] += values[j];
                observed_sum[j] += values[j];
            }
            observed_count += 1;
        }
        validate_values(&sum)?;
        unit_sums.push(sum);
    }
    validate_values(&observed_sum)?;
    let observed = observed_sum.map(|sum| sum / observed_count as f64);
    let mut rng = SplitMix64::new(seed);
    let mut sampled_statistics = Vec::with_capacity(replicates);
    let mut resampling_indices = Vec::with_capacity(replicates);
    for _ in 0..replicates {
        let mut sum = [0.0; 3];
        let mut count = 0usize;
        let mut indices = Vec::with_capacity(units.len());
        for _ in 0..units.len() {
            let index = rng.index(units.len())?;
            indices.push(index);
            count += units[index].len();
            for (j, value) in sum.iter_mut().enumerate() {
                *value += unit_sums[index][j];
            }
        }
        validate_values(&sum)?;
        sampled_statistics.push(sum.map(|value| value / count as f64));
        resampling_indices.push(indices);
    }
    let mut comparisons = Vec::with_capacity(3);
    for j in 0..3 {
        let mut statistics: Vec<_> = sampled_statistics.iter().map(|v| v[j]).collect();
        statistics.sort_by(f64::total_cmp);
        comparisons.push(BootstrapComparison {
            comparison_id: COMPARISON_IDS[j],
            observed_mean_loss_improvement: observed[j],
            percentile_95_lower: quantile_type7(&statistics, 0.025),
            percentile_95_upper: quantile_type7(&statistics, 0.975),
        });
    }
    Ok(BootstrapSummary {
        unit_kind,
        unit_count: units.len(),
        original_position_or_element_count: observed_count,
        replicates,
        seed,
        comparisons,
        resampling_indices,
        sampled_statistics,
        quantile_method: "Type7 linear interpolation on sorted values, q=0.025 and0.975",
        interpretation_limit: "Descriptive conditional resampling of frozen OOF losses, not guaranteed population coverage, simultaneous intervals, model-selection uncertainty or independent-study replication.",
    })
}

fn quantile_type7(sorted: &[f64], q: f64) -> f64 {
    let location = (sorted.len() - 1) as f64 * q;
    let lower = location.floor() as usize;
    let upper = location.ceil() as usize;
    let fraction = location - lower as f64;
    sorted[lower] + fraction * (sorted[upper] - sorted[lower])
}

fn validate_tails(tails: &[f64]) -> Result<()> {
    ensure!(!tails.is_empty(), "tail adjustment needs comparisons");
    ensure!(
        tails
            .iter()
            .all(|v| v.is_finite() && (0.0..=1.0).contains(v)),
        "tail fractions must be finite and between0 and1"
    );
    Ok(())
}

/// Holm step-down arithmetic applied to diagnostic tail fractions. Multiplicity
/// adjustment does not turn an invalid/exploratory null into calibrated inference.
pub fn holm_adjust(tails: &[f64]) -> Result<Vec<f64>> {
    validate_tails(tails)?;
    let mut order: Vec<_> = tails.iter().copied().enumerate().collect();
    order.sort_by(|(ai, a), (bi, b)| a.total_cmp(b).then(ai.cmp(bi)));
    let mut adjusted = vec![0.0; tails.len()];
    let mut previous: f64 = 0.0;
    for (rank, (original, tail)) in order.into_iter().enumerate() {
        previous = previous.max((tail * (tails.len() - rank) as f64).min(1.0));
        adjusted[original] = previous;
    }
    Ok(adjusted)
}

pub fn bonferroni_adjust(tails: &[f64]) -> Result<Vec<f64>> {
    validate_tails(tails)?;
    Ok(tails
        .iter()
        .map(|tail| (tail * tails.len() as f64).min(1.0))
        .collect())
}

#[derive(Clone, Debug, Serialize)]
pub struct FamilyTailAdjustment {
    pub raw: [f64; 3],
    pub holm_within_family: [f64; 3],
    pub bonferroni_across_both_families: [f64; 3],
}

#[derive(Clone, Debug, Serialize)]
pub struct AdjustedControlTails {
    pub comparison_ids: [&'static str; 3],
    pub guarded_hbb: FamilyTailAdjustment,
    pub leave_element_out: FamilyTailAdjustment,
    pub interpretation_limit: &'static str,
}

pub fn adjust_control_tails(
    guarded_hbb: [f64; 3],
    leave_element_out: [f64; 3],
) -> Result<AdjustedControlTails> {
    let both = [
        guarded_hbb[0],
        guarded_hbb[1],
        guarded_hbb[2],
        leave_element_out[0],
        leave_element_out[1],
        leave_element_out[2],
    ];
    let bonferroni = bonferroni_adjust(&both)?;
    let hbb_holm = holm_adjust(&guarded_hbb)?;
    let external_holm = holm_adjust(&leave_element_out)?;
    Ok(AdjustedControlTails {
        comparison_ids: COMPARISON_IDS,
        guarded_hbb: FamilyTailAdjustment {
            raw: guarded_hbb,
            holm_within_family: [hbb_holm[0], hbb_holm[1], hbb_holm[2]],
            bonferroni_across_both_families: [bonferroni[0], bonferroni[1], bonferroni[2]],
        },
        leave_element_out: FamilyTailAdjustment {
            raw: leave_element_out,
            holm_within_family: [external_holm[0], external_holm[1], external_holm[2]],
            bonferroni_across_both_families: [bonferroni[3], bonferroni[4], bonferroni[5]],
        },
        interpretation_limit: "Arithmetic multiplicity controls over six planned spatial-null tails; the artificial rotation distribution is not established as calibrated biological inference.",
    })
}

pub fn empirical_rotation_tail(observed: f64, null_statistics: &[f64]) -> Result<f64> {
    ensure!(
        observed.is_finite()
            && !null_statistics.is_empty()
            && null_statistics.iter().all(|v| v.is_finite()),
        "observed/null statistics must be finite and nonempty"
    );
    let extreme = null_statistics.iter().filter(|v| **v >= observed).count();
    Ok((1 + extreme) as f64 / (1 + null_statistics.len()) as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splitmix_stream_matches_independent_integer_fixtures() {
        let mut rng = SplitMix64::new(0);
        assert_eq!(rng.next_u64(), 0xe220_a839_7b1d_cdaf);
        assert_eq!(rng.next_u64(), 0x6e78_9e6a_a1b9_65f4);
        assert_eq!(rng.next_u64(), 0x06c4_5d18_8009_454f);
        let mut rng = SplitMix64::new(2_026_102_702);
        let indices: Vec<_> = (0..6).map(|_| rng.index(9).unwrap()).collect();
        assert_eq!(indices, [7, 0, 6, 4, 7, 7]);
        assert!(rng.index(0).is_err());
    }

    #[test]
    fn common_element_draws_preserve_three_paired_comparisons() {
        let report = paired_element_bootstrap(&[[1.0, 2.0, 3.0], [5.0, 10.0, 15.0]], 4, 0).unwrap();
        // SplitMix64 seed0 modulo2 starts 1,0,1,0,1,0,1,0: each draw averages
        // the two observed elements, with exactly the same indices for all arms.
        assert_eq!(report.resampling_indices, vec![vec![1, 0]; 4]);
        assert_eq!(report.sampled_statistics, vec![[3.0, 6.0, 9.0]; 4]);
        assert_eq!(report.comparisons[0].observed_mean_loss_improvement, 3.0);
        assert_eq!(report.comparisons[2].percentile_95_lower, 9.0);
        assert_eq!(report.comparisons[2].percentile_95_upper, 9.0);
    }

    #[test]
    fn physical_block_bootstrap_keeps_partial_blocks_and_weights_positions() {
        let positions = [
            PositionLosses {
                position_1based: 1,
                loss_differences: [1.0, 2.0, 3.0],
            },
            PositionLosses {
                position_1based: 2,
                loss_differences: [3.0, 6.0, 9.0],
            },
            PositionLosses {
                position_1based: 22,
                loss_differences: [8.0, 16.0, 24.0],
            },
        ];
        let report = paired_block_bootstrap(&positions, 21, 1, 0).unwrap();
        assert_eq!(report.unit_count, 2);
        assert_eq!(report.resampling_indices, vec![vec![1, 0]]);
        assert_eq!(report.sampled_statistics, vec![[4.0, 8.0, 12.0]]);
        assert_eq!(report.comparisons[0].observed_mean_loss_improvement, 4.0);
        // Equal weighting of block means would incorrectly give(2+8)/2=5.
    }

    #[test]
    fn bootstrap_rejects_empty_duplicate_or_nonfinite_units() {
        assert!(paired_element_bootstrap(&[], 10, 0).is_err());
        assert!(paired_element_bootstrap(&[[0.0; 3]], 0, 0).is_err());
        assert!(paired_element_bootstrap(&[[f64::NAN; 3]], 10, 0).is_err());
        let p = PositionLosses {
            position_1based: 1,
            loss_differences: [0.0; 3],
        };
        assert!(paired_block_bootstrap(&[p, p], 21, 10, 0).is_err());
        assert!(paired_block_bootstrap(&[p], 0, 10, 0).is_err());
    }

    #[test]
    fn type7_quantiles_interpolate_instead_of_selecting_nearest_rank() {
        assert_eq!(quantile_type7(&[0.0, 10.0], 0.025), 0.25);
        assert_eq!(quantile_type7(&[0.0, 10.0], 0.975), 9.75);
        assert_eq!(quantile_type7(&[7.0], 0.025), 7.0);
    }

    #[test]
    fn adjustment_preserves_order_and_enforces_both_planned_families() {
        assert_eq!(
            holm_adjust(&[0.04, 0.01, 0.03]).unwrap(),
            [0.06, 0.03, 0.06]
        );
        let report = adjust_control_tails([0.01, 0.2, 1.0], [0.5, 0.02, 0.4]).unwrap();
        assert_eq!(report.guarded_hbb.holm_within_family, [0.03, 0.4, 1.0]);
        assert_eq!(
            report.guarded_hbb.bonferroni_across_both_families,
            [0.06, 1.0, 1.0]
        );
        assert_eq!(
            report.leave_element_out.bonferroni_across_both_families,
            [1.0, 0.12, 1.0]
        );
        assert!(holm_adjust(&[]).is_err());
        assert!(bonferroni_adjust(&[f64::NAN]).is_err());
        assert!(holm_adjust(&[1.1]).is_err());
    }

    #[test]
    fn rotation_tail_uses_plus_one_and_counts_ties_as_extreme() {
        assert_eq!(empirical_rotation_tail(2.0, &[0.0, 1.0, 2.0]).unwrap(), 0.5);
        assert_eq!(
            empirical_rotation_tail(3.0, &[0.0, 1.0, 2.0]).unwrap(),
            0.25
        );
        assert!(empirical_rotation_tail(0.0, &[]).is_err());
    }
}
