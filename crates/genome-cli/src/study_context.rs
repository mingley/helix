//! Finite, matched-source context diagnostics; no model fitting or causal attribution.
use anyhow::{ensure, Context, Result};
use genome_regulatory::math::regression_metrics;
use genome_regulatory_data::{RegulatoryCohort, RegulatoryRecord};
use serde::Serialize;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;

#[derive(Serialize)]
struct Observation<'a> {
    identity: String,
    assay_id: &'a str,
    source_accession: &'a str,
    coefficient: f64,
}
pub fn quality_passes(row: &RegulatoryRecord) -> bool {
    row.score.is_some_and(f64::is_finite)
        && row.unique_tags.is_some_and(|x| x.is_finite() && x >= 10.0)
        && row.dna_sequences.is_some_and(|x| x.is_finite() && x > 0.0)
        && row.rna_sequences.is_some_and(|x| x.is_finite() && x > 0.0)
}

fn mean(values: &[f64]) -> f64 {
    if values.iter().all(|x| *x == values[0]) {
        values[0]
    } else {
        values.iter().map(|x| x / values.len() as f64).sum()
    }
}

/// The best unconstrained prediction for each observed identical input is its
/// mean label. This bound describes the supplied labels, not noise-free biology.
fn identity_floor(observations: &[Observation<'_>]) -> Result<Value> {
    ensure!(
        !observations.is_empty(),
        "identity floor needs observations"
    );
    let mut groups = BTreeMap::<&str, Vec<&Observation<'_>>>::new();
    for row in observations {
        ensure!(
            row.coefficient.is_finite(),
            "nonfinite identity-floor label"
        );
        groups.entry(&row.identity).or_default().push(row);
    }
    let mut sse = 0.0;
    let mut repeated_rows = 0;
    let mut repeats = Vec::new();
    for (identity, rows) in &groups {
        let ys: Vec<_> = rows.iter().map(|r| r.coefficient).collect();
        let center = mean(&ys);
        let group_sse: f64 = ys.iter().map(|y| (y - center).powi(2)).sum();
        sse += group_sse;
        if rows.len() > 1 {
            repeated_rows += rows.len();
            repeats.push(json!({"identity":identity,"mean_source_coefficient":center,
                "sum_squared_deviations":group_sse,"observations":rows}));
        }
    }
    Ok(json!({
        "observations":observations.len(),"distinct_inputs":groups.len(),
        "repeated_input_observations":repeated_rows,"sum_squared_deviations":sse,
        "minimum_empirical_mse_for_an_arbitrary_context_blind_map":sse/observations.len() as f64,
        "identity_definition":"Exact reference-target sequence SHA256 plus reference-validated genomic SNV identity; assay/cell/time/library identity excluded",
        "repeated_groups":repeats,
        "interpretation":"Finite observed-source squared-loss bound, conditional on identical provided sequence/variant inputs. It includes estimator noise and experimental context; it is not a biological irreducibility claim or a population Bayes-error estimate."
    }))
}

pub fn analyze(cohort: &RegulatoryCohort, protocol: &Value) -> Result<Value> {
    let eligible: Vec<_> = cohort
        .records()
        .iter()
        .filter(|r| quality_passes(r))
        .collect();
    let mut by_assay = BTreeMap::<&str, BTreeMap<&str, &RegulatoryRecord>>::new();
    let mut all_observations = Vec::new();
    for row in &eligible {
        let records = by_assay.entry(row.assay_id.as_str()).or_default();
        ensure!(
            records.insert(&row.variant_key, row).is_none(),
            "duplicate assay-variant input in context analysis"
        );
        all_observations.push(Observation {
            identity: format!(
                "{:x}|{}",
                Sha256::digest(row.target_sequence.as_bytes()),
                row.variant_key
            ),
            assay_id: &row.assay_id,
            source_accession: &row.id,
            coefficient: row.score.unwrap(),
        });
    }
    let declarations = protocol["validation_context_analysis"]["pairs"]
        .as_array()
        .context("missing fixed context-pair declarations")?;
    let mut pairs = Vec::new();
    for declaration in declarations {
        let left = declaration["left"].as_str().context("missing left assay")?;
        let right = declaration["right"]
            .as_str()
            .context("missing right assay")?;
        let a = by_assay
            .get(left)
            .context("left context has no eligible variants")?;
        let b = by_assay
            .get(right)
            .context("right context has no eligible variants")?;
        let mut matches = Vec::new();
        let mut observed = Vec::new();
        let mut predicted = Vec::new();
        let mut observations = Vec::new();
        for (key, l) in a {
            let Some(r) = b.get(key) else { continue };
            ensure!(l.target_sequence == r.target_sequence,
                "matched genomic allele has different reference-target inputs; context-blind bound is inapplicable");
            let y_l = l.score.unwrap();
            let y_r = r.score.unwrap();
            observed.push(y_l);
            predicted.push(y_r);
            let identity = format!("{:x}|{key}", Sha256::digest(l.target_sequence.as_bytes()));
            observations.push(Observation {
                identity: identity.clone(),
                assay_id: left,
                source_accession: &l.id,
                coefficient: y_l,
            });
            observations.push(Observation {
                identity,
                assay_id: right,
                source_accession: &r.id,
                coefficient: y_r,
            });
            matches.push(json!({"variant_key":key,"left_id":l.id,"right_id":r.id,
                "target_position_1based":l.target_position_1based,
                "left_coefficient":y_l,"right_coefficient":y_r,"left_minus_right":y_l-y_r,
                "left_tags":l.unique_tags,"right_tags":r.unique_tags,
                "left_cell":l.cell_context,"right_cell":r.cell_context}));
        }
        ensure!(
            !matches.is_empty(),
            "no exact genomic variants overlap a fixed context pair"
        );
        let metrics = regression_metrics(&observed, &predicted)?;
        let differences: Vec<_> = observed
            .iter()
            .zip(&predicted)
            .map(|(a, b)| a - b)
            .collect();
        let floor = identity_floor(&observations)?;
        let lower = floor["minimum_empirical_mse_for_an_arbitrary_context_blind_map"]
            .as_f64()
            .unwrap();
        ensure!(
            (lower - metrics.mse / 4.0).abs() < 1e-12,
            "paired squared-loss lower-bound identity failed"
        );
        pairs.push(json!({
            "left":left,"right":right,"relationship":declaration["relationship"],
            "matched_variants":matches.len(),"unmatched_left":a.len()-matches.len(),"unmatched_right":b.len()-matches.len(),
            "mse":metrics.mse,"mae":metrics.mae,"pearson":metrics.pearson,"spearman":metrics.spearman,
            "signed_mean_left_minus_right":mean(&differences),
            "finite_context_blind_mse_lower_bound":lower,"bound_observations":2*matches.len(),
            "identical_target_sequence_verified":true,"matches":matches,
            "inference_limit":"Repeated-library, time, cell and perturbation comparisons are source diagnostics. Differences include estimation noise; a technical repeat in another element is not an element-specific causal noise baseline."
        }));
    }
    Ok(json!({
        "schema_version":1,"quality_rule":"Finite score; tags>=10; finite positive DNA/RNA counts; retain partial positions",
        "eligible_source_rows":eligible.len(),"pairs":pairs,
        "deduplicated_all_assay_inventory_bound":identity_floor(&all_observations)?,
        "double_counting_policy":"Each pair is reported separately. The inventory calculation includes every source assay accession at most once; overlapping TERT pairs are not pooled as independent measurements.",
        "score_unit":"squared log2 reporter-effect coefficient for MSE; original log2 coefficient for means/MAE",
        "scientific_scope":"Same-study source estimates across selected elements and contexts, not independent laboratories or the entire human genome."
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn finite_input_identity_floor_matches_hand_calculation() {
        let data = [
            Observation {
                identity: "x".into(),
                assay_id: "a",
                source_accession: "1",
                coefficient: 1.0,
            },
            Observation {
                identity: "x".into(),
                assay_id: "b",
                source_accession: "2",
                coefficient: 3.0,
            },
            Observation {
                identity: "y".into(),
                assay_id: "a",
                source_accession: "3",
                coefficient: 4.0,
            },
        ];
        let r = identity_floor(&data).unwrap();
        assert_eq!(r["distinct_inputs"], 2);
        assert_eq!(r["sum_squared_deviations"], 2.0);
        assert!(
            (r["minimum_empirical_mse_for_an_arbitrary_context_blind_map"]
                .as_f64()
                .unwrap()
                - 2.0 / 3.0)
                .abs()
                < 1e-14
        );
    }
    #[test]
    fn constant_source_labels_have_exactly_zero_floor() {
        let data: Vec<_> = (0..7)
            .map(|_| Observation {
                identity: "x".into(),
                assay_id: "a",
                source_accession: "1",
                coefficient: 0.1,
            })
            .collect();
        let r = identity_floor(&data).unwrap();
        assert_eq!(r["sum_squared_deviations"], 0.0);
    }
}
