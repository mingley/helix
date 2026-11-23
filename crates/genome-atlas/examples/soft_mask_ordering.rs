//! Compute exact GC orderings from the primary soft-mask atlas report.
use anyhow::{ensure, Context, Result};
use genome_data::twobit::{BaseCounts, MaskedBaseCounts};
use num_rational::BigRational;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{cmp::Ordering, collections::BTreeMap, env, fs::File, io::BufReader};

const CLASSES: [&str; 4] = [
    "cds",
    "exon_outside_cds",
    "transcript_span_outside_exon",
    "outside_transcript_span",
];
const MODES: [&str; 4] = [
    "combined",
    "unmasked",
    "masked",
    "chromosome_common_mask_fraction",
];

#[derive(Deserialize)]
struct Report {
    schema_version: u32,
    operation: String,
    reference: Value,
    atlas: Value,
    chromosomes: Vec<Chromosome>,
}

#[derive(Deserialize)]
struct Chromosome {
    name: String,
    length: u64,
    classes: BTreeMap<String, MaskedBaseCounts>,
}

fn gc(counts: &BaseCounts) -> i128 {
    i128::from(counts.g) + i128::from(counts.c)
}

fn called(counts: &BaseCounts) -> i128 {
    gc(counts) + i128::from(counts.a) + i128::from(counts.t)
}

fn width(counts: &BaseCounts) -> i128 {
    called(counts) + i128::from(counts.n)
}

fn ratio(numerator: i128, denominator: i128) -> BigRational {
    BigRational::new(numerator.into(), denominator.into())
}

fn fraction(numerator: i128, denominator: i128) -> Option<BigRational> {
    (denominator > 0).then(|| ratio(numerator, denominator))
}

fn exact(value: &BigRational) -> Value {
    json!({"numerator": value.numer().to_string(), "denominator": value.denom().to_string()})
}

type Rate = Option<BigRational>;

fn mix(unmasked: &Rate, masked: &Rate, weight: &Rate) -> Rate {
    let w = weight.as_ref()?;
    if *w == ratio(0, 1) {
        return unmasked.clone();
    }
    if *w == ratio(1, 1) {
        return masked.clone();
    }
    Some((ratio(1, 1) - w) * unmasked.as_ref()? + w * masked.as_ref()?)
}

fn difference(first: &Rate, second: &Rate) -> Rate {
    Some(first.as_ref()? - second.as_ref()?)
}

#[derive(Serialize)]
struct PairCounts {
    first_class: &'static str,
    second_class: &'static str,
    greater: usize,
    reversed: Vec<String>,
    tied: Vec<String>,
    undefined: Vec<String>,
}

#[derive(Serialize)]
struct Ranking {
    chromosomes_with_defined_order: usize,
    chromosomes_in_strict_order: usize,
    adjacent: [PairCounts; 3],
}

impl Ranking {
    fn new() -> Self {
        Self {
            chromosomes_with_defined_order: 0,
            chromosomes_in_strict_order: 0,
            adjacent: std::array::from_fn(|i| PairCounts {
                first_class: CLASSES[i],
                second_class: CLASSES[i + 1],
                greater: 0,
                reversed: Vec::new(),
                tied: Vec::new(),
                undefined: Vec::new(),
            }),
        }
    }

    fn record(&mut self, name: &str, rates: &[Rate; 4]) {
        self.chromosomes_with_defined_order += usize::from(rates.iter().all(Option::is_some));
        let mut greater = 0;
        let zero = ratio(0, 1);
        for (pair, rates) in self.adjacent.iter_mut().zip(rates.windows(2)) {
            let gap = difference(&rates[0], &rates[1]);
            match gap.as_ref().map(|value| value.cmp(&zero)) {
                Some(Ordering::Greater) => {
                    pair.greater += 1;
                    greater += 1;
                }
                Some(Ordering::Less) => pair.reversed.push(name.to_owned()),
                Some(Ordering::Equal) => pair.tied.push(name.to_owned()),
                None => pair.undefined.push(name.to_owned()),
            }
        }
        self.chromosomes_in_strict_order += usize::from(greater == self.adjacent.len());
    }
}

#[derive(Serialize)]
struct Contrast {
    chromosome: String,
    first_class: &'static str,
    second_class: &'static str,
}

struct Bound {
    weight: BigRational,
    inclusive: bool,
    limiting_contrasts: Vec<Contrast>,
}

impl Bound {
    fn endpoint(weight: i128) -> Self {
        Self {
            weight: ratio(weight, 1),
            inclusive: true,
            limiting_contrasts: Vec::new(),
        }
    }

    fn tighten(&mut self, weight: BigRational, inclusive: bool, order: Ordering, who: Contrast) {
        let comparison = weight.cmp(&self.weight);
        if comparison == order || (comparison.is_eq() && self.inclusive && !inclusive) {
            self.weight = weight;
            self.inclusive = inclusive;
            self.limiting_contrasts = vec![who];
        } else if comparison.is_eq() && self.inclusive == inclusive {
            self.limiting_contrasts.push(who);
        }
    }

    fn report(&self) -> Value {
        json!({
            "masked_fraction": exact(&self.weight),
            "inclusive": self.inclusive,
            "limiting_contrasts": self.limiting_contrasts,
        })
    }
}

struct WeightRange {
    lower: Bound,
    upper: Bound,
    impossible: Vec<Contrast>,
    crossings: Vec<Value>,
}

impl WeightRange {
    fn new() -> Self {
        Self {
            lower: Bound::endpoint(0),
            upper: Bound::endpoint(1),
            impossible: Vec::new(),
            crossings: Vec::new(),
        }
    }

    fn intersect(&mut self, unmasked_gap: Rate, masked_gap: Rate, who: Contrast) {
        let zero = ratio(0, 1);
        let one = ratio(1, 1);
        match (unmasked_gap, masked_gap) {
            (Some(u), Some(m)) if u > zero && m > zero => {}
            (Some(u), Some(m)) if u > zero || m > zero => {
                // The linear gap is zero at u / (u - m). Strict order excludes the tie.
                let crossing = &u / (&u - &m);
                let above = m > u;
                if zero < crossing && crossing < one {
                    self.crossings.push(json!({
                        "contrast": who,
                        "masked_fraction": exact(&crossing),
                        "strict_above": above,
                    }));
                }
                if above {
                    self.lower.tighten(crossing, false, Ordering::Greater, who);
                } else {
                    self.upper.tighten(crossing, false, Ordering::Less, who);
                }
            }
            (Some(u), None) if u > zero => self.upper.tighten(zero, true, Ordering::Less, who),
            (None, Some(m)) if m > zero => self.lower.tighten(one, true, Ordering::Greater, who),
            _ => self.impossible.push(who),
        }
    }

    fn empty(&self) -> bool {
        let touching = self.lower.weight == self.upper.weight;
        let both_included = self.lower.inclusive && self.upper.inclusive;
        let disjoint = self.lower.weight > self.upper.weight || (touching && !both_included);
        !self.impossible.is_empty() || disjoint
    }

    fn report(&self) -> Value {
        json!({
            "empty": self.empty(),
            "lower_bound": self.lower.report(),
            "upper_bound": self.upper.report(),
            "individually_unsatisfiable": self.impossible,
            "interior_crossings": self.crossings,
        })
    }
}

fn analyze(report: &Report) -> Result<Value> {
    ensure!(
        report.schema_version == 1 && report.operation == "atlas_composition_by_soft_mask",
        "expected a schema v1 soft-mask atlas report"
    );
    for (identity, field) in [
        (&report.reference, "sha256"),
        (&report.atlas, "manifest_sha256"),
        (&report.atlas, "catalog_sha256"),
    ] {
        let hash = identity[field].as_str().context("missing source digest")?;
        ensure!(
            hash.len() == 64 && hash.bytes().all(|b| b.is_ascii_hexdigit()),
            "invalid source digest for {field}"
        );
    }
    let names = (1..=22)
        .map(|i| format!("chr{i}"))
        .chain(["chrX", "chrY", "chrM"].map(str::to_owned));
    ensure!(
        report.chromosomes.len() == 25,
        "expected 25 primary chromosomes"
    );

    let mut rankings: [Ranking; 4] = std::array::from_fn(|_| Ranking::new());
    let mut range = WeightRange::new();
    let mut observed_weights = BTreeMap::new();
    for (chromosome, name) in report.chromosomes.iter().zip(names) {
        ensure!(
            chromosome.name == name,
            "expected {name} in primary chromosome order"
        );
        let classes = &chromosome.classes;
        let has_labels = CLASSES.iter().all(|label| classes.contains_key(*label));
        ensure!(
            classes.len() == CLASSES.len() && has_labels,
            "invalid annotation classes for {name}"
        );
        let strata = CLASSES.map(|label| &classes[label]);
        let width: i128 = strata
            .iter()
            .map(|s| width(&s.unmasked) + width(&s.masked))
            .sum();
        ensure!(
            width == i128::from(chromosome.length),
            "counts do not cover {name}"
        );
        if name == "chrM" {
            continue;
        }

        let masked_called: i128 = strata.iter().map(|s| called(&s.masked)).sum();
        let unmasked_called: i128 = strata.iter().map(|s| called(&s.unmasked)).sum();
        let weight = fraction(masked_called, unmasked_called + masked_called);
        observed_weights.insert(name.clone(), weight.as_ref().map(exact));
        let combined = strata.map(|s| {
            let both_gc = gc(&s.unmasked) + gc(&s.masked);
            fraction(both_gc, called(&s.unmasked) + called(&s.masked))
        });
        let unmasked = strata.map(|s| fraction(gc(&s.unmasked), called(&s.unmasked)));
        let masked = strata.map(|s| fraction(gc(&s.masked), called(&s.masked)));
        let standardized = std::array::from_fn(|i| mix(&unmasked[i], &masked[i], &weight));
        for i in 0..(CLASSES.len() - 1) {
            let who = Contrast {
                chromosome: name.clone(),
                first_class: CLASSES[i],
                second_class: CLASSES[i + 1],
            };
            range.intersect(
                difference(&unmasked[i], &unmasked[i + 1]),
                difference(&masked[i], &masked[i + 1]),
                who,
            );
        }
        let all_rates = [combined, unmasked, masked, standardized];
        for (ranking, rates) in rankings.iter_mut().zip(all_rates) {
            ranking.record(&name, &rates);
        }
    }

    let modes: BTreeMap<_, _> = MODES.into_iter().zip(rankings).collect();
    Ok(json!({
        "schema_version": 1,
        "operation": "soft_mask_gc_ordering",
        "reference": report.reference,
        "atlas": report.atlas,
        "nuclear_chromosomes": 24,
        "class_order": CLASSES,
        "modes": modes,
        "chromosome_masked_called_fraction": observed_weights,
        "common_mask_fraction_range": range.report(),
    }))
}

fn main() -> Result<()> {
    let mut args = env::args_os().skip(1);
    let path = args
        .next()
        .context("usage: soft_mask_ordering REPORT.json")?;
    ensure!(args.next().is_none(), "expected one report path");
    let report: Report = serde_json::from_reader(BufReader::new(File::open(path)?))?;
    println!("{}", serde_json::to_string_pretty(&analyze(&report)?)?);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn who() -> Contrast {
        Contrast {
            chromosome: "chr1".into(),
            first_class: CLASSES[0],
            second_class: CLASSES[1],
        }
    }

    #[test]
    fn opposing_strict_crossovers_exclude_the_touching_point() {
        let u = Some(ratio(1, 5));
        let m = Some(ratio(-1, 10));
        assert_eq!(mix(&u, &m, &Some(ratio(3, 4))), Some(ratio(-1, 40)));
        let mut range = WeightRange::new();
        range.intersect(u, m, who());
        assert!(!range.empty());
        assert_eq!(range.upper.weight, ratio(2, 3));
        assert!(!range.upper.inclusive);
        range.intersect(Some(ratio(-1, 5)), Some(ratio(1, 10)), who());
        assert_eq!(range.lower.weight, range.upper.weight);
        assert!(range.empty());
    }

    #[test]
    fn empty_strata_require_zero_weight_and_ties_fail_strict_order() {
        let positive = Some(ratio(1, 4));
        let zero = Some(ratio(0, 1));
        let one = Some(ratio(1, 1));
        assert_eq!(mix(&positive, &None, &zero), positive);
        assert_eq!(mix(&None, &positive, &one), positive);
        assert_eq!(mix(&positive, &None, &Some(ratio(1, 10))), None);
        let mut range = WeightRange::new();
        range.intersect(positive.clone(), None, who());
        assert!(!range.empty());
        assert!(range.upper.inclusive);
        assert_eq!(range.upper.weight, ratio(0, 1));
        range.intersect(zero.clone(), positive, who());
        assert!(!range.lower.inclusive);
        assert!(range.empty());
        let mut tied = WeightRange::new();
        tied.intersect(zero.clone(), zero, who());
        assert_eq!(tied.impossible.len(), 1);
        assert!(tied.empty());
    }
}
