//! Sequence-only, reverse-complement-invariant features without fitted choices.
use anyhow::{bail, ensure, Result};
use serde::Serialize;

pub const FEATURE_COUNT: usize = 33;
const TRIMER_COUNT: usize = FEATURE_COUNT - 1;
const TRIMER_CODE_COUNT: usize = 64;
const TRIMER_CLASSES: [u8; TRIMER_CODE_COUNT] = trimer_classes();

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct SequenceFeatures {
    /// Reported for integrity; length is not a predictor in the 33-feature vector.
    pub length: usize,
    pub gc_fraction: f64,
    /// Lexicographic A<C<G<T order of min(kmer, reverse-complement(kmer)).
    pub canonical_trimer_counts: Vec<usize>,
    /// Counts divided by L-2, including every overlapping complete 3-mer.
    pub canonical_trimer_frequencies: Vec<f64>,
}

impl SequenceFeatures {
    /// All 32 canonical 3-mer frequencies in lexicographic order, then GC.
    pub fn vector(&self) -> Vec<f64> {
        let mut values = Vec::with_capacity(FEATURE_COUNT);
        values.extend_from_slice(&self.canonical_trimer_frequencies);
        values.push(self.gc_fraction);
        values
    }
}

pub fn feature_names() -> Vec<String> {
    let mut names = Vec::with_capacity(FEATURE_COUNT);
    for code in 0..TRIMER_CODE_COUNT {
        if code < reverse_complement_code(code) {
            let bases = b"ACGT";
            let trimer = [bases[code >> 4], bases[(code >> 2) & 3], bases[code & 3]];
            let name = std::str::from_utf8(&trimer).expect("ASCII DNA");
            names.push(format!("rc3mer_{name}"));
        }
    }
    names.push("gc_fraction".to_owned());
    names
}

/// Reject ambiguous bases and case changes rather than silently rewriting DNA.
pub fn canonical_reverse_complement(sequence: &str) -> Result<String> {
    validate_dna(sequence)?;
    let reverse = reverse_complement(sequence.as_bytes());
    Ok(sequence.min(reverse.as_str()).to_owned())
}

pub fn sequence_features(sequence: &str) -> Result<SequenceFeatures> {
    if sequence.len() < 3 {
        validate_dna(sequence)?;
        bail!("3-mer features require at least three bases");
    }
    let mut counts = vec![0; TRIMER_COUNT];
    let mut gc = 0;
    let mut trimer = 0;
    for (index, base) in sequence.bytes().enumerate() {
        let code = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => bail!("DNA must contain only uppercase A,C,G,T"),
        };
        gc += usize::from(matches!(code, 1 | 2));
        // Keep the most recent three bases, two bits per base.
        trimer = ((trimer << 2) | code) & (TRIMER_CODE_COUNT - 1);
        if index >= 2 {
            counts[usize::from(TRIMER_CLASSES[trimer])] += 1;
        }
    }

    let denominator = (sequence.len() - 2) as f64;
    Ok(SequenceFeatures {
        length: sequence.len(),
        gc_fraction: gc as f64 / sequence.len() as f64,
        canonical_trimer_frequencies: counts
            .iter()
            .map(|count| *count as f64 / denominator)
            .collect(),
        canonical_trimer_counts: counts,
    })
}

fn validate_dna(sequence: &str) -> Result<()> {
    ensure!(!sequence.is_empty(), "DNA sequence is empty");
    ensure!(
        sequence
            .bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T')),
        "DNA must contain only uppercase A,C,G,T"
    );
    Ok(())
}

fn reverse_complement(sequence: &[u8]) -> String {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => unreachable!("DNA validated before reverse complement"),
        })
        .collect()
}

// A=0, C=1, G=2 and T=3 make packed code order equal to DNA lexicographic order.
const fn reverse_complement_code(code: usize) -> usize {
    ((3 - (code & 3)) << 4) | ((3 - ((code >> 2) & 3)) << 2) | (3 - (code >> 4))
}

const fn trimer_classes() -> [u8; TRIMER_CODE_COUNT] {
    let mut classes = [0; TRIMER_CODE_COUNT];
    let mut class = 0;
    let mut code = 0;
    while code < classes.len() {
        let reverse = reverse_complement_code(code);
        if code < reverse {
            classes[code] = class;
            classes[reverse] = class;
            class += 1;
        }
        code += 1;
    }
    assert!(class as usize == TRIMER_COUNT);
    classes
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    #[test]
    fn canonical_counts_are_exact_and_reverse_complement_invariant() {
        let forward = sequence_features("AAACCC").unwrap();
        let reverse = sequence_features("GGGTTT").unwrap();
        assert_eq!(forward, reverse);
        assert_eq!(forward.gc_fraction, 0.5);
        assert_eq!(forward.canonical_trimer_counts.iter().sum::<usize>(), 4);
        assert_eq!(
            forward
                .canonical_trimer_counts
                .iter()
                .filter(|x| **x == 1)
                .count(),
            4
        );
        assert_eq!(forward.vector().len(), FEATURE_COUNT);
        assert_eq!(feature_names().len(), FEATURE_COUNT);
        assert_eq!(
            forward.canonical_trimer_frequencies.iter().sum::<f64>(),
            1.0
        );
        assert_eq!(canonical_reverse_complement("TTT").unwrap(), "AAA");
    }

    #[test]
    fn overlapping_windows_count_and_length_is_not_a_predictor() {
        let three = sequence_features("AAA").unwrap();
        let five = sequence_features("AAAAA").unwrap();
        assert_eq!(three.canonical_trimer_counts[0], 1);
        assert_eq!(five.canonical_trimer_counts[0], 3);
        assert_eq!(three.vector(), five.vector());
        assert_ne!(three.length, five.length);
        for invalid in ["", "AA", "acgt", "ACNT", "AC G"] {
            assert!(sequence_features(invalid).is_err());
        }
    }

    fn dna_strings(length: u32) -> impl Iterator<Item = String> {
        (0..4usize.pow(length)).map(move |mut value| {
            let mut bases = vec![b'A'; length as usize];
            for base in bases.iter_mut().rev() {
                *base = b"ACGT"[value % 4];
                value /= 4;
            }
            String::from_utf8(bases).unwrap()
        })
    }

    #[test]
    fn packed_counts_match_string_counting_and_lexicographic_order() {
        let trimers: Vec<_> = dna_strings(3)
            .map(|sequence| canonical_reverse_complement(&sequence).unwrap())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        assert_eq!(trimers.len(), TRIMER_COUNT);
        let mut names: Vec<_> = trimers
            .iter()
            .map(|label| format!("rc3mer_{label}"))
            .collect();
        names.push("gc_fraction".to_owned());
        assert_eq!(feature_names(), names);

        for sequence in (3..=7).flat_map(dna_strings) {
            let mut counts = vec![0; trimers.len()];
            for window in sequence.as_bytes().windows(3) {
                let window = std::str::from_utf8(window).unwrap();
                let label = canonical_reverse_complement(window).unwrap();
                counts[trimers.binary_search(&label).unwrap()] += 1;
            }
            let gc = sequence.bytes().filter(|base| matches!(base, b'C' | b'G'));
            let expected = SequenceFeatures {
                length: sequence.len(),
                gc_fraction: gc.count() as f64 / sequence.len() as f64,
                canonical_trimer_frequencies: counts
                    .iter()
                    .map(|count| *count as f64 / (sequence.len() - 2) as f64)
                    .collect(),
                canonical_trimer_counts: counts,
            };
            assert_eq!(sequence_features(&sequence).unwrap(), expected, "{sequence}");
        }
    }

    #[test]
    fn invalid_bases_keep_precedence_over_the_minimum_length_error() {
        for (sequence, message) in [
            ("", "DNA sequence is empty"),
            ("A", "3-mer features require at least three bases"),
            ("AG", "3-mer features require at least three bases"),
            ("N", "DNA must contain only uppercase A,C,G,T"),
            ("é", "DNA must contain only uppercase A,C,G,T"),
            ("AGGé", "DNA must contain only uppercase A,C,G,T"),
            ("ATg", "DNA must contain only uppercase A,C,G,T"),
        ] {
            assert_eq!(sequence_features(sequence).unwrap_err().to_string(), message);
        }
    }
}
