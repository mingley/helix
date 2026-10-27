//! Sequence-only, reverse-complement-invariant features without fitted choices.
use anyhow::{ensure, Result};
use serde::Serialize;

pub const FEATURE_COUNT: usize = 33;

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
    let mut names: Vec<_> = canonical_trimers()
        .into_iter()
        .map(|kmer| format!("rc3mer_{kmer}"))
        .collect();
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
    validate_dna(sequence)?;
    ensure!(
        sequence.len() >= 3,
        "3-mer features require at least three bases"
    );
    let trimers = canonical_trimers();
    let mut counts = vec![0; trimers.len()];
    for window in sequence.as_bytes().windows(3) {
        let forward = std::str::from_utf8(window)?;
        let reverse = reverse_complement(window);
        let canonical = forward.min(reverse.as_str());
        let index = trimers
            .binary_search_by(|kmer| kmer.as_str().cmp(canonical))
            .expect("every unambiguous trimer has a canonical class");
        counts[index] += 1;
    }
    let denominator = (sequence.len() - 2) as f64;
    Ok(SequenceFeatures {
        length: sequence.len(),
        gc_fraction: sequence
            .bytes()
            .filter(|base| matches!(base, b'C' | b'G'))
            .count() as f64
            / sequence.len() as f64,
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

fn canonical_trimers() -> Vec<String> {
    let mut trimers = Vec::with_capacity(32);
    for a in b"ACGT" {
        for b in b"ACGT" {
            for c in b"ACGT" {
                let bases = [*a, *b, *c];
                let kmer = std::str::from_utf8(&bases).expect("ASCII DNA");
                if kmer <= reverse_complement(&bases).as_str() {
                    trimers.push(kmer.to_owned());
                }
            }
        }
    }
    trimers
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
