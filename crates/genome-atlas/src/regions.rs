use crate::{merge_owned_intervals, Atlas, Gene};
use anyhow::{Context, Result};

/// Exclusive classes in descending annotation priority across both strands.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
#[repr(usize)]
pub enum AnnotationClass {
    Cds = 0,
    ExonOutsideCds = 1,
    TranscriptSpanOutsideExon = 2,
    OutsideTranscriptSpan = 3,
}

impl AnnotationClass {
    pub const fn index(self) -> usize {
        self as usize
    }
}

impl Atlas {
    /// Partition a primary chromosome into maximal runs of one annotation class.
    ///
    /// Each (exclusive end, class index) starts at the previous end, or zero.
    /// A coordinate covered by several curated models takes the highest priority.
    pub fn annotation_runs(&self, chromosome: &str) -> Result<Vec<(u64, usize)>> {
        let chromosome = self
            .chromosomes()
            .iter()
            .find(|entry| entry.name == chromosome)
            .with_context(|| format!("chromosome absent from the primary atlas: {chromosome}"))?;
        let genes = self
            .genes()
            .iter()
            .filter(|gene| gene.chromosome == chromosome.name);
        classify(chromosome.length, genes)
    }
}

fn classify<'a>(length: u64, genes: impl Iterator<Item = &'a Gene>) -> Result<Vec<(u64, usize)>> {
    let mut cds = Vec::new();
    let mut exons = Vec::new();
    let mut spans = Vec::new();
    for gene in genes {
        cds.extend_from_slice(&gene.cds_intervals);
        exons.extend_from_slice(&gene.exon_intervals);
        spans.extend_from_slice(&gene.transcript_intervals);
    }

    // Within each merged layer only one interval can cover the current position.
    use AnnotationClass::{Cds, ExonOutsideCds, TranscriptSpanOutsideExon};
    let mut layers = [
        (Cds, merge_owned_intervals(cds)?),
        (ExonOutsideCds, merge_owned_intervals(exons)?),
        (TranscriptSpanOutsideExon, merge_owned_intervals(spans)?),
    ]
    .map(|(class, intervals)| (class, intervals.into_iter().peekable()));
    let mut runs: Vec<(u64, usize)> = Vec::new();
    let mut start = 0;
    while start < length {
        let mut end = length;
        let mut class = AnnotationClass::OutsideTranscriptSpan;
        for (candidate, intervals) in &mut layers {
            let _ = intervals.next_if(|interval| interval[1] == start);
            if let Some(&[left, right]) = intervals.peek() {
                if left <= start {
                    end = end.min(right);
                    class = class.min(*candidate);
                } else {
                    end = end.min(left);
                }
            }
        }
        if let Some(last) = runs.last_mut().filter(|last| last.1 == class.index()) {
            last.0 = end;
        } else {
            runs.push((end, class.index()));
        }
        start = end;
    }
    Ok(runs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn locus(
        id: &str,
        strand: &str,
        spans: &[[u64; 2]],
        exons: &[[u64; 2]],
        cds: &[[u64; 2]],
    ) -> Gene {
        Gene {
            gene_id: id.into(),
            symbol: id.into(),
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            strand: strand.into(),
            start: spans[0][0],
            end: spans.last().unwrap()[1],
            transcript_count: spans.len() as u64,
            transcript_intervals: spans.to_vec(),
            exon_intervals: exons.to_vec(),
            cds_intervals: cds.to_vec(),
        }
    }

    #[test]
    fn opposite_strands_use_global_priority_and_preserve_transcript_gaps() {
        let forward = locus(
            "1",
            "+",
            &[[1, 7], [12, 15]],
            &[[1, 4], [6, 7], [12, 15]],
            &[[2, 4], [13, 14]],
        );
        let reverse = locus("2", "-", &[[3, 9]], &[[3, 5], [7, 9]], &[[4, 5], [7, 9]]);
        let expected = [
            (1, 3),
            (2, 1),
            (5, 0),
            (6, 2),
            (7, 1),
            (9, 0),
            (12, 3),
            (13, 1),
            (14, 0),
            (15, 1),
            (17, 3),
        ];
        for order in [[&forward, &reverse], [&reverse, &forward]] {
            assert_eq!(classify(17, order.into_iter()).unwrap(), expected);
        }
        assert_eq!(classify(8, std::iter::empty()).unwrap(), [(8, 3)]);
        assert!(classify(0, std::iter::empty()).unwrap().is_empty());
    }
}
