//! Forward-reference nucleotide composition across the curated primary atlas.
use crate::{AnnotationClass, Atlas};
use anyhow::{ensure, Context, Result};
use genome_data::twobit::{BaseCounts, TwoBitReference};
use serde::Serialize;
use std::{
    collections::BTreeMap,
    io::{Read, Seek},
};

#[derive(Debug, Serialize)]
pub struct AtlasIdentity {
    pub assembly: String,
    pub annotation: String,
    pub coordinate_system: String,
    pub manifest_sha256: String,
    pub catalog_sha256: String,
}

/// All positions, including N, are assigned to exactly one of these classes.
#[derive(Debug, Serialize)]
pub struct ClassComposition {
    pub cds: BaseCounts,
    pub exon_outside_cds: BaseCounts,
    pub transcript_span_outside_exon: BaseCounts,
    pub outside_transcript_span: BaseCounts,
}

impl From<[BaseCounts; 4]> for ClassComposition {
    fn from(counts: [BaseCounts; 4]) -> Self {
        Self {
            cds: counts[AnnotationClass::Cds.index()],
            exon_outside_cds: counts[AnnotationClass::ExonOutsideCds.index()],
            transcript_span_outside_exon: counts[AnnotationClass::TranscriptSpanOutsideExon.index()],
            outside_transcript_span: counts[AnnotationClass::OutsideTranscriptSpan.index()],
        }
    }
}

#[derive(Debug, Serialize)]
pub struct ChromosomeComposition {
    pub name: String,
    pub length: u64,
    pub classes: ClassComposition,
}

#[derive(Debug, Serialize)]
pub struct AtlasComposition {
    pub atlas: AtlasIdentity,
    pub chromosomes: Vec<ChromosomeComposition>,
    pub total: ClassComposition,
}

impl Atlas {
    /// Count A/C/G/T/N on each primary chromosome in annotation priority order.
    /// The caller supplies a reference in the atlas coordinate system. Literal
    /// contig names and lengths are checked before the first packed DNA read.
    pub fn sequence_composition<R: Read + Seek>(
        &self,
        reference: &mut TwoBitReference<R>,
    ) -> Result<AtlasComposition> {
        let sizes: BTreeMap<_, _> = reference.contigs().map(|c| (c.name, c.length)).collect();
        for chromosome in self.chromosomes() {
            let name = &chromosome.name;
            let length = sizes
                .get(name.as_str())
                .with_context(|| format!("reference lacks primary atlas contig {name}"))?;
            ensure!(
                *length == chromosome.length,
                "reference length {length} for {name} differs from atlas length {}",
                chromosome.length
            );
        }

        let mut chromosomes = Vec::with_capacity(self.chromosomes().len());
        let mut total = [BaseCounts::default(); 4];
        for chromosome in self.chromosomes() {
            let runs = self.annotation_runs(&chromosome.name)?;
            let counts = reference.partition_base_counts::<4>(&chromosome.name, &runs)?;
            for (sum, current) in total.iter_mut().zip(counts) {
                *sum += current;
            }
            chromosomes.push(ChromosomeComposition {
                name: chromosome.name.clone(),
                length: chromosome.length,
                classes: counts.into(),
            });
        }
        let catalog_sha256 = self
            .manifest
            .files
            .get("genes.json")
            .context("atlas manifest lacks normalized gene catalog")?;
        Ok(AtlasComposition {
            atlas: AtlasIdentity {
                assembly: self.manifest.assembly.clone(),
                annotation: self.manifest.annotation.clone(),
                coordinate_system: self.manifest.coordinate_system.clone(),
                manifest_sha256: self.manifest_sha256().to_owned(),
                catalog_sha256: catalog_sha256.clone(),
            },
            chromosomes,
            total: total.into(),
        })
    }
}
