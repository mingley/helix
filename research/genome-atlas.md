# Genome context: what the annotation establishes

The first systems-level goal is to connect reference sequence, transcript
architecture, molecular changes, and phenotype while preserving which links are
observed, annotated, inferred or still unknown. A whole-genome annotation atlas
provides a map of where to ask those questions; it is not a completed explanation
of the genome.

## Measured properties of the frozen annotation

The [dataset card](../datasets/genome-atlas/README.md) and its manifest define the
source, inclusion policy and counting units. The 25 primary reference sequences
span 3,088,286,401 bases. Within this curated RefSeq set, merged transcript spans
cover 46.0215%, exons 3.01410%, and coding exon portions 1.11911%.

The gap between transcript-span and exon coverage is largely sequence inside
transcript models that is not exonic in any included isoform. These intervals
include introns; the atlas does not establish their complete functions. The gap
between exon and CDS coverage includes untranslated and noncoding RNA exons.
Sequence outside the curated models can contain regulatory elements, repetitive
sequences, unmodeled transcripts and other functional or unresolved features.
It must not be labeled nonfunctional from this atlas.

The 29,949 gene/chromosome/strand groups represent 29,922 NCBI Gene identities,
including 19,427 protein-coding loci. Gene number is not equivalent to transcript
number, protein number, number of regulatory programs, or an explanation of
phenotype. The 93,455 curated models already demonstrate the importance of
alternative transcript structure, but they are not all possible cell states or
all human transcripts.

## HBB sits in a structured locus

The atlas places HBB, NCBI Gene 3043, on the minus strand of chromosome 11 at
`[5225463, 5227071)` with three exons. The nearest mapped locus is HBD, separated
by 5,766 reference bases between gene envelopes. Other nearby loci include HBBP1
(14,882 bp), BGLT3 (17,482 bp), HBG1 (21,197 bp), HBG2 (26,116 bp), and HBE1
(41,273 bp), all toward higher genomic coordinates. OR51V1 lies 24,781 bp toward
lower coordinates.

Because HBB is on the minus strand, higher genomic coordinate is its transcript
5-prime direction. Genomic left/right and transcriptional upstream/downstream
must be distinguished before making a regulatory or splice-site claim. These
proximity observations alone do not establish enhancer contacts, developmental
switching mechanisms or disease causality.

This locus is a useful mechanistic starting point: an explicit genomic SNV can be
mapped through a checked transcript model to a coding or splice-related change.
That molecular annotation can then be compared with clinical assertions, rather
than predicting clinical labels directly from those same assertions.

## Next questions and discriminating observations

1. Which sequence changes alter the protein, canonical splice signals, or the
   untranslated transcript? Reconstruct the transcript and translation from
   checked reference data; keep all predicted effects distinct from measurements.
2. Where do recomputed molecular consequences differ from database consequence
   labels? Examine transcript choice, coordinate conventions, start/stop treatment
   and evidence provenance before treating mismatches as biological discoveries.
3. Which synonymous, intronic or UTR variants have credible functional evidence?
   These annotation categories alone do not imply neutrality. Relevant follow-up
   measurements include splicing, transcript abundance, translation and protein
   behavior under specified cellular conditions.
4. Which claims depend on developmental stage or erythroid context? Proximity in
   this reference atlas provides candidates, but independent expression and
   regulatory-contact observations are needed to test those relationships.
5. Can a sequence-derived explanation predict an independent measured endpoint?
   Freeze an independent assay protocol before model selection. ClinVar curation
   and a single reference sequence cannot establish that success themselves.

## Known blind spots

The primary reference is a coordinate system, not a representation of every human
haplotype. Excluding 438 non-primary contigs removes alternate/fix/unplaced models
from this counting exercise. Duplicated loci, structural variation, population
variation, chromatin state, tissue-specific expression, environment and
interactions remain outside the atlas. The curated mitochondrial source contains
13 protein models but omits its RNA genes. Annotation absence is an observation
about the source set, not an observation that a biological component is absent.

The appropriate next output is therefore a set of reference-grounded mechanistic
claims, explicit contradictions, and experiments capable of falsifying them.
