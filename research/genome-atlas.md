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

## Reference sequence composition by annotation

The 25 primary hg38 chromosomes contain 3,088,286,401 positions. The Rust scan
assigned each position once, using the union across both strands and all
93,455 curated transcript models. A position within any coding exon takes
priority over other exons, and any exon takes priority over transcript spans.
The counts use the forward reference sequence. GC is (G + C) / (A + C + G + T).
The N column counts unresolved reference positions.
[Per-chromosome counts and source hashes](atlas-composition/summary.json)

| Annotation class | Reference positions | N | GC among called bases |
| --- | ---: | ---: | ---: |
| Coding sequence | 34,561,429 | 0 | 52.29% |
| Exon outside coding sequence | 58,522,466 | 1 | 45.29% |
| Transcript span outside exons | 1,328,191,871 | 720,427 | 41.18% |
| Outside curated transcript spans | 1,667,010,635 | 149,910,292 | 40.16% |

The same strict ordering of GC fractions holds within all 24 nuclear
chromosomes. The smallest adjacent difference is on chr2, where the transcript
span outside exons exceeds the region outside transcript spans by 0.11436
percentage points. The two intermediate classes have zero width on chrM.
[Exact comparison from raw-source counts](atlas-composition/source-check.json)
[Comparison within masked and unmasked sequence](atlas-soft-mask.md)

Across all primaries, 149,910,292 of the 150,630,720 N positions (99.52%)
fall outside curated transcript spans. The sole exonic N in the UCSC twoBit
is at chr21 position 8,217,207 (1-based), or `[8217206, 8217207)`. It lies
inside two overlapping plus-strand exons, the 45S ribosomal RNA precursor
RNA45SN2 (`NR_146144.1`) and the 28S rRNA RNA28SN2 (`NR_146148.1`).
Both models place it outside coding sequence. The source rows and exon
boundaries are retained for replay.
[Source records for the site](atlas-composition/source-check.json) ·
[NCBI RNA45SN2](https://www.ncbi.nlm.nih.gov/gene/109864279/) ·
[NCBI RNA28SN2](https://www.ncbi.nlm.nih.gov/gene/109864282/)

NCBI's GRCh38.p14 genomic FASTA for chromosome 21 (`NC_000021.9`) also
records a literal N at this position. All 17 bases in the window centered on
it agree between FASTA and twoBit. Across the full 46,709,983-base chromosome,
the files differ at three positions. The FASTA records two M symbols (A or C)
and one R (A or G), and UCSC represents each as N. All three are outside the
curated exon union. The comparison uses the forward sequence with case folded.
[Chromosome FASTA comparison and input hashes](atlas-composition/chr21-fasta-control.json) ·
[Versioned NCBI genomic FASTA](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr21.fna.gz) ·
[UCSC ambiguity-code conventions](https://github.com/ucscGenomeBrowser/kent/blob/fa1ec8f7a7071301796d3050f7bb566679b55bbe/src/hg/makeDb/trackDb/human/hg38/gold.html)

An independent Python implementation read the raw UCSC genePred file,
constructed the transcript, exon and CDS unions, counted bases in each union,
and subtracted cumulative counts. It reproduced all 500 chromosome, class
and nucleotide counts. Three hosted runs of the Rust command took 1.21,
1.23 and 1.25 seconds. Peak memory stayed below 96 MiB. The timer covered
reference hashing, atlas loading and counting after the reference download.
[Independent replay](atlas-composition/verify.py) ·
[1.21 s](https://github.com/mingley/helix/actions/runs/44679046007/job/136097757860) ·
[1.23 s](https://github.com/mingley/helix/actions/runs/44676686407/job/136090076363) ·
[1.25 s](https://github.com/mingley/helix/actions/runs/44677803852/job/136093712546)

The curated chrM transcript spans cover exactly the 11,341 positions in the
GenBank CDS union. Of the 5,228 positions outside those spans, 4,017 (76.84%)
belong to GenBank RNA features on NC_012920.1. The 22 tRNAs cover 1,504
positions and the two rRNAs cover 2,513. The control checks interval membership
base by base and matches all 16,569 nucleotides between the GenBank record and
the retained UCSC sequence. This quantifies one omission of the curated track.
[Mitochondrial source control](atlas-composition/mitochondrial-control.json) ·
[Offline replay](atlas-composition/mitochondrial_control.py)

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
haplotype. The raw annotation has models on 438 non-primary contigs. The hg38 twoBit
image has 455 contigs, of which the sequence scan counts the 25 primaries. Duplicated loci, structural variation, population
variation, chromatin state, tissue-specific expression, environment and
interactions remain outside the atlas. The curated mitochondrial source contains
13 protein models but omits its RNA genes. Annotation absence is an observation
about the source set, not an observation that a biological component is absent.

The appropriate next output is therefore a set of reference-grounded mechanistic
claims, explicit contradictions, and experiments capable of falsifying them.
