# HBB promoter MPRA: acquisition and descriptive evidence protocol

## Scope and disclosure

Acquire a frozen, source-attributed copy of the HBB promoter reporter assay from
MaveDB score set `urn:mavedb:00000018-a-1` and establish its exact mapping to our
GRCh38 reference. This experiment adds a measured regulatory endpoint beyond
the existing clinical classifications and codon-consequence calculations.

Availability checks have already retrieved the CSV header and a few example
rows. No score-pattern exploration, fitting, genotype-score association,
threshold selection or performance comparison has been performed before this
specification. This is not an untouched or independently preregistered holdout.
The immediate deliverable is validated descriptive evidence, not a predictor.

## Authoritative source and permission gate

- Metadata:
  `https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000018-a-1`
- Score table:
  `https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000018-a-1/scores`
- Study: Kircher et al., Nature Communications (2019),
  DOI `10.1038/s41467-019-11526-w`; PMID `31395865`.
- Retain a source snapshot, retrieval time, byte count, response type and SHA-256
  for every persisted input. Preserve the raw metadata and CSV unchanged.
- Do not persist or redistribute score data unless the score-set data license
  and its scope are verified. A software repository license or a metadata-only
  license is insufficient. If terms conflict, retain only permissible metadata
  and report the unresolved requirement.

## Endpoint and context

MaveDB identifies the assay as saturation-mutagenesis massively parallel reporter
assay of the HBB promoter in human HEL 92.1.7 cells. Mutated regulatory fragments
drive a reporter in a plasmid context. The published pipeline uses a multiple
linear model of DNA and RNA barcode counts.

The primary score is the linear-model coefficient for the log2 variant expression
effect. It is a relative reporter-expression effect, not an endogenous HBB RNA
concentration, hemoglobin abundance, oxygen-binding measurement or clinical
pathogenicity probability. Retain coefficient p-values and counts as provided;
do not convert them into a new binary clinical label.

Metadata describes a minimum of ten unique barcode tags for inclusion. Do not
infer replicate number or independence from unique-tag counts: barcodes,
constructs, transfections and biological replicates are different units. Resolve
the study's replicate structure from the primary methods and record any
unresolved detail rather than inventing it.

## Reference and variant reconciliation

The source target is a 187-nucleotide DNA fragment. Its metadata names hg38
chromosome 11, start `5227022` and end `5227208`; the source's coordinate convention
must be resolved by exact sequence reconciliation, not guessed.

1. Compare the complete source target sequence and its reverse complement with
   the independently frozen GRCh38 HBB interval.
2. Require exactly one matching position/orientation and record a zero-based,
   half-open mapped interval. Keep target-sequence orientation separate from the
   HBB gene's minus-strand transcription orientation.
3. Parse source HGVS in its own target-relative coordinate system. For each
   unambiguous SNV, verify the target reference base and map that allele back to
   GRCh38. If reverse-complemented, complement both alleles.
4. Distinguish single substitutions, deletions, equality/reference rows,
   multi-variants, malformed coordinates and other unsupported descriptions.
   Preserve all raw rows; use explicit exclusion reasons for SNV-only summaries.
5. Verify row count against `numVariants` and reject duplicate accessions,
   duplicate SNV representations, reference mismatches and out-of-bounds rows.

The complete source fragment must map before any gene-context conclusion is
made. Spanning the promoter/5-prime UTR boundary is legitimate and must be
reported, rather than forcing every position into one mechanism category.

## Missingness and descriptive summaries

Before computing any summary, disclose denominator and variant scope:

- All source rows and counts by representation class.
- Numeric, finite-score rows versus missing/non-numeric values, separately from
  unsupported genomic representations.
- Reference-validated SNVs and excluded rows, with reasons.
- Source-provided tag/count availability and any stated filtering rule.
- Number of unique target positions; all alternate alleles at one position share
  local sequence context.

Do not silently drop a row because its coefficient, p-value or read count looks
inconvenient. Do not invent a coefficient cutoff, significance threshold,
calibration target or clinical classification. Scores and provided p-values may
be displayed descriptively with their original units; association testing and
model selection require a later separately specified analysis.

## Follow-on prediction study, not executed here

Any trained regulatory predictor needs a new frozen plan before outcomes are
used. Candidate partitions must group all substitutions at a genomic position
and use contiguous promoter blocks with an exclusion buffer at least as large as
the model's sequence window, rather than random row splitting. With one 187-nt
element, that is a limited within-element test, not cross-gene or
cross-experiment generalization. Independent loci, cell contexts and endogenous
perturbation measurements remain necessary.

Audit overlap with the existing ClinVar variants and any source annotations used
as features. Clinical curation can include this functional evidence; such overlap
cannot be presented as independent confirmation. Current public model training
sets may already contain the assay, so pretrained-model data exposure must be
declared when known and marked unknown otherwise.

The result of this acquisition can establish that an independent experimental
endpoint exists and that its coordinates can be reproduced. It cannot establish
that our mechanism engine predicts regulatory effects or that the human genome
has been comprehensively understood.
