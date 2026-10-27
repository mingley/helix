# Research objective

Build an inspectable, source-grounded partial model of how human genomic sequence
produces molecular effects. Extend the model only through reference-validated
representations, explicit assumptions, quantitative observations and falsifiable
questions. Genome annotation, logical consistency and clinical classification
are different kinds of evidence.

## Current milestone: understand a complete, bounded HBB example

1. Place HBB in a frozen annotation of the primary human chromosomes, distinguishing
   genes, transcripts, coding sequence, noncoding annotations and uncovered sequence.
2. Recover its named DNA → spliced RNA → translated protein chain exactly.
3. Recompute all local single-base sequence consequences without using ClinVar
   consequence labels as inputs, and enumerate the complete disclosed SNV universe.
4. Reconcile existing clinical assertions while exposing missing context,
   ascertainment bias, review strength, nearby variants and unestablished independence.
5. Add a separately produced public HBB reporter-expression assay with exact target
   mapping, measured units, cell context and independently checked reuse terms.
6. Evaluate the computational parts of five explicit questions and preserve the
   biological extensions as untested where relevant observations are absent.
7. Save deterministic machine-readable outputs, a readable synthesis, an interactive
   locus explorer, source commitments and scoped Lean proofs; run integrity and
   implementation checks before publishing code.

This milestone is a reference mechanism map and measured-evidence integration,
not reverse engineering of the entire human genome or a new biological discovery.
The reporter assay tests a fragment in a named cellular context; it does not
establish endogenous regulation, organism physiology or clinical outcomes.

## Research decisions

- The original 54 ClinVar records, including the former heldout partition, have
  been explored. They are development evidence now, not an untouched test set.
- The prior 96% versus 84% clinical result does not establish superiority
  (five discordant pairs, exact paired p = 0.375 under unestablished assumptions).
- Prespecified reference and assay-mapping checks can fail. Do not change
  expected alleles, target orientation, source hashes or biology to make them pass.
- Reporter data are not converted into clinical pathogenicity labels. Culture
  replicates from one study are not independent scientific replication.
- An unmodeled region, synonymous change, intact canonical splice motif or absent
  clinical record does not imply lack of biological function.

## Following milestone

Freeze an independent functional-prediction task before fitting a richer model:
define the endpoint, relevant cellular context, genomic/experimental leakage
groups, comparison methods, missingness policy and uncertainties. Validate
sequence-to-assay mapping, then test whether the extra mechanistic features
improve predictions on genuinely withheld units. Expand to additional genes and
regulatory programs through the same contracts rather than claiming whole-genome
understanding from one locus.
