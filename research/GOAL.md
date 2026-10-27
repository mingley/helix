# Research objective

Build an inspectable, source-grounded partial model of how human genomic sequence
produces molecular effects. Extend it through reference-validated representations,
explicit assumptions, measured endpoints and falsifiable questions. Annotation,
logical consistency, clinical classification and molecular measurements remain
different kinds of evidence.

## Completed foundation

The workspace now reconstructs HBB's reference DNA → spliced RNA → translated
protein chain, enumerates a bounded SNV universe, audits public clinical
assertions, maps reporter-assay sequences and places the locus in a frozen
primary-chromosome gene atlas. The mechanism map and five explicit questions are
documented in [the HBB synthesis](hbb-mechanisms.md).

These are computational consistency checks and evidence integration. They do not
establish the downstream effect of every variant or a complete model of the
human genome.

## Current result: the frozen regulatory study is negative

The [study report](regulatory-study.md) tests six fixed models across nine
held-out promoter elements. None of the five nonconstant models beats the
training-mean baseline on primary macro position MSE. HBB's small row-random
full-model advantage reverses under the whole-position, 20-base-guarded test.
No model or feature was selected to replace those negative results.

Source forensics reconcile all 62 HBB database equality records to original
deletions. They produce 53 edited identities before QC, but predefined support
filtering leaves 35 rows and 35 identities. The apparent raw duplicate-label
constraint is not a biological noise floor.

Independent implementations reproduce 131,112 predictions and the primary
losses, fits and controls. The complete frozen comparison still fails:
26 of 240 exact-Spearman checks disagree. A separately specified post-hoc
sensitivity analysis matches all metric pairs at every declared prediction
grid; it changes the rank estimand and is not a retroactive v1 pass.
See the [numerical review](rank-sensitivity-review.md).

The [novelty assessment](regulatory-novelty-review.md) distinguishes established
motif/context biology and known source-problem classes from this analysis's
specific quantification. No new biological mechanism or global-first claim is
established.

## Next falsifiable goal: predict across experimentally defined contexts

Test whether adding explicitly measured context improves regulatory-effect
prediction on independent data beyond context-appropriate mean and sequence-only
baselines. Keep the negative v1 study unchanged.

Before fitting the next round:

1. Freeze a new endpoint, source selection, licensing, sequence mapping and
   experimental-unit split. Prefer independent studies and matched same-element
   measurements over treating nearby alleles or library records as independent.
2. State which context is observed: cell state, construct, time, perturbation
   and available factor measurements. Do not infer a reporter RNA sequence from
   native annotation alone.
3. Fix features, models, baselines, inclusion/QC, uncertainty and reporting
   precision before inspecting evaluation outcomes. A finite rank grid, if used,
   needs a scientific justification and a separately declared estimand.
4. Evaluate the predefined primary error on withheld experimental units and
   retain all comparisons, exclusions and failed controls. A claimed gain must
   survive the declared baseline and uncertainty checks, not one favorable split.
5. Use independent endogenous or other discriminating observations to test any
   proposed biological mechanism. Reporter coefficients alone do not establish
   clinical outcomes or native-chromatin function.

## Standing research rules

- The original ClinVar partition and HBB reporter outcomes are explored
  development evidence, not untouched validation.
- Preserve raw sources and failed attempts. Do not change alleles, orientations,
  hashes, QC or metric definitions to make a result pass.
- Software agreement, same-study culture replication and independent biological
  replication are distinct.
- A synonymous variant, intact core splice motif, unmodeled region or absent
  record does not imply no biological function.
- Expand across genes and regulatory programs through these contracts; the
  amount of annotation loaded is not the amount of biology explained.
