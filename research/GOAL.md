# Research objective

Build an inspectable, source-grounded partial model of how human genomic sequence
produces molecular effects. Extend it through reference-validated representations,
explicit assumptions, measured endpoints and falsifiable questions. Annotation,
logical consistency, clinical classification and molecular measurements remain
different kinds of evidence.

## Current result: predict a measured activity component

The [published-activity benchmark](context-activity-study.md) reduces grouped
held-out MSE by **8.7996%** beyond cell means using 32 reverse-complement-canonical
3mer frequencies and GC fraction. The 55,276 reference-valid 200nt elements form
2,761 groups and were tested in three measured cell contexts. Independent
calculations pass all 258 comparison categories and reproduce 829,140 predictions.
The primary conditional component-resampling interval for the MSE improvement is
0.0598440–0.0687169, on the original published activity scales.

The endpoint is unchanged published library-normalized log2(RNA/DNA), not
absolute transcription. The common component is the mean of these three
published scales. Its predictability does not identify a shared transcription
factor or a causal endogenous mechanism. The panel and experiment are shared
with the earlier specificity study; this is a new endpoint analysis, not a new
independent study.

[Official source recovery](context-endpoint-recovery.md) supplies the previously
missing activity table and validates a precise ESS-table to model-table
standardization relationship. Earlier negative arithmetic tests and historical
normalization uncertainties remain recorded. The old v1 source gate and v2
source-formula rule are not rewritten as successes. The
[identifiability analysis](context-identifiability.md) proves why specificity
alone cannot uniquely reconstruct common activity, even with column moment
constraints.

## Next falsifiable goal: predict perturbations in a separate study

The [external transfer contract](../experiments/variant-transfer-v1.md) now
freezes a test of signed Siraj 2026 reporter effects in K562 and HepG2 using the
five existing activity fits without retraining. Source semantics, exact DNA,
coordinate/reference validation and exclusions from all explored loci must pass
before outcomes can be evaluated. The protocol was committed before acquiring
the new outcome asset; a second exact schema/cohort freeze remains required.

The [outcome-free mutation analysis](variant-locality-study.md) quantifies the
model's limits before that test: 745,052 proposed interior edits preserve every
feature, and distant substitutions have zero ideal predicted interaction. These
facts describe the representation and do not count as measured transfer success.

Test whether a model learned from explicit sequence/context measurements
predicts variant effects at untouched regulatory loci in a genuinely separate
experimental study. Separate composition-based prediction from hypotheses about
specific regulatory interactions using observations that can distinguish them.

1. Establish source independence, assay units, actual inserted DNA, perturbation
   identity, replicate support and licensing from metadata before acquisition.
   Prior HBB/Kircher outcomes and the current Agarwal panel are explored evidence.
2. Freeze a small model and outcome contract, whole-locus/homology groups,
   context-matched baselines and negative controls before inspecting new outcomes.
   Do not select a source or threshold for a favorable effect.
3. Test signed perturbation predictions and uncertainty under those conditions,
   retaining source failures and negative results. Whole-fragment activity gains
   do not establish accurate single-base effect predictions.
4. Require directed perturbations and endogenous or other discriminating evidence
   before claiming a biological regulatory mechanism. Continue tracing the
   unresolved upstream ESS processing as a separate source-provenance question.

## Completed foundation

The workspace now reconstructs HBB's reference DNA → spliced RNA → translated
protein chain, enumerates a bounded SNV universe, audits public clinical
assertions, maps reporter-assay sequences and places the locus in a frozen
primary-chromosome gene atlas. The mechanism map and five explicit questions are
documented in [the HBB synthesis](hbb-mechanisms.md).

These are computational consistency checks and evidence integration. They do not
establish the downstream effect of every variant or a complete model of the
human genome.

## Retained negative regulatory study and earlier limitations

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
