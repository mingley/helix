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

## Completed external transfer and exact model limits

The [exploratory external test](variant-transfer-v2-study.md) uses a separate
Siraj experiment and the five old activity fits without retraining. Across
2,153 complete-case SNVs in 172 components passing the frozen prior-source exclusion guards, fixed I macro MSE is
0.2033556 versus 0.2033201 for zero effect. Improvement is −0.00003545, with
conditional 95% range [−0.00144180, 0.00140062]. The useful-transfer rule fails;
modest gains or harm remain compatible with that range.

All 17 independent numerical comparison categories pass, including exact source
floats, 852,588 feature numerators and 344,000 draw indices. Global grouping
separately passes 4,689 checks. The [strict v1 source gate](variant-transfer-study.md)
remains failed on an undocumented NA. V2 explicitly follows that exposure and
excludes literal NA operationally; it does not infer author missingness semantics
or claim an outcome-naive analysis. All Siraj main-panel values used or internally
decoded are explored evidence for future work.

The [outcome-free mutation analysis](variant-locality-study.md) enumerates
745,052 proposed interior edits that preserve every feature. Twenty-one Lean
theorems show the local-context restriction and zero ideal predicted interaction
between substitutions at distance three or greater. These are model properties,
not measured biological neutrality or proof of why external transfer failed.

The [exhaustive nearby-pair analysis](variant-locality/pair-geometry-v1/README.md)
shows nonzero feature interactions for every interior two-substitution class at
distances one and two. Together with the distant-zero proof, this describes the
interaction range of these particular features. The separate
[collision analysis](variant-locality/collision-v1/README.md) finds 92 nonzero
feature-delta groups shared across different reverse-complement edit orbits;
they cover 4,942,664 of 32,502,288 proposed panel slots. These are exact model
limitations, not measured frequencies or biological effect claims. Independent
full-sequence and saved-coefficient checks reproduce the enumerated results.

## Next falsifiable goal: test a predicted absence of interaction

Find source-compatible four-construct measurements: reference, each of two
single substitutions, and their double substitution. For edits at least three
bases apart, the current model predicts exactly zero ideal epistasis. Test the
measured interaction `activity(double) − activity(single1) − activity(single2)
+ activity(reference)` under a common within-cell normalization.

1. Establish actual inserted DNA, allele identity, common activity units,
   replicate support and terms from metadata before new values. Previously used
   Agarwal, Kircher/HBB and Siraj main-panel outcomes are explored.
2. Freeze all four-member groups, separation, endpoints, source/QC rules,
   baselines, uncertainty and failure conditions before the selected outcomes.
   An unresolved source sign or normalization ambiguity stops dependent work.
3. Distinguish count/replicate uncertainty from reproducible nonadditivity.
   Retain every eligible group and negative result; do not select positive
   interactions as a new untouched validation panel.
4. Require further endogenous or other discriminating evidence before claiming
   a biological regulatory mechanism. Continue the historical ESS-source trace
   as a separate provenance question.

The [study-design note](epistasis-design.md) states a candidate replicate-aware
estimand and its independence assumptions. It is not a frozen source-specific
protocol; source feasibility and provenance must be settled before applying it.
The [Siraj design audit](epistasis-source-audit/README.md), independently checked
against the original XLSX, identifies 7,542 sequence-valid windows in 2,156
pair identifiers. Those are design opportunities, not measured outcomes or
independent biological replicates. The [endpoint review](epistasis/source/siraj-s20-semantics.md)
records unresolved table export, common-library and replicate questions; the
biological test remains open until a suitable endpoint and its limitations are
established. The [prior-art note](epistasis-prior-art.md) describes related
published work and existing exposure.

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
