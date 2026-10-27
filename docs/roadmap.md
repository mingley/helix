# Research milestones

## Implemented: reference-grounded HBB research and measured regulatory evidence

- Reconstruct the GRCh38 minus-strand transcript, coding sequence and RefSeq
  protein, cross-check reviewed UniProt, and verify known HbS/HbC controls.
- Enumerate 4,824 local SNVs with coding, canonical splice and boundary context.
- Recompute all 54 frozen coarse consequences; report the 17 coding records
  separately and preserve the limited interpretation of agreement.
- Audit selection, source review, missing conditions, nearby train/test loci
  and the now-exposed clinical benchmark.
- Freeze a primary-chromosome RefSeq atlas with merged coverage and HBB neighbors.
- Integrate a separately produced HBB promoter/5-prime-UTR MPRA: 561 SNV effect
  estimates, exact plus-strand target mapping, typed log2 reporter units, source
  context and 62 equality rows preserved without guessing their meaning.
- Evaluate computational components of five source-backed hypotheses; retain
  biological extensions as untested.
- Generate deterministic reports, a standalone interactive explorer, a revision
  ledger and scoped Lean properties without equating them with biological truth.

See [the current goal](../research/GOAL.md) and
[research synthesis](../research/hbb-mechanisms.md).

## Implemented: auditable computational vertical slice

- Nine-crate Rust workspace with assembly-indexed variants and reference checks.
- Frozen public data, provenance commitments and deterministic locus partitions.
- Typed evidence support, replication contracts and contextual contradictions.
- Bounded deterministic research loop with sparse/mixed-group questions.
- Conditional Lean proofs and concrete probability bounds.
- Frozen predictions, common baseline, calibration and disagreement reports.
- Repeatable artifacts, audit command, pinned setup and CI.

## Implemented: frozen regulatory prediction and source forensics

The sixteen-crate workspace now evaluates six fixed sequence models across nine
held-out promoter elements, guarded HBB positions, HBB transfer and an explicitly
easier row-random diagnostic. A metadata-first cohort includes 15 score sets,
with duplicate/context assays reserved for validation-only comparisons. All
derive from one published study, not independent studies. No fitted model beats
the mean baseline in the primary evaluation.

Independent Rust/Python computation reproduces all 131,112 predictions and
primary errors/controls. The full frozen comparison remains failed on 26 of 240
secondary exact-rank metrics. A separate post-hoc precision-grid supplement
diagnoses numerical tie sensitivity without changing that status.

Original-source reconciliation maps all 62 HBB equality descriptions to author
deletions. Their apparent duplicate-identity prediction constraint disappears
under the preexisting tag-quality rule. This is source-estimate accounting,
not a biological noise floor. See the [documented findings](../research/regulatory-study.md).

## Next: independent functional benchmark and context-aware prediction

Before fitting a further model, freeze a new question that can fail: does a
context-aware model improve over sequence-only and training-mean baselines on
an independently produced assay, with whole elements and experimental contexts
held out? Select endpoint-compatible data and inspect source independence first.
The existing same-study contexts cannot establish this result.

Choose any additional public functional assay and separate data by independent experimental
unit. One possible distinct target is BRCA1 saturation genome editing from
[Findlay et al., Nature 2018](https://doi.org/10.1038/s41586-018-0461-z).
Inspect licensing, reference assembly, isoform and assay endpoint first. The
BRCA1 dataset is not included or claimed as validated; the
[acquisition review](../research/independent-assay-plan.md) documents unresolved
historical/current reuse terms and why the HBB assay was selected instead.

Acceptance criteria:

1. Freeze data commitments, leakage groups, endpoint, baseline and metrics before
   selecting models.
2. Add reference-validated allele mapping with explicit assay units and context.
3. Compare with an appropriate published predictor on identical held-out targets.
4. Report uncertainty, missingness, calibration, annotation sensitivity and every
   failed hypothesis.
5. Require new independent observations before claiming discovery.

## Next: interoperable variants and richer evidence

Add a versioned GA4GH VRS adapter with sequence identity, normalization and
computed identifiers checked against authoritative fixtures. Keep unsupported
representations explicit. Add noodles formats and rust-bio algorithms when a
concrete experiment requires them.

Represent measured effect sizes, units, assay conditions, confidence intervals,
cohorts and individual opposing submissions. Keep disease-specific pathogenicity
separate from molecular function.

## Next: evidence-retrieving LLM agents

Implement `ResearchAgent` with a typed proposal schema, bounded retrieval,
source snapshots, study-level provenance and an isolated holdout store. Add
active retrieval of opposing studies. Persist stochastic outputs and compare
with the deterministic baseline.

Select a provider and credential binding when this work begins. The existing
workflow needs no keys. Do not grant arbitrary shell/Lean execution or treat
generated prose as accepted evidence.

## Later: larger scientific models

Add experimentally grounded causal edges and foundation-model adapters with
independent protocols. Lean remains the proof authority until alternatives'
compatibility and trust assumptions are evaluated. A Rust-native proof kernel is
an integration choice, not a substitute for empirical validation.
