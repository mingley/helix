# Research milestones

## Current: numerical prediction passes; source-formula validation does not

The **17-crate workspace** now evaluates matched-context sequence prediction
on 55,276 reference-valid regulatory elements. Context-specific slopes achieve
13.434% lower MSE than additive DNA/context on the fixed published columns.
All 243 independent comparison categories pass, including exact core features,
memberships, fits, 829,140 predictions and declared loss/resampling checks.

The source columns contradict the paper's unweighted ESS formula. Their
weighted rank-two structure is compatible with extra scaling, but historical
normalization remains unresolved. The full workflow preserves that failure
and exits 2 despite complete numerical reproduction. It neither repairs the
targets nor claims raw-activity prediction or a new regulatory mechanism.
Read the [study](../research/context-specificity-study.md),
[semantic review](../research/context-specificity-semantic-review.md) and
[snapshot](../research/results/context-specificity-v2/README.md).

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
  context and 62 legacy equality descriptions retained; the original-source
  deletion crosswalk is documented separately.
- Evaluate computational components of five source-backed hypotheses; retain
  biological extensions as untested.
- Generate deterministic reports, a standalone interactive explorer, a revision
  ledger and scoped Lean properties without equating them with biological truth.

See [the current goal](../research/GOAL.md) and
[research synthesis](../research/hbb-mechanisms.md).

## Implemented: auditable computational vertical slice

- Assembly-indexed variants and reference checks in the Rust workspace.
- Frozen public data, provenance commitments and deterministic locus partitions.
- Typed evidence support, replication contracts and contextual contradictions.
- Bounded deterministic research loop with sparse/mixed-group questions.
- Conditional Lean proofs and concrete probability bounds.
- Frozen predictions, common baseline, calibration and disagreement reports.
- Repeatable artifacts, audit command, pinned setup and CI.

## Implemented: frozen regulatory prediction and source forensics

The earlier regulatory study evaluates six fixed sequence models across nine
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

## Next: recover authoritative normalization or compatible measurements

The raw-activity context v1 stopped before fitting because its specified
activity endpoint was absent. V2 openly changed to the available score columns
under a new frozen protocol; numerical success does not erase either source
failure.

1. Resolve the frozen S9 release against authoritative preprocessing code,
   operation order and version-linked tables. A fitted plane cannot identify
   historical normalization.
2. Recover underlying activity/replicate measurements if their reuse terms and
   source mapping permit it, or select an independently compatible assay.
3. Freeze any new endpoint, source, exclusions, units, groups, baselines and
   uncertainty before new outcome inspection. Preserve the existing targets and
   failed checks; do not project labels to force a formula.
4. Distinguish fixed-column prediction, endpoint validation, unseen-cell
   transfer, source-study validation and independent laboratory replication.
5. Use new directed or endogenous observations to discriminate proposed
   mechanisms. Sequence associations and computation alone are insufficient.

Other assays remain future candidates. The
[BRCA1 acquisition review](../research/independent-assay-plan.md) records its
different endpoint and unresolved historical/current reuse terms; it is not
included or claimed as validated.

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
