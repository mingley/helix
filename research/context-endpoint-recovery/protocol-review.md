# Endpoint-recovery protocol adherence review

**Status: primary implementation and empirical artifact review complete;
independent numerical comparison pending.** The reviewer ran all 14 synthetic
tests successfully and verified all 56 primary output files against their
manifest. One numerical edge-case bug was corrected before the empirical run.
The reviewer does not edit active primary code or frozen protocols.

The authority is
[context-endpoint-recovery-v1.json](../../experiments/context-endpoint-recovery-v1.json),
SHA-256 `00175e41cf9fa0cbf177ce8fc7ef143eb508aec8b11d058f7acb4458284dc868`.
The expected implementation is
[audit_recovered_context.py](../../scripts/audit_recovered_context.py).

## Required checks

| Area | Required behavior | Review status |
| --- | --- | --- |
| Inputs and schema | Verify frozen protocol/source commitments; preserve exact member/header/ID binding before quantitative checks | Source/schema/evidence hashes checked; official archive receipt verified |
| Variant enumeration | Exactly 36 distinct tuples: three activity bases × two pre-z-score populations × two SD conventions × three post-steps; keep every nonestimable record | Exact Cartesian tuple set verified; all 36 estimable |
| Replicate aggregation | Require three finite values per element/context; never use partial-row means or impute | Static/synthetic pass |
| Replicate medians | Calculate separately on all finite source values for each cell/replicate before requiring complete triples | Static/synthetic pass |
| Pre-z-score populations | Distinguish complete-three-context source rows from each cell's own finite source rows | Static/missing-context fixture pass |
| SD conventions | Population `ddof=0` and sample `ddof=1`; zero SD or insufficient population is nonestimable | Corrected exact-constant guard; tests pass |
| Row centering | Subtract all-three-cell mean only from complete transformed rows | Static/hand-computed fixture pass |
| Post-z-score populations | Use the full complete source population for that pipeline, before old-S9/admitted intersections | Static/no-renormalization fixture pass |
| Fixed cohorts | Full appropriate source population, exact old/new intersection and fixed old-admitted intersection; no favorable residual filtering | Static pass; empty-pair status clarified |
| Exact arithmetic | Exact Decimal equality only for unchanged source cells; derived-float comparisons have no source decimal lexeme | Direct-versus-derived fixture pass |
| General tolerance | `abs(candidate-source) <= 1e-8 + 1e-8*abs(source)`; no widening or reference-direction ambiguity | Static/reference-direction fixture pass |
| Separate code hypothesis | Every exact S12 ID must have a unique finite three-context S11 ESS match; sample z-score on the entire S12 subset, `ddof=1`, maximum error tolerance `1e-12` | Full 55338-ID result passes; maximum error 5.15e−14 |
| Missingness and duplicates | Complete accounting; conflicting empirical keys block affected comparisons rather than averaging or silently choosing | Static/duplicate/missing-key fixtures pass |
| Output completeness | Every variant, cohort, context, missing-field reason and normalizer population retained | 56 files verified; normalizer counts and all route records checked |
| Scope | No predictive fitting, target replacement, free affine calibration, optimized weights or historical-certification claim from compatibility | Static pass |

## Static findings

The initial implementation accepted a mathematically constant decimal context
as having a tiny positive SD. For seven rows `[0.1, i, 6-i]`, `i=0..6`, the first
column is exactly constant, yet NumPy computes its population SD as roughly
`1.39e-17`. The pre-z-score pipeline and separate sample-SD code hypothesis
therefore incorrectly marked this fixture estimable. The implementer added
explicit exact-all-equal checks before pre/post normalization and the separate
hypothesis, with no epsilon or newly selected threshold. Decimal-constant
regressions now pass. This is an arithmetic implementation correction, not a
protocol or endpoint change.

Two reporting clarifications were requested and implemented before the empirical
run:

- A comparison cohort with zero finite paired cells should be explicitly
  nonestimable. Endpoint-method availability must not be mistaken for an
  agreement result merely because the parent comparison object is present.
- The final output should separately state replicate-aggregate arithmetic and
  published-ESS-formula status, with references to every fixed comparison rather
  than selecting a unique historical reconstruction. The current comparison
  tables contain the needed arithmetic, and the final summary now references
  them in separate status objects.

The reviewed static implementation SHA-256 is
`2472660613a5f63f41394eb91aaf99b8e5606fddeca0b23e2ee30c372c75a4f4`.
The schema-binding SHA-256 is
`a59592816f57f7254b780c83ba3e929e3b38f4fd06b2b196a0cddda1d9bfb98e`.
The reviewer verified the S11 source bytes/hash and all five binding-evidence
hashes. The official archive receipt reports the exact 342234089 bytes,
complete ZIP CRC validation and SHA-256
`aefb1974e9ebe5c196a89a00f4e868c0abaf69e3d3952e79c74bfb42f217f6a4`.
It explicitly treats the matching downstream prediction as a cross-source
check, not a publisher-supplied digest.

## Primary empirical review

The reviewed primary summary is
[runs/recovered-context-primary/summary.json](../../runs/recovered-context-primary/summary.json),
SHA-256 `41fc2c0ca904c8320e0e7921d269754599cdf45107366a50cc310fbb01078b67`.
Its archived implementation matches the static-review hash above. The exact
36 Cartesian tuples are present, all estimable. Every post-normalization count
equals the complete source population for its base, and every complete-population
pre-normalizer uses that same count. Old/admitted comparison intersections are
applied afterward.

| Population/check | Observed result |
| --- | --- |
| Published activity aggregate | 56981 source IDs; 56977 complete three-context rows |
| Three-replicate mean bases | 56981 source IDs; 56941 complete three-context rows |
| Published S11 ESS | 55977 complete IDs; maximum row-sum magnitude 1.02e−14 |
| Exact S11 ESS → S12 subset sample-z-score test | All 55338 IDs / 166014 cells pass `1e-12`; maximum 5.15e−14 |
| Same test on old admitted IDs | All 55276 IDs / 165828 cells pass |
| Fixed activity-to-ESS routes | None of the 36 reconstruct S11 ESS or old S9 across all compared cells within the frozen tolerance |

The published aggregate versus the strict three-replicate mean has 170823
complete paired cells. Only five satisfy the frozen tolerance; maximum absolute
difference is `0.0004999977481365825`, MAE `0.0002503983803933819`. This magnitude
must not be converted into a new rounding tolerance or a passed v1 contract.
No replacement aggregate was selected.

The normalization populations differ: the activity aggregate has 1000 more
complete rows than S11 ESS. Therefore the bounded conclusion is **none of the
36 prespecified routes reconstructs the source scores**, not that every possible
implementation of the paper's prose has been disproved. A different population
requires its own explicitly post-hoc contract. Historical raw-activity-to-ESS
processing remains unidentified.

For the separately frozen activity-v3 gate, the reviewer checked all 55276 fixed
admitted IDs against the unchanged published aggregate: none is absent or
nonfinite in any of the three contexts. Their strict three-replicate comparison
also retains all 55276 IDs. This validates coverage of the literal published
endpoint; it does not remove the aggregation uncertainty or independently
validate a predictive model.

## Mathematical interpretation

Some fixed routes can be observationally equivalent. For example, a common
positive scale introduced by an earlier SD convention can disappear after a
subsequent column standardization. Equivalent variants must remain in the
output. No best match is selected as the uniquely identified history.

The separate S11-to-S12 subset test must not be confused with the 36 method-grid
routes: the grid uses full source normalization populations, whereas the code
hypothesis explicitly uses the complete, exact S12 ID subset. Both are fixed
before the new payload values are inspected.

Source-to-source agreement and reconstruction agreement are also distinct.
Exact source-cell Decimal equality can establish that two published tables carry
the same numbers; a floating reconstruction within tolerance establishes an
arithmetic compatibility result. Neither alone creates an independent biological
validation or producer historical documentation.

## Gated activity v3 design review

The parent separately froze
[context-activity-v3.json](../../experiments/context-activity-v3.json), SHA-256
`64e14037fb1f42ddc1ecaca1d23bb8096fabac1f1e892d105eaa20f525fbbb20`,
before inspecting the new headers/activity values. This is a distinct analysis
contract, not a change to endpoint recovery or a repair of v1. The reviewer
made no JSON edits.

No design blocker was identified **conditional on every source, schema, coverage
and design-identity gate passing**. The literal published aggregate is an
operational endpoint; ambiguous empirical meaning, a different worksheet/header,
missing admitted ID or changed source bytes must stop this version rather than
trigger a convenient mapping or smaller cohort.

The primary C-versus-A gain has a useful exact interpretation. For held-out
element `i`, let `a_i` be mean activity across three cells, `m` the grand
training mean, `m_c` each training cell mean and `f_i=z_i beta`. Then

`mean_c[(y_ic-m_c)^2-(y_ic-m_c-f_i)^2]`
`= (a_i-m)^2-(a_i-m-f_i)^2`.

Thus the primary improvement is precisely the shared-DNA improvement on the
mean of the three published activity scales. It is not a gain created merely by
adding context intercepts. The secondary A-versus-I contrast and balanced
pairwise-contrast identity retain their previous interpretation. These are
correlated structural checks, not independent confirmations.

Implementation recommendations sent before fitting:

1. Reuse only independently validated DNA, groups and outcome-independent
   features. Reconstruct v3 observations directly from the pinned empirical
   activity cells. Do not reuse any ESS-derived response mean, coefficient,
   prediction or loss cache.
2. The v3 names `prediction_abs_rel_tolerance` and `metric_abs_rel_tolerance`
   leave the absolute/relative combination unspecified. Enforce the earlier,
   stricter absolute limits `1e-7` and `1e-8`, respectively; report relative
   differences separately. This avoids silently accepting a looser combined
   tolerance. Keep all failures.
3. Report unresolved replicate/normalization arithmetic before fitting. A
   literal published-aggregate binding does not retroactively satisfy v1's
   specific three-replicate-normalization requirement or identify raw count
   preprocessing.
4. “Common activity” concerns these published, library-normalized reporter
   scales. It is not absolute transcripts per cell, universal expression or a
   causal sequence mechanism. The new endpoint is from the same selected source
   experiment, so independent biological replication remains absent.

If a source normalization difference consists only of a fixed additive offset
per cell, the free context intercepts absorb it: C/A/I losses and the primary
C-versus-A gain are invariant. G/S have only a global intercept, so their
individual absolute losses need not be invariant to different cell offsets.
This is an algebraic property, not evidence that the producer used a particular
offset or that the original preprocessing was prospective. Other
transformations can change the endpoint and context weighting and must not be
silently substituted.

**V3 numerical implementation/results review is pending.** This design review
does not authorize a fit before the frozen recovery and endpoint gates pass.
