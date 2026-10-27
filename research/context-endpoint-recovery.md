# Round 3: endpoint recovery and identifiability

## Result and remaining boundary

**A bounded source relationship is recovered numerically; upstream activity-to-ESS
normalization remains unresolved.** The official archive supplies uncentered
activity aggregates and centered ESS. Sample-standardizing that ESS on the
exact S12 subset reproduces every old S9 score within 5.152e-14. All 36 frozen
activity-to-ESS routes fail their tolerance. A separate frozen activity study
now predicts the mean of the published cell scales with **8.7996% lower MSE**
beyond context means and **258/258 independent numerical checks passing**.
That result uses newly recovered observations rather than inversion of ESS;
it does not resolve upstream normalization or establish a new mechanism.

Round 3 separates three questions:

1. Can the permitted source assets and their reuse/byte commitments be verified?
2. Can authoritative version-linked evidence explain how frozen S9 was made?
3. Which activity information would the documented transform retain, even if
   its history were fully resolved?

The existing [v2 benchmark](context-specificity-study.md) is retained unchanged:
13.434% lower held-out MSE on its fixed published columns and complete independent
numerical reproduction, with unresolved scientific source-formula validation.
The earlier [raw-activity v1 source gate](context-v1-source-gate.md) remains
blocked before fitting. This round does not project, rescale or reinterpret
their targets to make either protocol pass.

## Frozen recovery contract and exposure

| Item | Recorded evidence |
| --- | --- |
| Recovery authority | [context-endpoint-recovery-v1](../experiments/context-endpoint-recovery-v1.json), SHA-256 `00175e41cf9fa0cbf177ce8fc7ef143eb508aec8b11d058f7acb4458284dc868` |
| Separate activity authority | [context-activity-v3](../experiments/context-activity-v3.json), SHA-256 `64e14037fb1f42ddc1ecaca1d23bb8096fabac1f1e892d105eaa20f525fbbb20`, frozen before new headers/activity values were inspected |
| Official archive integrity | [Verification receipt](context-endpoint-recovery/source/archive-validation/verification.json); all thirteen members pass |
| Joint-data member | `jointdata.xlsx`, 53,557,805 bytes, SHA-256 `9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e` |
| Exact activity schema | `all_cell_types_summary`, row 1: `name`, `HepG2 [log2(rna/dna)]`, `K562 [log2(rna/dna)]`, `WTC11 [log2(rna/dna)]`; [header-only receipt](context-endpoint-recovery/source/headers-before-values.json) |
| Complete required-ID coverage | Independent source parser verifies all 55,276 fixed IDs have ordinary finite numeric aggregates in all three cells; no required formula cells |
| Aggregation status | Literal aggregate versus unrounded replicate mean differs at rounding-compatible magnitude; historical normalization still unresolved |

The original v2 columns, zero-sum failure, weighted-plane result and model
performance are already known. Any new normalization hypothesis motivated by
them is post-hoc and must remain labeled as such. Literature/code metadata
inspection is prior exposure, not an untouched holdout. Only public methods, code and source metadata were inspected before the recovery
freeze; new numerical comparisons followed its recorded source/schema gates.

## Completed mathematical and synthetic evidence

Across-context centering `P(x)=x-mean(x)*1` removes the common component:
`P(x+h*1)=P(x)`. Exact equal-centered profiles form common-shift equivalence
classes. Fixed downstream processing cannot recover a component absent from
its input. This is established noninjectivity, formalized for the ideal map in
[ContextIdentifiability.lean](../proofs/ContextIdentifiability.lean), not new
biology or an assertion about the actual producer's code.

![Synthetic centering input-stage example](context-identifiability-figures/identifiability.svg)

The inputs `(2,4,3)` and `(7,9,8)` have different means and identical centered
outputs. This first example concerns the input to centering; arbitrary shifts
are not assumed to preserve a recomputed upstream normalizer.

A stronger four-row counterexample satisfies the upstream moment constraints.
Both distinct panels have zero column means and column sums of squares fourteen,
so the same sample normalization `sqrt(14/3)` applies. Their centered outputs
are identical while their common row means differ.

![Synthetic normalization-controlled counterexample](context-identifiability-figures/normalization-controlled.svg)

The [standalone figure bundle](context-identifiability-figures/README.md)
retains exact CSVs, fixed synthetic inputs, SVG/PNG and implementation/source
hash metadata. These examples neither use S9 data nor identify its history.
The figure manifests have been regenerated and verified against the final
formal source/dependency lock and every retained artifact hash.

## Endpoint hierarchy and interpretation

The recovered **library-normalized log2 RNA/DNA activity** provides a measured
endpoint; it does not calibrate absolute transcripts per cell.
Barcode counts, replicate identities and calibration assumptions are additional
information. Recovering a centered transform's scale alone does not uniquely
recover its removed input-stage common component. The
[critical literature review](context-endpoint-literature.md) gives source
definitions, exact algebra and relevant MPRA/relative-data prior art.

ESS is signed centered data, not automatically a positive composition. The
many-to-one and induced-dependence warnings from compositional literature are
useful analogies; they do not justify logging negative scores, adding arbitrary
pseudocounts or calling them proportions.

## Actual source-recovery findings

The official archive contains a new joint-data workbook with the exact activity
headers. Its design workbook and fold/score workbook are byte-identical to the
old frozen S7 and S9 respectively, according to the
[member identities](context-endpoint-recovery/source/member-source-identities.json).
This is a specific source-lineage result: the old score bytes also occur in the
official final archive. It does not identify the score creation operation or
validate the paper's formula. Mere release change cannot be assumed to explain
away the discrepancy.

| Gate | Required evidence | Current finding |
| --- | --- | --- |
| Asset integrity/access | Permitted source bytes and checksums, or documented access failure | Official archive/member verification recovered |
| Semantic role | Explicit measured activity versus model prediction/ESS, exact IDs and version | Literal measured aggregate fields, source coverage and centered ESS distinguished |
| Historical transform | Authoritative code/version and intermediate/output linkage | Upstream history unresolved; a precise S11-ESS-to-S12 numerical map passes |
| Numerical reproduction | All predeclared candidates/comparisons and mismatches | All 36 fixed activity routes reported and fail; separate fixed code hypothesis passes |
| Independent audit | Independently computed identity/transform checks | 7,008/7,008 comparison checks pass; source lexemes exact and all 36 arrays reproduced |

### Centered final ESS and old model labels are different representations

The [primary recovery summary](results/context-endpoint-recovery-v1/primary/summary.json)
has **55,977 complete S11 ESS triples**. Their maximum absolute sum is
**1.0214e-14**, unlike the old S9/S12 model labels. Restricting S11 ESS to all
**55,338 exact S12 IDs**, then subtracting each column's subset mean and dividing
by its subset sample SD (`ddof=1`), reproduces **166,014/166,014 source cells**
at the separately fixed 1e-12 tolerance, maximum absolute difference
**5.1514e-14**. The admitted 55,276-element subset also passes without selecting
IDs by agreement.

This is a numerical source-representation bridge supported by observed values.
It explains how a zero-sum source representation can yield the weighted-plane
model labels without retaining an unweighted zero sum. It does not certify the
producer's unique historical operation order or recover the discarded common
component. The candidate was fixed before new values were inspected and was
motivated by existing code/metadata, not invented as a globally new method.

### Available activity does not reproduce upstream ESS under the fixed routes

The literal published aggregate has **56,977 complete three-cell IDs**.
The 36 frozen routes vary the permitted activity base, normalization population,
SD convention and post-centering scaling. All are estimable; none reproduces
S11 ESS or old S9 at the declared tolerance. This retains the failed candidates
and limits the source claim to the successful downstream map. Other untested
histories are not ruled out, and no best-looking route replaces the targets.

The complete-three-context comparison between the arithmetic mean of the three
published unrounded replicate values and the literal aggregate uses 56,941 IDs.
Its maximum absolute discrepancy is **0.0004999977** and MAE **0.0002503984**.
These magnitudes are compatible with rounding to 0.001, but are not proof of an
exact historical aggregation/rounding rule. Missingness differs between source
aggregates, replicate summaries, final ESS and the model-label subset; these
populations must not be silently interchanged.

The [independent recovery comparison](results/context-endpoint-recovery-v1/comparison/comparison.json)
passes **7,008/7,008 checks**: 796,645 selected workbook rows/lexemes match
exactly, and all 36 derived arrays reproduce within maximum difference
**4.352e-14**. Two earlier comparison attempts stopped on schema-adapter
assumptions and are retained; they were not scientific data mismatches.
Computational agreement confirms that both implementations reproduce the same
bounded negative routes and successful downstream map, not that the failed
activity-to-ESS formulas are valid.

## Separately frozen activity test

V3 uses the literal published activity aggregates without recentering or
rescaling, with the same 55,276 pre-existing elements, groups, DNA features and
fixed models. Its primary comparison is context mean (C) versus additive
DNA/context (A), asking whether sequence predicts the **mean of the three
published library-normalized activity scales**. This mean is not uniquely
identified by specificity alone. It need not be independent of specificity or
have all its predictive signal removed, especially with unequal upstream cell
gains. It is not absolute expression or universal shared TF strength.
The I-versus-A comparison is secondary.

The independent [Python adapter](../scripts/replicate_context_activity.py)
reuses the inspected independent NumPy fitting/bootstrap functions, independently
parses the raw workbook with openpyxl, and ignores old ESS outcomes and outcome
caches. Seventeen new synthetic controls pass: strict numeric/header/ID binding,
fail-closed missingness/formulas, a hand-calculated common-slope ridge example,
and the C-minus-A versus G-minus-S mean-activity loss identity. No additional
activity transforms or ESS-derived fits are used in the replay.

The mean-activity identity is a correlated algebraic diagnostic, not independent
biological confirmation. V3's target is library-normalized reporter activity,
not absolute transcription or independent-source replication. The authorized
direct-source replay and Rust comparison pass **258/258 categories**, including
829,140 predictions and 1,824,108 exact feature fractions. Context-mean MSE is
**0.729254313**, versus additive DNA/context **0.665082737**, primary gain
**0.064171577**, conditional range **[0.059844002, 0.068716904]**. See the
[activity report](context-activity-study.md) for every model/contrast, exposure,
failed implementation attempt and the generic weighting-note qualification.
Equal-element weighting is the v3 authority.

Access and numerical compatibility are source-evidence checks. Neither certifies
the historical production pipeline. The original activity observations are
retained separately rather than reconstructed by inverting specificity.

## Separate post-hoc ESS-cohort diagnostic

After the activity and ESS populations were found to differ by 1,000 IDs,
we froze the separate
[six-route cohort audit](../experiments/context-ess-cohort-audit-v1.json) before
performing those normalizations. It is post-hoc, distinct from the original
36-route recovery contract, and does not modify v3 activity targets or fits.

All six routes use the full 55,977 ESS-ID cohort with the three unchanged
activity bases and two SD conventions. Every route fails the S11 ESS tolerance;
none of 167,931 cells is within tolerance. S11 MAE is approximately 0.07810,
maximum error 1.57638–1.57690. Their six standardized S12 derivatives likewise
fail the separate 1e-12 criterion, with maximum error approximately 1.20917.
The [recovery snapshot](results/context-endpoint-recovery-v1/README.md) retains
the six-route addendum. All **1,162 independent numerical, ID and metadata
comparison checks pass**, reproducing these six negative results. An earlier
comparison differed only in the reserved label for two unmapped-category IDs;
the explicit label adapter and failed receipt are retained. No empirical
category, cohort, activity or normalization was changed.

The 1,000-ID gap is mixed: enhancer classes contribute 193, 323 and 322 IDs,
promoters 118, synthetic positives 33, shuffled controls six, synthetic negatives
four and K562 positives one. It is not solely the extra variable-promoter class
or a single control class. These bounded negative tests do not prove that every
possible normalization or population could fail. No further transform is
selected to make a result agree.

## Novelty and next discriminating evidence

The general loss of common scale under centering, the need for calibration,
and count/precision-aware MPRA analysis are established mathematics and methods.
The exact source-version crosswalk, explicit numerical discrepancy or bounded
recovery attempt can be useful case-specific computational findings. Global
novelty and new biology remain unestablished.

The exact ESS subset/sample-standardization hypothesis was present in
downstream public code before this investigation. The pinned
[Dalal processing code](https://github.com/TaykhoomDalal/LentiMPRA-Enhancer-Processing/blob/d2751a6ddff67ce4406288e9f6c3d932b0c12cdd/LentiMPRA/process.py#L633-L661)
and reported initial 2026-08-21 commit have the same function AST. See the
[code trace](context-normalization-trace/README.md). That is concrete prior art,
not original-producer certification or an authenticated historical timestamp.
The bounded official-byte comparisons and documented limitations are this
analysis's specific contribution; they are not a global-first claim.

Further work should distinguish an authoritative upstream reconstruction from
mere numerical compatibility, retain all failed routes and preserve the
well-defined literal-activity benchmark. New independent, appropriately
calibrated measurements are needed for broader biological claims. The ideal
identifiability theorem limits inversion of relative scores; it does not make
the measured biological system unknowable.
