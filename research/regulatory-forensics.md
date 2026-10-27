# HBB regulatory-assay source forensics

## Confirmed identity reconciliation

All **62** literal equality descriptions in the frozen MaveDB HBB table match
**single-base deletion rows in the author's original deposit**. This conclusion
comes from source identities and exact numeric fields, not effect-size patterns.

The original table contains 623 rows: 561 substitutions and 62 alternate alleles
represented as `-`. Every database row matches exactly one original record by
position/allele, coefficient and all three count columns. P-values agree within
five-decimal rounding. The author's file-format documentation explicitly defines
`-` as a one-base deletion; its application code removes `Alt == "-"` when the
deletion option is disabled.

The exact database import, normalization or migration step responsible remains
**unresolved**. No frozen file or generic HGVS parser was changed. No external
issue or message was posted.

All 62 crosswalks and precise calculations are in
[regulatory-forensics.json](regulatory-forensics.json). The reproducible offline
implementation is
[reconcile_hbb.py](forensic-sources/reconcile_hbb.py).

## Commitments and chronology

| Evidence | Retained identity |
| --- | --- |
| Original author table | [OSF GRCh38_HBB.tsv](https://osf.io/download/zpv8r/) |
| Original project | [10.17605/OSF.IO/75B2M](https://doi.org/10.17605/OSF.IO/75B2M) |
| Table version | 1, created November 26, 2018 |
| Table size | 33,879 bytes |
| Original SHA-256 | `705694a6d45cd61159883250393b165a5813b40cc06e1dcdbe3270faad58af07` |
| Frozen database CSV SHA-256 | `ec8b189f70c740f0ab2a11863a0a464fb421d93fefb9da1240cd3831872f57c6` |
| Author code/documentation commit | `16a48622a20754b16ee356d2c4e214459ffe1476` |
| Publication | Kircher et al.2019, [doi:10.1038/s41467-019-11526-w](https://doi.org/10.1038/s41467-019-11526-w), PMID31395865 |
| Frozen prediction protocol SHA-256 | `c8e97fabef2afe8ccaa5a01a7ef7bb952f9c9f15857e5c540af112a66e78cb58` |

OSF's service-reported SHA-256 equals the downloaded digest. The 2018 version-one
table predates the 2019 MaveDB deposit, without establishing when equality
notation first appeared.

[retrieval.json](forensic-sources/retrieval.json) retains URLs, times, byte counts
and hashes for each downloaded response. Sources total approximately one
megabyte. Article XML, original data, OSF metadata and pinned author documentation
remain available for inspection.

## Matching procedure

The 187-nt target is plus-genomic at zero-based GRCh38 `[5227021,5227208)`.
Source one-based genomic position for target `n.k` is `5227021 + k`.
HBB transcription is minus; no unvalidated reversal of target coordinates is used.

For every database row the verifier:

1. Checks source hashes and OSF's independent service commitment.
2. Checks the original reference allele against the exact target.
3. Maps substitutions by position and alternate base.
4. Tests equality records against the original deletion key.
5. Requires exact decimal score, Tags, DNA and RNA agreement.
6. Requires p-value difference at most `0.000005`.
7. Enforces a bijection covering all 623 rows.

The largest p-value difference is `0.00000498055046681`. There are no unmatched
rows. No subset selection, score-nearest-neighbor matching, fitted relationship
or source-order assumption is used.

| Example | Matched value |
| --- | --- |
| Database | `urn:mavedb:00000018-a-1#1`, `n.1=` |
| Original | chr11:5227022, `G → -` |
| Coefficient | `-0.51` in both |
| Tags / DNA / RNA | `68 / 231 / 538` in both |
| Database p-value | `0.0029` |
| Original p-value | `0.00289870880628956` |

Derived `n.kdel` descriptions are **source-anchored**, not complete HGVS
3-prime-normalized identifiers. This crosswalk must not automatically reinterpret
arbitrary equality notation in unrelated score sets.

## 62 source rows encode only 53 different altered targets

Reconstructing each deletion's complete 186-nt edited target yields **53 distinct
sequences**. Nine pairs delete equivalent bases within repeats:

| One-based target positions | Original coefficients | Original tag counts |
| --- | --- | --- |
| 9, 11 | −1.16, 0.03 | 3, 4977 |
| 83, 85 | 0.06, −0.13 | 3, 1770 |
| 92, 93 | 1.05, 0.05 | 2, 1799 |
| 96, 97 | −0.07, 0.19 | 1, 132 |
| 114, 116 | −0.31, −0.08 | 3, 1844 |
| 137, 138 | −0.86, −0.22 | 4, 336 |
| 152, 153 | 1.12, −0.49 | 3, 304 |
| 157, 158 | −1.45, 0.04 | 3, 98 |
| 186, 187 | −0.43, 0.65 | 1, 14 |

Equivalent edited sequences cannot become different unique sequence inputs merely
because an alignment places a deletion at different positions. Different fitted
coefficients do not prove different intrinsic effects of identical isolated DNA.
Representation, multi-mutant background assignment, estimation and correlated
measurement remain possible causes.

For identical input x with observed coefficients y_i and mean ybar:

```text
sum_i (y_i - f(x))^2
  = sum_i (y_i - ybar)^2 + n * (f(x) - ybar)^2
  >= sum_i (y_i - ybar)^2
```

The nine groups' within-identity sum of squares is **4.480450**:

- Minimum observed-label MSE across 62 raw deletion rows:
  **0.07226532258064516** squared log2-effect units.
- Minimum observed-label MSE across 18 repeated-identity rows:
  **0.24891388888888889**.
- Mean group sum of squares across nine groups:
  **0.49782777777777778**; this is group SSE, not per-row MSE.

These are finite algebraic constraints on deposited **point estimates**. They
are not fitted predictor performance, confidence intervals, a biological noise
floor or the latent physical effects.

### Predefined QC falsifies a biological interpretation of that raw bound

**Every equivalent pair has one row below ten tags and one row at or above ten.**
The primary paper's preexisting minimum-ten-tag rule retains 35 deletion rows
with 35 distinct edited targets. All nine repeated-identity groups vanish; this
particular duplicate-identity MSE constraint becomes **zero**.

That does not establish noise-free qualifying measurements. It demonstrates why
a bound from unfiltered source labels must not be presented as irreducible
biology. Ten tags was specified by the paper and frozen protocol, not chosen
after inspecting these differences.

The fixed predictor uses only 561 substitutions, all of which meet ten tags.
Its selection is unaffected. Future deletion modeling must resolve identities,
apply declared QC and group equivalent sequences before evaluation.

## P-values, tags and the additive library

The original table retains 27 deletion rows below ten tags and no substitution
row below ten. Raw retention does not itself violate the downstream-analysis
threshold; the author application filters support at display/analysis time.

MaveDB has **84 literal zero p-values**; the original table has **none**. They
match five-decimal rounding. Database zero therefore means neither mathematical
zero nor infinite significance, and cannot justify a new significance filter.

Error-prone-PCR constructs can contain several changes. Coefficients were fitted
jointly in an additive model using log RNA and log DNA terms across three
transfections of the same constructed library. They do not isolate every allele
in an independent single-mutant clone.

The primary paper explicitly warns that coefficient p-values are statistically
inflated by tag-independence assumptions, correlated counts and unmodeled
epistasis. The retained XML includes this warning in model-fitting methods,
the minimum-tag rationale in paragraph `Par17` and regression in `Par58`.

## Clinical evidence ancestry

The two selected ClinVar SNVs overlapping the fragment were fetched read-only
as full XML. Versions match the frozen ESummary records:

| Record | Explicit PubMed citations in retained XML |
| --- | --- |
| `VCV004910425.1`, c.−43C>T | PMID25741868 |
| `VCV003051702.2`, g.5227091C>T | None |

Neither explicitly cites Kircher PMID31395865. This concerns **recorded
citations**, not independence of organizations, cohorts, biological knowledge or
uncited sources. Reporter expression remains a different endpoint from aggregate
clinical classification.

## Prior knowledge and novelty limits

Read-only public searches found relevant prior acknowledgment:

- [VEP-bench issue32](https://github.com/Open-Athena/VEP-bench/issues/32),
  created September 2, 2026, explicitly warns that a legacy **F9** MaveDB deletion
  appears as target equality and recommends source/CADD validation.
- [Issue62](https://github.com/Open-Athena/VEP-bench/issues/62) discusses MPRA
  equivalent-allele normalization and source-QC ordering.

Thus a legacy deletion/equality mismatch is **already publicly known**.
This exact HBB crosswalk, 62-to-53 sequence accounting and support-stratified
source inconsistency are reproducible analyses here; their global novelty is
**unestablished**.

The bounded Europe PMC query and its two hits are retained. A quoted-numeric
GitHub query returned unrelated numeric-string hits and cannot establish an
absence of prior reports. No comprehensive novelty claim is made.

## Adversarial review of the frozen prediction protocol

Radius-ten supports and strict greater-than-twenty-base guards produce disjoint
train/test input windows. Whole-position grouping is appropriate. They do not
remove covariance inherited from source coefficients jointly fitted to one
mutagenized library.

Nine leave-element-out targets test transport across loci, but remain one
selected study with differing cell/construct contexts. Fixed models, no
score/p-value/count covariates, no test normalization, and validation-only
context records avoid several easier shortcuts.

Interpretation must retain:

- Tags are support, not independent studies or separate effect estimates.
  Source coefficient covariance is unavailable.
- HBB was already exposed. All its coefficients were read for source matching;
  the equivalent-pair values above were read for this identity calculation.
  This work did not inspect new-promoter outcomes or select predictors.
- Bootstrap intervals condition on overlapping fitted OOF predictions and are
  descriptive, not guaranteed population-confidence intervals.
- Chemistry-stratified rotations break some within-position/library dependence.
  Their tail fractions are falsification controls, not biological p-values.
- Context-blind algebraic bounds mix technical and potential biological
  variation; they do not identify either causal contribution.

## Reproduce and source terms

```bash
python3 research/forensic-sources/reconcile_hbb.py --verify
python3 research/forensic-sources/reconcile_hbb.py --self-test
```

Verification checks the complete derived JSON without changing raw scores.
Four self-tests reject changed numeric fields, wrong positions/reference,
non-bijective matching and failures to account for predefined QC.

The article is CC BY4.0; MaveDB's score-set metadata declare CC0. The original
OSF table is a public author deposit with its historical producer data-use
statement. No separate OSF Creative Commons grant was established, so it is not
relabeled CC0. Pinned author code/documentation are GPL-3.0 reference snapshots
with the license retained separately, not incorporated into Helix's MIT matcher.
