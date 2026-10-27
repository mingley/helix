# Regulatory research: novelty and claim review

## Status and standard

This review distinguishes **new to this workspace** from **new biological
knowledge**. The search is bounded and cannot establish global novelty.
A novel model coefficient, a source representation correction or a rediscovered
motif effect is not by itself a new mechanism.

The reference/assay study is Kircher et al. 2019,
[doi:10.1038/s41467-019-11526-w](https://doi.org/10.1038/s41467-019-11526-w).
The exact study protocol is
[regulatory-prediction-v1.json](../experiments/regulatory-prediction-v1.json).
HBB was already explored; the additional cohort and fixed features were selected
before its new score acquisition. The protocol is retrospective and exploratory,
not preregistered clinical inference.

## Claim ledger

| Candidate claim | Current assessment | What can be stated |
| --- | --- | --- |
| Allele-aware motif changes affect expression | Established methodological and biological prior | Our four-profile feature is a fixed, interpretable implementation, not a newly discovered general mechanism. |
| Motif predictors transfer differently across elements/cell contexts | Established qualitative prior; exact reanalysis results may add dataset-specific quantification | The original study already reported substantial element/context heterogeneity. New metrics must identify their narrower contribution. |
| Legacy MaveDB equality rows can encode source deletions | This problem family was publicly documented before this work | Exact HBB source reconciliation is a reproducible curation result; a global novelty claim is unsupported. |
| HBB's 62 source deletion rows encode 53 distinct edited targets | A verified dataset-specific computational observation here; global novelty unestablished | Complete altered-sequence identity is a stronger grouping key than an unnormalized deletion position. |
| Repeated-identity disagreements imply an irreducible biological noise floor | Rejected by source semantics and predefined QC | A finite bound on raw deposited estimates is not a bound on latent biological effects; the duplicate-identity bound vanishes in the qualifying subset. |
| Motif/full models improve primary held-out regulatory prediction | Unsupported by the frozen primary MSE result | All five nonconstant models are worse than the mean baseline; a tiny motif-versus-chemistry contrast is uncertain and does not establish a new mechanism. |
| Tiny prediction differences cannot affect rank conclusions | Falsified by 26/240 frozen exact-rank comparisons | Primary losses reproduce, while exact ranks remain numerically sensitive at ties. |
| All-grid matching repairs the original study | Unsupported | The post-hoc supplement changes the rank estimand; it does not override frozen failures or establish biology. |

## 1. Motif associations are prior art

[The literature assessment](regulatory-literature.md) records relevant sources
and access scope. Important antecedents include:

- EKLF/KLF1 regulation of beta-globin CACCC elements: Miller/Bieker 1993
  [PMID 7682653](https://pubmed.ncbi.nlm.nih.gov/7682653/) and Marini et al. 2004
  [PMID 15384985](https://pubmed.ncbi.nlm.nih.gov/15384985/).
- Globin CCAAT/NF-Y regulation and developmental context: Martyn et al. 2017
  [PMID 27718361](https://pubmed.ncbi.nlm.nih.gov/27718361/), Fang et al. 2004
  [PMID 14645237](https://pubmed.ncbi.nlm.nih.gov/14645237/).
- Allele-aware motif scoring tools: motifbreakR 2015
  [PMID 26272984](https://pubmed.ncbi.nlm.nih.gov/26272984/) and atSNP 2015
  [PMID 26092860](https://pubmed.ncbi.nlm.nih.gov/26092860/).
- Direct mutant-TERT GABP recruitment: Bell et al. 2015
  [PMID 25977370](https://pubmed.ncbi.nlm.nih.gov/25977370/).
- Kircher 2019's own comparison with JASPAR predictions, discussion of
  reference-scan failure to model gained motifs, cross-element heterogeneity,
  and the SORT1 orientation/minimal-promoter interpretation.

A motif association in this assay can nominate a candidate explanation, but
neither a PFM score nor a fitted regression verifies occupancy or causality.
The four-profile model also omits many factors: for example, it has no
ETS/GABP matrix. Failure to explain TERT's published mechanism with this feature
set would not be a surprising new property of biology.

## 2. A source-forensic observation survives direct reconciliation

The [forensic analysis](regulatory-forensics.md) and its retained source records
show an exact bijection between all 623 database rows and the original author
table. Specifically, every one of the 62 equality-labelled database rows matches
a source single-base deletion by position and exact effect/count fields, with
p-values agreeing within five-decimal rounding.

The original table predates the database deposit, and the producer documentation
defines its deletion symbol explicitly. The analysis does **not** identify the
particular import or migration operation responsible for the database notation.
It therefore supports “original deletions represented as database equalities”
rather than an attribution to an unexamined software component or person.

This result is useful for reproducibility: arbitrary equality notation must not
be globally reinterpreted as deletion, and the original file remains unchanged.
A source-anchored crosswalk can be reviewed separately from a general HGVS
parser.

### Prior public acknowledgment prevents a broad novelty claim

[VEP-bench issue 32](https://github.com/Open-Athena/VEP-bench/issues/32), created
September 2, 2026, explicitly warns of a legacy F9 source deletion represented as
target equality in MaveDB and recommends source validation.
[Issue 62](https://github.com/Open-Athena/VEP-bench/issues/62), created
September 4, 2026, discusses equivalent MPRA allele normalization and QC ordering.
Their public timestamps predate this workspace analysis.

These sources acknowledge the **problem family**; they do not by themselves
establish that every precise HBB count or bound below was previously reported.
Conversely, not locating those exact numbers in this bounded search does not
prove they are globally novel.

## 3. Sequence identity, source estimates and QC

Reconstructing the complete altered HBB target converts 62 source deletion rows
into 53 distinct sequences. Nine pairs are repeated sequence identities with
different deposited effect estimates.

For identical sequence input x, the elementary squared-error identity gives:

```text
sum_i (y_i - f(x))^2
  = sum_i (y_i - mean(y))^2 + n*(f(x)-mean(y))^2
  >= sum_i (y_i - mean(y))^2
```

For these raw deposited rows, the verified within-identity sum of squares is
4.480450, implying an observed-estimate MSE bound of approximately 0.0722653
over the 62 raw deletion rows. This is a finite accounting result. The
identity is standard algebra; the dataset-specific input identities and exact
arithmetic are the contribution of the reconciliation.

Crucially, **each of the nine pairs has one row below ten tags and one at or
above ten**. The study's predefined minimum-ten-tag criterion leaves 35 rows with
35 distinct edited targets, so this particular duplicate-identity bound becomes
zero after QC. This does not show that the retained measurements are noise-free.

The bound therefore cannot support any of these stronger statements:

- The physical response of identical DNA is intrinsically different.
- Sequence alone has an irreducible biological prediction error of 0.0722653.
- Tag count identifies an independent biological replicate.
- The raw-row result applies to the frozen SNV evaluation.
- The remaining context, barcode or library effects have been causally resolved.

The study's coefficients come from a joint additive library model and their
covariance is unavailable. Representation, support, background mutations and
estimation remain relevant. This is an instructive falsification of a tempting
biological interpretation, not evidence for a new biological law.

## 4. Additional source precision result

The database has 84 literal zero p-values while the source has none; exact source
matching shows consistency with five-decimal rounding. “0” in the deposited
database field is not mathematical zero or infinite statistical significance.
This is a source-precision audit, not a new statistical method. The original
paper additionally discusses statistical-inflation limitations of its coefficient
p-values. The frozen predictors do not use them as features.

## 5. Frozen model results: a negative primary benchmark

The [retained study snapshot](results/regulatory-v1/README.md) and
[final independent comparison](independent-regulatory-replication/comparison-v1.json)
are the release evidence. Earlier attempt receipts remain historical records.
The frozen protocol hash is
`c8e97fabef2afe8ccaa5a01a7ef7bb952f9c9f15857e5c540af112a66e78cb58`.

Independent Rust/Python calculations reproduce all 131,112 predictions within
`7.729e-12`, inside the frozen `1e-7` tolerance. Primary errors, fitted
preprocessing, folds, resamples and spatial controls also match at their declared
tolerances. This verifies the computations actually checked, not independent
biological replication.

**The full frozen comparison remains failed:** 26 of 240 secondary raw
exact-Spearman checks disagree beyond `1e-8`; 214 comparisons match, including
36 jointly undefined cases. Corrected constant-vector/reduction artifacts do not
remove all near-machine-precision splitting of mathematical ties. The same-input
rank algorithms agree. This numerical failure is not a regulatory observation.

A [separately specified post-hoc sensitivity protocol](../experiments/numerical-rank-sensitivity-v1.json)
reports every prediction grid: `1e-12`, `1e-10`, `1e-8`, `1e-6` and `0.01`.
All five reproduce 240/240 comparisons, preserving the 36 undefined statuses.
No grid was selected. Coarsening changes the rank estimand and **does not
retroactively pass frozen v1**; at the coarsest grid some correlations change
by about 0.109. See the [numerical review](rank-sensitivity-review.md).

### Primary leave-element-out result

The primary test comprises nine selected promoter elements, 10,365 qualifying
SNV rows and 3,485 positions. Errors are averaged over alleles within positions,
then positions within elements, then equally over elements. The endpoint is the
deposited signed log2 reporter-expression coefficient.

| Frozen model | Primary macro position MSE |
| --- | ---: |
| Training-mean baseline | 0.11639219 |
| Mutation chemistry | 0.11689921 |
| Chemistry + coordinate | 0.11652876 |
| Chemistry + delta-3mer | 0.11837090 |
| Chemistry + four motif deltas | 0.11682949 |
| Full fixed feature set | 0.11798887 |

**None of the five nonconstant models beats the mean baseline on this primary
metric.** This is a negative result for these fixed features, ridge fit and
selected cohort. It does not show that regulatory sequence is uninformative,
that the known TF mechanisms are false, or that a better context-aware model
cannot predict the endpoint.

This dataset-specific result is a controlled reanalysis rather than a novel
model method or a new biological mechanism. The original paper already reported
poor and heterogeneous computational predictability, under different models and
metrics; its reported correlations must not be numerically equated to our MSE.

### Why the tiny motif-versus-chemistry contrast is not a discovery

Adding four motif deltas reduces MSE relative to chemistry alone by
`0.0000697228`, approximately **0.060%** of the chemistry MSE.
The motif model still performs worse than the mean baseline.

The prespecified descriptive element-bootstrap interval for that improvement
is approximately `[-0.00058044, 0.00064285]` and crosses zero.
The chemistry-stratified rotation tail is `0.01` raw, `0.03` after within-family
Holm adjustment and `0.06` after adjustment across all six planned tails.

These numbers describe different conditional controls. The artificial rotation
distribution is not a calibrated biological null, the bootstrap conditions on
the fixed OOF fits and nine selected elements, and no positive result becomes
a mechanism merely because one adjusted tail is small. With 99 rotations,
the smallest raw tail is 0.01 and the smallest six-tail Bonferroni value is 0.06.
Changing the number of rotations after observing this threshold would require
a separately declared analysis rather than a claim that v1 established novelty.

### The row-random appearance does not survive the guarded HBB test

| HBB evaluation | Mean MSE | Full-model MSE |
| --- | ---: | ---: |
| Explicitly leaky row-random diagnostic | 0.02146213 | 0.02043621 |
| Whole-position, 20-base-guarded CV | 0.02193389 | 0.02960965 |

The small row-random advantage reverses under the prespecified harder
whole-position/guarded evaluation. That blocks a transportable-mechanism claim
from the easier split. The gap must not be attributed solely to leakage:
the split also changes training proximity, available data and test difficulty,
while source-library dependence can remain after guarding.

A secondary HBB-to-external diagnostic gives a modest coordinate-baseline
advantage (MSE `0.11557964` versus its mean comparator's `0.11672352`).
This is retained as a secondary observation; it does not reverse the negative
primary result or establish a motif mechanism. No model was selected from it.

## 6. Context differences: finite observed-label constraints

The independent-Python report includes exact-identity comparisons between
validation-only records. Two examples are:

| Same-variant source comparison | Matched SNVs | Pair MSE | Pearson | Identical-input two-context MSE bound |
| --- | ---: | ---: | ---: | ---: |
| LDLR independently constructed library repeat | 943 | 0.03393340 | 0.97265 | 0.00848335 |
| TERT GBM versus HEK293T cell context | 770 | 0.10288039 | 0.75903 | 0.02572010 |

For two observed coefficients `y1,y2` for the same variant, any predictor forced
to return the same value in both contexts has average squared error at least
`(y1-y2)^2/4`. Averaging over matched variants gives the displayed pair-MSE/4
quantity. This is standard algebra applied to the deposited estimates.

It does **not** establish a latent biological noise floor or identify the causal
source of the differences. LDLR and TERT have different targets, effect ranges,
contexts and estimation properties; comparing these two MSEs is not a controlled
decomposition of “technical noise” versus “cell biology.” The LDLR repeat is also
not a separate study. A context-aware predictor is not subject to the same
identical-output restriction.

Targeted re-reading of the original Results confirms that this is prior art,
not merely a general expectation. Kircher 2019 already reported an approximately
0.97 LDLR repeat correlation after its ten-tag selection, tested the same TERT
library in HEK293T and GBM, and discussed cell-dependent effects around an E2F
repressor site. Its TERT/GABPA analysis also explicitly examined allele-induced
PWM score changes. Our exact eligibility and metrics differ, so the published
correlation is not treated as an identical numerical replication.

The original study and earlier TERT work already establish the importance of
assay/cellular context and GABP-related regulation. The exact paired counts and
bounds are new-to-this-analysis quantifications; their global novelty is not
established and they should not be marketed as newly discovered biology.

## 7. What remains a candidate, and what does not

The strongest surviving contributions are reproducible **computational and
curation observations**:

1. The exact HBB source crosswalk and its support-stratified identity accounting.
2. The failure of these fixed simple sequence models to beat the primary mean
   baseline, despite an easier HBB split suggesting a small advantage.
3. The quantified mismatch between a context-blind single-value prediction
   contract and multiple context-dependent source estimates.

None currently identifies a new TF binding event, a new regulatory element,
a causal clinical mechanism or a previously unknown biological law.
The remaining exact-rank failure and separately labelled resolution sensitivity
are numerical-methods results, not genomic observations or evidence for a new
biological association.

A future biological claim would need a specific allele/context hypothesis,
new independent measurements, suitable controls and a new frozen evaluation.
Useful follow-ons include matched-context regulatory data, direct occupancy or
endogenous-expression evidence, and reporter transcript mapping when an RNA
mechanism is proposed. Expanding models or motifs after seeing v1 should be a
separate research round; the negative v1 result must remain retained.

## Reproduce the assessed source claim

```bash
python3 research/forensic-sources/reconcile_hbb.py --verify
python3 scripts/fetch_regulatory_motifs.py --validate-existing datasets/regulatory-motifs
```

Both were executed during this review. The forensic verifier checks the complete
derived source JSON; the motif validator checks all frozen count/source/hash and
scoring contracts. This review did not perform new model fitting or claim formal
Lean verification of the displayed real-valued squared-error identity.
