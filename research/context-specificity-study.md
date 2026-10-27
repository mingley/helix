# Predicting published context scores with an unresolved source transform

## Result and semantic status

**The numerical benchmark improves prediction of the supplied source columns,
but their biological/processing meaning remains unresolved.** Every admitted
element fails the unweighted zero-sum implication of the published ESS formula
by more than 1e-8. The maximum absolute row sum is **1.0987267**, far beyond
rounding. This is a source-semantic audit failure, not a successfully validated
implementation of the final paper's stated element specificity score (ESS).

On the unchanged published columns, the context-specific DNA model
has MSE **0.8653892**, versus **0.9996902** for additive DNA/context prediction:
an improvement of **0.1343010**, about **13.434%**. The conditional component
resampling range is **[0.1298962, 0.1389274]**. These numbers describe prediction
of the frozen source columns. They do not establish a new biological mechanism,
raw-activity prediction, or a resolved source preprocessing history.

The complete corrected independent numerical comparison passes all **243
categories**, including exact integer features and **1,824,108 bitwise-equal
feature fractions** on the 200-base cores, all **829,140 model predictions**,
fits, memberships, losses and resampling checks. Numerical agreement cannot
repair endpoint semantics. The complete workflow exits **2** for the unresolved
scientific source-formula validation.

The preceding [context-prediction-v1](context-v1-source-gate.md) remains blocked
before fitting: its required observed activity means were absent. The distinct
[v2 protocol](../experiments/context-specificity-v2.json), frozen at
`2026-10-27T13:39:07.406725+00:00`, has SHA-256
`8f0576845fecb5cfcf9bfa354fab446d849244448a3489905a6ce33bdeb21cec`.
It is retained unchanged, including the formula now contradicted by source
arithmetic. V2 does not retrospectively pass v1 or validate the disputed formula.

## Source, exposure and admitted panel

The source publication is Agarwal et al., *Massively parallel characterization
of transcriptional regulatory elements*, Nature (2025), DOI
[10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9),
PMID 39814889, [PMC11903340](https://pmc.ncbi.nlm.nih.gov/articles/PMC11903340/).
Its joint library was tested in HepG2, K562 and WTC11 human iPS cells.
The paper already predicts cell specificity with multi-task sequence models
and interprets cell-specific motifs. The general idea is prior art.

The exact frozen source is the preprint-associated Zenodo 10558183 S7/S9
release. The workbooks were acquired under v1; schemas, IDs, sequences and
coordinates were inspected before v2. Numeric score analysis/fitting followed
the separate v2 freeze. The source-formula discrepancy was then observed and a
[separate post-hoc arithmetic audit](../experiments/context-source-arithmetic-audit-v1.json)
was frozen before its diagnostic execution. No untouched-payload, wholly blind
history or outcome-independent selection of that arithmetic question is claimed.

Final-paper Methods describe an endpoint class; they do not prove that old S9
is byte-identical to final Supplementary Table 10 or underwent the same
processing. The preprint's ESS paragraph gives the same zero-sum construction,
so a simple final-versus-preprint wording difference does not resolve the
discrepancy. The preprint also names its ESS table Supplementary Table 7,
whereas the frozen Zenodo S7 is design-only: version/table mapping must not be
assumed from numbering alone.

| Item | Observed accounting |
| --- | --- |
| S7 design records | 60,000 oligos, each 230 bases: 15-base adaptor + 200-base core + 15-base adaptor |
| S9 finite, nonempty unique score IDs | 55,338; exact S7 ID matches |
| Earlier 59,107 count | Superseded worksheet maximum-row accounting, not physical record count |
| Reference-valid admitted elements | 55,276, with all three source columns |
| Admitted measurements | 165,828; correlated contexts, not independent biological replicates |
| Root groups | 2,761; largest group 129 elements |
| Guarded outer folds | Five; minimum retained train/test genomic gaps exceed 10,000 bases |

The source core is the feature/reference input; source IDs are joined exactly.
Fixed invalid/missing-coordinate and reference exclusions are retained in the
import ledger, without coordinate rescue or score-dependent selection.
The panel was originally selected using earlier activity bins and variable
promoters. Grouped evaluation describes this selected, eligible source panel,
not arbitrary human DNA or the genome-wide prevalence of specificity.
See [the literature/design review](context-literature.md) for detailed scope.

## The stated formula and the supplied scores disagree

The final paper's Methods, **“Calculation of element specificity scores”**,
states:

> To compute element specificity scores (ESSs) using activity scores from the
> joint library, log2(RNA/DNA) values from each cell line were first
> z-score-transformed. Then an ESS for each element was computed by subtracting
> the element’s score in each cell type by the mean element score across cell
> types. A full table of ESSs is provided (Supplementary Table 10).

This is a CC BY4.0 excerpt. The official
[XML source](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11903340/fullTextXML)
had SHA-256
`41a554f31ea67adb1ce983b05ed9194e8fa1596bd4e6d2208df1854de9f33320`.
The [preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC10028905/), DOI
`10.1101/2023.03.05.531189`, gives the same sequence of operations in
“Calculation of element specificity scores (ESSs).”; its reviewed XML hash is
`bc7d21c09f401f3d95dee58d4f53abc84670590a36b86740b8d3899e85626f25`.

With source activity `a[i,c]` and source per-cell z-score `z[i,c]`, that
description implies:

```text
e[i,c] = z[i,c] - (z[i,HepG2] + z[i,K562] + z[i,WTC11]) / 3
sum_c(e[i,c]) = 0
```

The mean is across all three cells, not the other two. SD convention and exact
source normalization population are not specified here. Those omissions do
not explain a failure of row centering: subtracting a row mean gives zero
row sum regardless of the preceding column scales.

In the supplied S9 columns, **all 55,338 source rows**, including **all 55,276
admitted rows**, have absolute unweighted sums above 1e-8. The source RMS row
sum is **0.1394458**, median absolute sum **0.0665197**, and maximum absolute
sum **1.0987267**. An exact-decimal witness has sum
`1.0987267147462894`, excluding ordinary binary summation as the explanation.

V2 retains the values unchanged: no projection, rescaling, repair, residual-based
filtering or predictive refit. The formula/column relationship is **unresolved**.
The benchmark may describe the supplied scores, but cannot be represented as a
fully passed validation of the stated final-paper ESS endpoint.

## A post-hoc arithmetic observation narrows the possibilities

The separately frozen arithmetic diagnostic used all 55,338 source rows and
the 55,276-element admitted subset. It reports:

- Numerical rank two in each three-column matrix under its fixed SVD threshold.
- Whole-source column means near zero and sample SDs near one.
- A fitted positive null vector, normalized to mean one, of approximately
  **(1.024591106, 0.907411107, 1.067997787)** in HepG2/K562/WTC11 order.
- A maximum absolute weighted row sum of **3.055e-14** for the full source,
  versus **1.0987267** for the unweighted sum. The admitted subset has the same
  plane to numerical precision.

These facts are **consistent with** row-centered contrasts followed by unequal
column standardization. For example, if zero-sum contrasts `t[i,c]` were
rescaled as `y[i,c]=t[i,c]/s[c]`, then `sum_c(s[c]*y[i,c])=0`, while the
unweighted sum need not vanish. This is one compatible algebraic description;
it is not proof of the author's historical operation order, scale parameters,
code version or source table lineage.

The fitted weights are forensic diagnostics, not targets or model covariates.
They do not enter the frozen benchmark. Both paper versions give the
unweighted formula, and inspected author model-loading/prediction code has not
identified an authoritative creation history for S9. Access to the relevant
source preprocessing code and version-linked intermediate table would
discriminate the proposed scaling explanation from other equivalent histories.
The [semantic review](context-specificity-semantic-review.md) retains the
exact-decimal witness, all-source moments, pinned code checks and audit receipts.

The rank-two observation means the supplied columns do not furnish three
independent numerical dimensions. It does not establish independent biological
replicates, recover raw activity, or validate the published formula. A zero
value must not be read as inactivity. Under the intended contrast endpoint it
means a relative reference point, and the exact reference point in old S9
remains insufficiently traced.

## Fixed models and conditional algebra

The five models are global mean (G), context means (C), shared DNA (S),
additive shared DNA/context (A), and context-specific DNA (I). Their 33 features
are 32 reverse-complement-canonical overlapping trinucleotide frequencies and
GC fraction, computed only from the 200-base core. Fixed ridge alpha is ten;
feature means/population SDs use unique training elements only, with
unpenalized response means. No motif, score, source prediction, count, activity
bin or test-derived statistic is a fitted feature.

Under the **ideal stated zero-sum formula**, the shared per-element mean target
would be zero, S would collapse to G, and A to C. This is conditional endpoint
theory, **not an observed property of S9**. Its nonzero shared mean can contain
a mixture of relative contrasts after unequal column scaling; it must not be
called recovered universal activity. The observed A–C difference is therefore
not evidence for a common biological activity component.

Balanced linear ridge nevertheless implies that the mean context-specific
slope equals the shared slope, so I and A have the same across-cell prediction
mean. This identity does not require zero-sum source labels. The initial saved
model checks pass, including a per-element squared-loss contrast identity
residual at most **3.065e-14**.

For three errors with mean `r_bar`:

```text
sum_over_pairs (r_c-r_d)^2 = 3 * sum_c (r_c-r_bar)^2
pair_contrast_MSE = 3 * (cell_MSE - mean_error^2)
```

With I/A's shared prediction mean and identical element/group weights:

```text
pair_contrast_MSE(A) - pair_contrast_MSE(I)
    = 3 * (cell_MSE(A) - cell_MSE(I))
```

This is an implementation/weighting diagnostic, not a second independent
biological validation or an additional discovery.

## Predictive result on the supplied columns

| Model | Equal-element, three-cell MSE |
| --- | ---: |
| Fixed zero audit | 0.999925820 |
| G: global mean | 0.999925882 |
| C: context means | 0.999969298 |
| S: shared DNA | 0.999646789 |
| A: additive DNA/context | 0.999690205 |
| I: context-specific DNA | **0.865389189** |

The declared primary A-minus-I improvement is **0.134301016**, approximately
**13.434% of A's MSE**. I's MAE is **0.6669891**, versus A's **0.7146734**.
All models and the fixed-zero audit are retained. Equal-component MSE, a fixed
sensitivity rather than the primary weighting, is 1.0279388 for A and 0.9028599
for I.

The 2,000 paired component resamples of fixed out-of-fold losses give a
conditional 95% percentile range **[0.1298962, 0.1389274]** for primary
improvement. This ignores target-normalization uncertainty, refitting,
source-study variation and possible residual component dependence. It is not a
calibrated significance result or a population confidence guarantee.

Corrected numerical reproduction reports **243/243 checks passing**, including
all five models' 165,828 predictions each, with maximum prediction difference
**8.746e-12** under the frozen 1e-7 tolerance. All 55,276 source integer feature
vectors and 1,824,108 feature fractions match exactly. Fits, full grouped/guarded
memberships, pair diagnostics, losses and common resampling tables reproduce.
No frozen tolerance was widened.

The initial comparison had 234 checks and lacked the final exact raw-feature
comparison. A separate initial artifact serialization-order integrity failure
was corrected without changing the targets or model objectives. Earlier
attempts are retained; the corrected independent comparison establishes the
numerical pass, while the scientific source-formula failure remains unresolved.

## Remaining evidence and retained record

The final [study snapshot](results/context-specificity-v2/README.md) and
[figures](results/context-specificity-v2/figures/) will preserve the frozen
protocols, source exclusions, grouping/guard receipts, every model/contrast,
failed attempts, independent comparisons and source-arithmetic audit.

| Check | Current status |
| --- | --- |
| V1 raw-activity source gate | Blocked, unchanged; no fit |
| V2 source headers/values and reference eligibility | Bound to frozen source; 55,276 admitted |
| V2 stated ESS formula versus supplied columns | Failed/unresolved; not rounding |
| Post-hoc weighted-plane arithmetic | Completed; no historical normalization proof |
| Fits, resampling and numerical comparison | All 243 categories pass |
| Independent exact-core raw features | All counts and 1,824,108 fractions exact |
| Authoritative S9 preprocessing/version provenance | Unresolved |
| New biological mechanism or global novelty | Not established |

The bounded result is predictable variation in fixed published columns for
held-out groups in three known experimental contexts. The assay does not
isolate cell identity from delivery/integration/culture, restore native
regulatory contacts, or identify causal TFs. TF/motif claims require directed
sequence changes, occupancy, perturbation/rescue and endogenous validation.

The source already studies sequence-based specificity prediction. This
benchmark's observed improvement, source discrepancy and reproducibility
controls can be useful new-to-analysis findings. They do not establish new
global methodology or biology. The semantic failure remains visible even if
the numerical benchmark is reproducible.

## Reproduction and next scientific goal

After the pinned setup described in the [README](../README.md), run with a new
output prefix:

```bash
bash scripts/check_context_study.sh runs/context-review
```

This rederives reference checks/groups, runs both implementations and produces
all diagnostics. It obtains the checksum-verified 835 MB hg38 reference for
full rederivation; `bash scripts/setup_reference.sh` can cache it separately.
Exit **2** retains the failed/unresolved scientific semantic gate even when
numerical reproduction passes. The published snapshot supports review without
that reference download or refitting.

The next goal is authoritative, version-linked normalization/raw activity or
independently compatible measurements, rather than further fitting to force
the old formula. No projection or rescaling of v2 targets is authorized as a
repair. The earlier negative [promoter study](regulatory-study.md), its 26/240
raw exact-rank failures and the blocked context v1 remain unchanged.
