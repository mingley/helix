# Published specificity columns: source arithmetic and semantic limits

The frozen source columns are numerically well defined, but their relationship
to the paper's stated ESS formula remains **unresolved**. They do not obey the
expected unweighted zero-sum identity. They instead occupy a positive-weighted
two-dimensional plane and have individually standardized source columns. This
is consistent with an additional column-scaling operation; it does not identify
the historical processing pipeline.

The [v2 protocol](../experiments/context-specificity-v2.json) remains unchanged,
SHA-256 `8f0576845fecb5cfcf9bfa354fab446d849244448a3489905a6ce33bdeb21cec`.
Its zero-sum audit explicitly requires reporting semantic/audit failures and
prohibits projecting targets, snapping model terms, filtering by residual or
inventing a new endpoint. It does not authorize repairing the labels. The
previous [v1 raw-activity gate remains blocked](context-v1-source-gate.md).

## What the source establishes

The final paper, DOI [10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9),
Methods “Calculation of element specificity scores”, describes per-cell
z-score transformation followed by subtracting each element's mean across cell
types. Thus its stated formula is

`E_ic = Z_ic - mean_d(Z_id)`.

For three cells this entails `E_iH + E_iK + E_iW = 0`, irrespective of the
across-element SD convention. No subsequent per-column standardization is
described in that paragraph. The independently retrieved official
[XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11903340/fullTextXML)
has SHA-256
`41a554f31ea67adb1ce983b05ed9194e8fa1596bd4e6d2208df1854de9f33320`.

The linked [preprint, DOI 10.1101/2023.03.05.531189](https://doi.org/10.1101/2023.03.05.531189),
PMCID PMC10028905, states the same operation under its ESS calculation methods.
Its [XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10028905/fullTextXML)
has SHA-256
`bc7d21c09f401f3d95dee58d4f53abc84670590a36b86740b8d3899e85626f25`.
Its Supplementary Table 7 description promises activity/replicate summaries and
ESS, whereas the frozen Zenodo workbook called “Table S7 - joint lib design”
contains design information. This reinforces the need to distinguish releases
and table names; it does not supply the missing transformation for the frozen
S9 columns.

Code-only inspection found no ESS creation step in the inspected loaders:

| Pinned source | Inspected behavior |
| --- | --- |
| `visze/sequence_cnn_models`, snapshot `7638fce137db0445123efe8d7e2c35e248fafc5f`, `workflow/rules/regression_input.smk` | Joins externally supplied labels to FASTA sequences and partitions folds |
| `autosome-ru/human_legnet`, snapshot `a20dbc7d0b08daae3af1aebcfb8a71d9f198ab5d`, `datamodule.py` and `dataset.py` | Loads supplied targets as `mean_value`, partitions folds and returns them |
| Same human_legnet snapshot, `utils.py` | Seed, sequence, reverse-complement and encoding utilities |
| Same human_legnet snapshot, `analyze_concordance.R` | Allelic-binding variant comparison, not ESS construction |

These limited negative checks do not rule out upstream, external or unpublished
processing. They cannot certify a hypothetical post-centering rescale. The
relevant code URLs are the corresponding paths under
`https://raw.githubusercontent.com/<repository>/<snapshot>/`; source hashes for
the additional checked files are `feab3e55fbc8aa52b473cbf6be790affdc8d3318e221ca19a4732d3968a668a5`
for `datamodule.py`,
`88304f6235a39eb83dad805d30175e4d97ef8a695eabf032303d6536a1db0fae`
for `utils.py`, and
`bdb03695e0f0704f4b1b74ac9cb11aa942c459d1a79ffb849c9fb22531553228`
for `analyze_concordance.R`. No additional outcome table was acquired for this
provenance review.

## Separately frozen arithmetic diagnostic

The initial v2 fit exposed the row-sum discrepancy. Before inspecting source
moments or singular values, the parent froze the separate post-hoc
[arithmetic audit](../experiments/context-source-arithmetic-audit-v1.json),
SHA-256 `8f092c7ff075d7f455791e5a72f189781e611e3a6bda69b0fd9d8a7947990eae`.
This diagnostic examines the unchanged 55338 finite source rows and 55276
reference-valid admitted rows. It does not refit a predictive model or supply
fitted weights to the benchmark.

The [offline implementation](../scripts/audit_context_source_arithmetic.py)
verifies source/panel/protocol hashes, exact worksheet headers and IDs, and that
every admitted target equals its raw source-cell float value. Four synthetic
checks cover exact-decimal cancellation, a rescaled zero-sum plane, a full-rank
matrix and empty-input rejection. They passed before the real audit. The full
report is [summary.json](../runs/context-source-arithmetic-audit-v1/summary.json),
SHA-256 `acf227b68396006a2934aedae4f088f867f4cbcf937d4586668d20d0d589be3a`;
the archived implementation SHA-256 is
`9a9779f80ea3002eab0ba904530cf1ec5c1bc3d7bb2f96155018d0369a3dbf42`.

All-source moments show per-column sample SD one to numerical precision:

| Context | Mean | Population SD | Sample SD |
| --- | ---: | ---: | ---: |
| HepG2 | −8.01e−16 | 0.999990964576623 | 1.000000000000001 |
| K562 | −1.27e−15 | 0.999990964576629 | 1.000000000000006 |
| WTC11 | −6.79e−16 | 0.999990964576625 | 1.000000000000002 |

All 55338 source rows have absolute unweighted sums above `1e-8`; 55337 exceed
`1e-6`. The RMS sum is 0.1394458170 and maximum is 1.0987267147. The admitted
subset has RMS 0.1394609830 and the same maximum. This is a source-score
property, not an error introduced by reference filtering or predictive fitting.

The exact source witness is ID `K562_peak46638`, Excel row 28170:

| Context | Raw XML decimal value |
| --- | ---: |
| HepG2 | −2.1243548522111699 |
| K562 | 7.8814064486789697 |
| WTC11 | −4.6583248817215104 |
| **Exact decimal sum** | **1.0987267147462894** |

Raw and column-centered matrices both have numerical rank two under the fixed
relative cutoff `1e-10`:

| Cohort/matrix | Largest singular value | Second | Smallest | Rank |
| --- | ---: | ---: | ---: | ---: |
| All source, raw | 300.6652535 | 274.9752813 | 4.82e−13 | 2 |
| All source, centered | 300.6652535 | 274.9752813 | 6.09e−13 | 2 |
| Admitted, raw | 300.4520467 | 274.8531727 | 4.78e−13 | 2 |
| Admitted, centered | 300.4520433 | 274.8531656 | 6.00e−13 | 2 |

The all-source raw right-null weights normalized to mean one are:

| HepG2 | K562 | WTC11 |
| ---: | ---: | ---: |
| 1.024591105964359 | 0.907411107171192 | 1.067997786864450 |

Every weight is positive. The weighted raw-source residual has RMS
`3.55e-15`, maximum absolute `3.05e-14`; the other cohort/centering variants
give essentially the same weights and similarly tiny residuals. The full
report retains every moment, minimum/maximum, row-sum quantile, singular value,
signed weight and residual summary without choosing a favorable variant.

These observations are **compatible with** taking a centered specificity matrix
and subsequently standardizing each column separately. If `Y w ≈ 0` for
positive weights `w`, then `E = Y diag(w)` is approximately zero-sum. When the
Y columns already have unit sample SD, standardizing this constructed E returns
Y. This proves that such a representation exists, not that it was the producer's
history; other latent-score constructions can generate the same matrix. The
overall scale of a preceding E is also unidentified. Neither raw expression
nor historical normalization can be recovered by fitting the null vector.

## Reporting boundary

Keep the original target cells, fits and out-of-fold losses unchanged. Report
distinct statuses rather than collapsing them into one passing flag:

| Check | Defensible status |
| --- | --- |
| Exact frozen workbook/column binding | Verified by hashes, schema and cell comparison |
| Expected unweighted ESS formula | Failed for these source values |
| Historical extra normalization | Unresolved; weighted-plane explanation is compatible, not proven |
| Arithmetic diagnostic | Completed under its separately frozen post-hoc contract |
| Predictive numerical replication | Determined by the separate benchmark comparator, not this audit |

A result may be described narrowly as prediction of the three published S9
columns, with the primary A-versus-I comparison retained as specified. It must
carry the failed/unresolved semantic flag and cannot be presented as complete
v2 semantic validation or final-paper centered-ESS prediction. This is neither
a repaired v1 raw-activity study nor an independent biological validation.

Balanced-model identities remain valid even when observed targets are not
zero-sum: ridge linearity still gives mean I slopes equal to S slopes, so the
predeclared pair-contrast/MSE identity remains an implementation check. Those
checks do not certify the source endpoint. Three correlated cell scores are
not three independent experimental units; component resampling keeps their
dependence together.

No inferred source weights may enter features, targets, folds, model fitting or
reported v2 primary losses. Re-centering or rescaling would define a new
endpoint requiring a separately versioned exploratory protocol. The unresolved
provenance limits scientific interpretation even if prediction reproduction
passes. It is a documented methods/source finding, not genomic novelty.
