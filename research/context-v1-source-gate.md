# Context prediction v1: blocked at the semantic source gate

**Status: blocked before any model fit.** The immutable authority remains
[context-prediction-v1.json](../experiments/context-prediction-v1.json), SHA-256
`a089a6668ad73e943797df405b3ba510b1377190b416da1bd8df12bc368b77c2`.
Neither that JSON nor its human protocol was changed after freezing.

The two permitted Zenodo 10558183 workbooks were acquired with verified source
checksums. Header-only inspection found:

| Asset | Observed schema | Gate result |
| --- | --- | --- |
| Table S7 - joint lib design.xlsx | Element IDs, categories, genomic coordinates/strand and 230-nt adapter/element sequences | Design information; no observed raw-activity endpoint |
| Table S9 - folds and performance.xlsx | Fold, Sequence ID and three Element Specificity Score columns | Transformed specificity; explicitly forbidden as a v1 raw-activity label |

The required arithmetic mean of three measured, median-normalized log2 RNA/DNA
replicate activities could not be established from these sources. Published ESS
or CNN predictions cannot fill that gap. No v1 fit or claimed activity-prediction
result is produced. This is the protocol's intended fail-closed behavior, not a
negative predictive result.

Workbook hashes:

| Asset | SHA-256 | Published MD5 |
| --- | --- | --- |
| S7 | `f9ec4c6c2d230b0ef507a39e5d5c73b3b36f90ccc5758c33b442c7561e20e7e5` | `5d9b71452e9648e8473edd74a4ae26fc` |
| S9 | `0f9ff22d4b46fab36d419399e6ab743dc08a0813f812dc7bd96f90eb7e0f2b8a` | `ef967c8388d84da9ec81d22e11c868bd` |

A separately versioned `context-specificity-v2` may target the actual published
ESS columns after freezing its own semantic contract. It must disclose that
the files were already acquired under v1, their schemas were inspected and the
endpoint change follows a source-availability failure. It must not claim an
untouched payload, raw activity, prospective normalization or a successful v1.
The protocol author did not inspect numeric ESS values while preparing v2.
