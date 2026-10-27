# Independent v2 verification

The independent implementation reproduces the fixed five-model benchmark on
55,276 elements across HepG2, K562 and WTC11, with 2,761 sequence/genomic
components and five guarded family folds. Exact source, feature, fold and
resampling checks pass. Numerical agreement does not resolve a discrepancy
between the supplied values and the cited ESS formula.

The endpoint is unchanged **published element specificity score**, not raw
activity. The original raw-activity protocol was blocked because its required
endpoint was unavailable. V2 was frozen separately before ESS numeric inspection
or model fitting.

| Model | Equal-element MSE | Equal-element MAE |
| --- | ---: | ---: |
| Global mean | 0.999925882 | 0.714638588 |
| Context mean | 0.999969298 | 0.714648089 |
| Shared DNA | 0.999646789 | 0.714664352 |
| Additive context/shared DNA | 0.999690205 | 0.714673354 |
| Independent context DNA | 0.865389189 | 0.666989056 |

The primary independent-versus-additive comparison reduces MSE by 0.134301016,
or 13.4343% relative to the additive model. The paired component-bootstrap range
is [0.129896207, 0.138927404]. MAE improves by 0.047684298, with range
[0.045739191, 0.049716105]. These are conditional fixed-prediction ranges for a
selected panel, not calibrated biological significance, model-refit uncertainty
or independent-study evidence. Equal-component weighting gives MSE improvement
0.125078893 as a fixed sensitivity.

The other comparisons remain diagnostics. Additive versus shared DNA worsens
MSE by 0.000043416. Additive versus context means improves MSE by 0.000279094
while slightly worsening MAE. They were not used to select a model, endpoint,
feature, context or reporting subset.

## Source and numerical checks

S7 supplies 60,000 design IDs and S9 supplies 55,338 actual score IDs. Worksheet
maximum row is not a record count. Every score ID exactly matches a design ID;
no suffix transform is used. Among joined rows, 52 fail core/reference agreement
and 10 have non200 spans. A blank score-ID row and 4,662 S7-only IDs are retained
in accounting. Every admitted score triple is finite.

The independent ZIP/XML parser reconstructs eligibility from every source row,
including exclusions. Each admitted core200 matches pinned hg38 at the unshifted
declared interval and strand, including eligible alternate contigs. Normalized
observations equal independently parsed raw cells exactly.

Different Hamming and guard implementations match all 372,850 candidates, every
component/member hash and every train/test/guard partition. No family is split
to improve balance. All 1,768,832 integer trimer counts and 1,824,108 feature
fractions match exactly.

All 829,140 predictions pass the fixed 1e-7 tolerance; maximum difference is
8.746e-12. All fitted statistics, sequence/group/context/pair losses, five model
bootstrap tables, three contrast tables and percentile ranges pass. The complete
2,000 by 2,761 draw table matches exactly. Rank endpoints are absent.

## Source-formula discrepancy

The cited final-paper formula describes per-element centering across three
contexts, predicting row sums near zero. Exact rational arithmetic on original
XML literals contradicts that prediction for the frozen preprint-associated
release:

- All 55,276 retained rows exceed absolute sum 1e-8; 55,275 exceed 1e-6.
- Mean absolute sum is 0.097278087, RMS sum 0.139460983, and maximum sum
  1.0987267147462894.
- The maximum is source ID `K562_peak46638`, physical S9 row 28,170, with
  HepG2/K562/WTC11 literals `-2.1243548522111699`,
  `7.8814064486789697` and `-4.6583248817215104`.

Raw headers and cell mapping are independently verified. The discrepancy exceeds
floating rounding. No recentering, snapping, row selection or target change is
performed. A source-version or normalization difference may explain it, but this
audit does not establish another formula. Fixed-zero MSE is approximately
0.999925820.

The source already predicts specificity, so this benchmark does not discover
context dependence. Its panel was selected using prior activity. Published
whole-source-panel normalization includes future held-out elements; grouped
fitting cannot undo it. This supports a model-relative predictive relationship
for three known cells on these frozen numbers. It does not establish raw-activity
accuracy, unseen-cell/genome-wide generalization, causality or novel biology.

## Attempts

The first and only independent Python fit succeeded unchanged. The initial Rust
output agreed numerically but its separate ledger audit failed because struct
serialization and sorted JSON produced different commitments. Root retained
that core and fixed the integrity path.

The corrected core added a complete feature dump and explicit unresolved
source-semantics status; models and targets did not change. Both comparisons
are retained. The separate exact-rational endpoint audit involves no fitting.
Paths, hashes, commands and statuses appear in `run-summary-v2.json`.
