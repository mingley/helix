# Published cell-specificity prediction

The target is source ESS, not raw log2 reporter activity. Source normalization used the full selected panel, including future held-out elements. Model fits alone are grouped and training-only. This is not a prospective target-preprocessing test, an unseen-cell test or a new biological mechanism.

| Model | Equal-element MSE | Equal-element MAE |
| --- | ---: | ---: |
| GlobalMean | 0.999925882 | 0.714638588 |
| ContextMean | 0.999969298 | 0.714648089 |
| SharedDna | 0.999646789 | 0.714664352 |
| AdditiveContextDna | 0.999690205 | 0.714673354 |
| IndependentContextDna | 0.865389189 | 0.666989056 |

The maximum absolute sum of three published scores is 1.098726714746. This source check must be reconciled with the stated centering formula before claiming full semantic validation; original scores are retained unchanged. Under ideal unweighted centering, shared activity is removed, but that formula must not be assumed for this workbook. Mean/shared baselines do not establish genomic inactivity. Primary comparison is independent-context slopes versus additive context/shared-DNA slopes. All source, grouping, model and bootstrap outputs are retained. Independent numerical comparison is required separately; an integrity audit does not establish empirical truth or reproduction.
