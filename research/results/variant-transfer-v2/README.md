# Exploratory external variant-effect transfer

The fixed five-fit model does not beat zero-effect prediction on 2,153 complete
variants in 172 components passing the frozen prior-source exclusion guards. Macro MSE is 0.2033556 versus 0.2033201;
zero-minus-model gain is -0.0000354532, conditional 95% range
[-0.00144180, 0.00140062]. Useful-transfer rule: **false**.

This is explicitly exploratory after the preserved v1 source-gate failure.
Exactly 20 literal NA cells in 19 variants are operationally unreported, with no
claim about why they are missing. No new-outcome fitting or calibration occurred.
The metadata cohort, source/model/context/grouping/metrics and resampling remain
fixed. V1 is retained separately; its original source gate remains failed.

All 17 independent evaluation comparison categories pass, including 852,588 exact
feature numerators, all 51,672 fold/ensemble predictions and 344,000 bootstrap
indices. Of 4,344 selected source fields, 4,324 are bit-exact numeric values and 20 are
the literal NA/missing records. Full global grouping separately passes 4,689 checks.
Raw assembly/source-class preparation is an authenticated handoff; the numerical
comparison does not claim full independent reconstruction of every source step.

Every original source ledger, missingness record, prediction, input, loss,
bootstrap draw, implementation, independent comparison and standalone figure is
retained in its named directory. `publication.json` hashes every byte except
itself and distinguishes numerical agreement from predictive failure.

See the [study report](../../variant-transfer-v2-study.md) and
[original failed-v1 report](../../variant-transfer-study.md).

Verify with standard-library Python, without source download or refitting:

```bash
python3 -S research/results/variant-transfer-v2/packaging-implementation.py \
  --verify research/results/variant-transfer-v2
```
