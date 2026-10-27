# Frozen-model prediction envelope

This outcome-free audit applies the five saved activity-model fits to every one
of 3,072 interior five-base SNV classes, for the shared model and three matched
cell-specific models. It retains 12,288 ensemble predictions, 61,440 fold
contributions, every tied extremum and all 64 zero-feature edits.

The independent calculation agrees exactly with the separately written transfer
predictor on 86,016 comparisons. All 101,376 local-versus-full integer count
checks pass. The full source artifact manifest is pinned to
`2ea476bc5413dc13486dec9f3d7e9162a50a0e47f1d83ce29efabc5f5cc5af9d`.
The original eight-file computation manifest and all nine original audit files
remain unchanged; this README is explanatory material outside that manifest.

The maximum absolute five-fit prediction is 0.0915981 for the shared model,
0.1029483 for HepG2, 0.1323450 for K562 and 0.1288580 for WTC11. These are
**interior, fixed-200nt model prediction bounds**, not bounds on biological
effects or calibrated cross-study measurements. The separate
[numeric review](../numeric-review-v1/README.md) includes boundary substitutions,
which can exceed these maxima.

The implementation contains its own exact feature counts and coefficient
calculation. `comparison-implementation.py` preserves the separately written
predictor used for comparison. `saved-fits.json` binds the extracted old fits;
no response fitting or new experimental outcome reading occurred. The script
supports a fresh `--out` directory; inspect `--help` for the required source
artifact and comparison-script arguments. Replay into a new directory
reproduced all nine original files byte-for-byte.

See [summary.json](summary.json), [protocol.json](protocol.json) and
[all-snv-predictions.csv](all-snv-predictions.csv) for the complete results.
