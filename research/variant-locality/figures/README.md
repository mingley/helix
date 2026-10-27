# Feature-equivalence distribution

![Exact feature-equivalence histogram](panel-feature-equivalence.svg)

All 55,276 fixed 200nt panel records are included. The histogram counts exact feature-equivalent interior single-base substitution opportunities per record, out of 588 opportunities (196 positions times three alternate bases). There are 745,052 such opportunities among 32,502,288, or 2.2923%. Feature equality is a representation property and does not establish biological neutrality. No measured responses or fitted model predictions enter the plot.

`feature-equivalent-counts.csv` retains every integer bin from zero to the observed maximum 33, including small tail counts. SVG and 180dpi PNG are standalone exports. `manifest.json` binds both figures, the CSV, source scan summary/per-element data/manifest, plotting script and pinned research dependency lock. This explanatory README is outside that figure manifest.

From the repository root, use a new output directory:

```bash
source scripts/research-env.sh
python3 scripts/plot_variant_locality.py --out /tmp/helix-locality-figure-review
```

The plotting script verifies every committed scan artifact before reading the full per-element distribution, and refuses an existing output directory. A fresh replay reproduced SVG, PNG, CSV and manifest byte-for-byte with the pinned environment. The PNG was visually inspected for complete labels, caption, tails and biological-interpretation limits. The immutable scan and exact local-space audit were not modified.
