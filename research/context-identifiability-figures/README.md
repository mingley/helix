# Synthetic identifiability illustrations

These are original, fixed **synthetic mathematical examples**. They contain no
S9 values, reconstructed experimental activity, inferred source transformation
or newly discovered biological mechanism.

## Common shift at the centering input

[SVG](identifiability.svg), [PNG](identifiability.png),
[exact CSV](identifiability-exact.csv), [frozen inputs](inputs.json),
[provenance manifest](manifest.json).

`(2,4,3)` and `(7,9,8)` differ by a common shift of five and both center to
`(-1,1,0)`. This illustrates the input-stage ideal map `P(x)=x-mean(x)*1`.
It does not recompute upstream column standardization or establish that this
map produced the source workbook.

## Identical normalization moments do not identify the input

[SVG](normalization-controlled.svg), [PNG](normalization-controlled.png),
[exact CSV](normalization-controlled-exact.csv),
[frozen inputs](normalized-panel-inputs.json),
[provenance manifest](normalization-controlled-manifest.json).

Two different four-element, three-context panels each have column sums zero
and sums of squares fourteen. Their sample SD is the same `sqrt(14/3)`.
Applying the same column standardization and across-context centering gives
identical outputs, despite different common row components. Exact CSV fields
retain rational inputs/moments and symbolic normalized values such as
`1/sqrt(14/3)`, rather than presenting a rounded float as exact.

Both constructions are covered by the ideal-map
[Lean development](../../proofs/ContextIdentifiability.lean) and
[technical note](../context-identifiability.md). The formal integer operator
`center3` is three times ideal rational centering. A source hash in the figure
manifest is provenance, not itself a Lean-check receipt or empirical evidence.
The [literature review](../context-endpoint-literature.md) distinguishes this
established mathematics from case-specific source findings.

## Regenerate

Use the pinned research environment and a new output directory:

```bash
source scripts/research-env.sh
python3 scripts/plot_context_identifiability.py --out runs/identifiability-review
python3 scripts/plot_context_identifiability.py \
  --example normalization-controlled --out runs/identifiability-review
```

The script verifies each frozen synthetic-input hash, computes exact checks
with rational arithmetic, and refuses to overwrite existing outputs. It records
its implementation, formal-source and dependency-lock hashes, runtime versions,
artifact hashes and interpretation limits. SVG dates are omitted and fixed
hash salts avoid random IDs. Rendering is deterministic with the pinned
Matplotlib/font/runtime; alternative rendering stacks may differ.
