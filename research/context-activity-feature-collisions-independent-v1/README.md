# Independent check of the v3 activity feature collision

This is a **post-hoc check of already-explored labels** from the same reporter
experiment, not an independent biological validation. It reconstructs the
[primary audit](../context-activity-feature-collisions-v1/README.md) without
importing its code. It counts raw DNA windows into the 32 canonical trimer
classes and GC count, checks all 55,276 feature vectors against the published
v3 feature snapshot, independently derives group and test-fold membership, and
compares every primary assignment and every field of its collision and summary
reports. It computes within-class squared-error floors as exact rational
arithmetic on the retained JSON-decoded binary64 values.

The one full-vector collision contains `ENSG00000144566` and
`ENSG00000163576`. Both occupy hg38 `chr3:19946928-19947128`, share their
component and test fold 4, and have distinct 200nt cores that are reverse
complements. For each cell context, all five saved model predictions are
exactly equal between the two records. The three previously published activity
differences are approximately `0.548`, `0.328`, and `0.274` log2 units in
HepG2, K562, and WTC11. The source-selected collision-only pooled squared-error
floor is approximately `0.040247`, and its equal-element, equal-cell
contribution to the full 55,276-element panel is approximately
`0.00000145622`. These are finite observed-label arithmetic bounds, not a
biological noise floor or a claim of generalization to unseen experiments.

## Original source checks and precision

The verifier directly reads the pinned original S7 and S11 workbooks using a
separate ZIP/XML parser. S7 rows 58011 and 58298 contain the respective `+`
and `-` designs. Their 230nt source oligos have identical 15nt left adapters
(`AGGACCGGATCAACT`) and identical 15nt right adapters
(`CATTGCGTGAACCGA`), with reverse-complemented 200nt cores. The complete 230nt
oligos are not themselves global reverse complements. These design facts do not
determine why the previously observed activity values differ; assay and
measurement effects remain possible.

S11 rows 55546 and 55266 contain the selected activity aggregates. Some
original XML numeric lexemes are longer than the JSON spellings: for example,
`1.8380000000000001` and `1.838` decode to the same binary64 value. The
verifier compares decoded binary64 values, then applies
`Fraction.from_float` for exact floors. Thus the rational results use the
retained binary64 values exactly; they do not treat either decimal spelling as
an exact biological measurement. Numerators and denominators in the primary
artifacts are decimal strings so they remain exact in generic JSON consumers.

The S11 streaming XML parser instantiates other row contents and the workbook's
shared-string table while interpreting only the six selected aggregate numeric
cells. It reads IDs while locating the selected rows. The verifier also loads
all previously explored v3 activity labels and the panel JSON; the latter
instantiates old specificity fields that are not used as activity labels. No
new assay source is acquired or analyzed.

## Reproduction

The [receipt](receipt.json) binds the source workbooks, v3 panel, activity,
features, grouping, saved study, primary artifacts and independent code by
SHA-256. From the repository root, with the retained sources present:

```bash
python3 -m unittest discover -s research/context-activity-feature-collisions-independent-v1 -p 'test_*.py'
collision_review_dir=$(mktemp -d)
python3 research/context-activity-feature-collisions-independent-v1/verify.py \
  --receipt "$collision_review_dir/receipt.json"
cmp research/context-activity-feature-collisions-independent-v1/receipt.json "$collision_review_dir/receipt.json"
```

The receipt path must be new and is created exclusively. Synthetic controls
cover reverse-complement and non-reverse-complement feature collisions, exact
binary64 floor arithmetic, XML lexemes and formula rejection, altered
assignment/summary comparisons, and output-path safety under optimization.
Reusing the same source data is computational verification, not another
experiment or a causal test of orientation.
