# Exhaustive local pair interactions in the fixed model

For an interior pair of genuine single-base substitutions, let the model's
feature interaction be `F(double) + F(reference) - F(first) - F(second)`.
With reverse-complement-canonical trimer counts and GC count, every one of
2,304 oriented adjacent-edit classes and every one of 576 classes at distance
two has a **nonzero feature interaction**. At distances of three or more the
existing Lean proof establishes a zero interaction for any fixed trimer and
base weights; the canonical feature vectors are a special case of those sums. Thus, for
interior substitutions, this specific representation's feature interaction is
zero exactly when separation is at least three bases. The nonzero result is a
finite exhaustive calculation; the separate Lean file proves the zero direction,
not the exhaustive nonzero direction.

| Separation | Minimal shared-window reference contexts | Edit classes | Distinct interaction vectors | Count-vector L1 distribution |
| --- | ---: | ---: | ---: | --- |
| 1 nt | 256 four-base words | 2,304 | 502 | 96 at 4; 768 at 6; 1,440 at 8 |
| 2 nt | 64 three-base words | 576 | 144 | 576 at 4 |

All classes assume both replacements differ from their original bases. A class
includes both directed substitutions for one oriented reference. The minimal
context is the union of trimer windows containing **both** edited sites: two
windows for distance one and one for distance two. Windows containing zero or
one of the sites cancel in the four-way difference; GC cancels site by site.
Reverse-complementing a class or reversing both substitutions preserves its
interaction vector; reversing only one substitution changes its sign. The code
checks those symmetries for every class. It also independently subtracts four
full feature vectors over every combination of two external bases on each side
(737,280 comparisons) and one 200-base sequence per class (2,880 comparisons).
These counts agree with the previous full-context enumeration: 16 contexts per
minimal adjacent class and 256 per distance-two class. They are not estimates of
how often pairs occur in the genome or in a population.

The saved five-fit ensemble has these exhaustive absolute interaction maxima in
its predicted activity units:

| Separation | Shared model | HepG2 | K562 | WTC11 |
| --- | ---: | ---: | ---: | ---: |
| 1 nt | 0.0857051 | 0.1553290 | 0.1402557 | 0.1088221 |
| 2 nt | 0.0663751 | 0.0917066 | 0.0843153 | 0.0732702 |

For these maxima, the calculation treats each saved binary64 raw slope as its
exact rational value, averages the five fitted slopes, and divides trimer counts
by 198 for a 200-base insert. Exact rational numerators, denominators, and all
tied witnesses are in `summary.json`; per-class displayed predictions are
rounded binary64 conversions. The feature property does **not** require any
particular coefficients; a nonzero feature vector can still produce a zero
prediction for suitable weights. The maxima are properties of the fixed fitted
models and are not upper bounds on observed or true biological interactions.
The proof and exhaustive calculations do not establish endogenous function,
causality, or clinical significance.

`all-pairs.csv` contains every class and the 33-element signed integer feature
interaction vector, in the canonical trimer order listed in `summary.json`,
followed by GC. No reporter outcomes were read for this calculation; the saved
fits were previously trained on reporter activity. Reproduce the byte-identical
output from the repository root using a fresh destination:

```bash
python3 research/variant-locality/pair-geometry-v1/implementation.py \
  --out /tmp/helix-pair-geometry-replay
```

The output `manifest.json` hashes the generator, CSV and summary; this README is
explanatory material outside that computation manifest. The source model's
SHA-256 is recorded in the summary. Nothing in this calculation validates the
model's interaction predictions against four-construct measurements.

An independent [review](../pair-geometry-review-v1/README.md) recomputes all
four feature vectors on a different full 200nt context for every class, and
checks all predictions and exact-rational extrema against the saved fits.
