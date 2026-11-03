# GC and trimer features: an exact boundary identity

For every A/C/G/T sequence of length `L >= 4`, let `G` be the number of G or C
bases and let `C_j` count the overlapping three-base words in reverse-complement
class `j`. Let `g_j` be the number of G/C bases in a word of that class; reverse
complementation preserves this number. Define

```text
T = sum_j g_j C_j
B = 2 GC(first) + GC(second) + GC(penultimate) + 2 GC(last)
3 G = T + B.
```

Here `GC(base)` is 1 for G or C and 0 otherwise. The identity follows by counting
each base once for every three-base window containing it: an endpoint occurs in
one window, its neighbor in two, and every other position in three. Subtracting
these multiplicities from three gives exactly `B`. It uses integer counts and
does not depend on a fitted model, outcomes or floating-point approximations.

For the actual feature denominators, `f_j = C_j/(L-2)` and `f_GC = G/L`, so

```text
f_GC = (L-2)/(3L) * sum_j g_j f_j + B/(3L).
```

Thus, among equal-length sequences with the same `B`, the GC coordinate is
determined by the 32 canonical trimer frequencies. In particular, edits at
zero-based positions `2` through `L-3` leave `B` unchanged. Two such constructs
with identical trimer vectors must also have identical GC fractions, for any
deterministic predictor using these features. For a linear model with original
unstandardized trimer coefficients `beta_j` and GC coefficient `beta_GC`, the
same exact score can be written using trimer coefficients
`beta_j + beta_GC*(L-2)*g_j/(3L)` and the single boundary term
`beta_GC*B/(3L)`. A fitted standardized model first needs its centering and
scaling absorbed into unstandardized coefficients and intercept.

This is a reparameterization of a fixed score, not a claim that removing GC and
refitting gives the same model: ridge penalties and fitted coefficients depend
on representation. Sequences with different boundary statistics also need not
have GC determined by their trimer frequencies. For example, the two 200-base
sequences `AAC` repeated 66 times plus `AA`, and `ACA` repeated 66 times plus
`AC`, have exactly the same three canonical-trimer counts (66 each) but have
66 and 67 GC bases respectively. Their `B` values are 0 and 3. The GC coordinate
does distinguish this pair.

The separate [Lean proof](../../../proofs/GCTrimerIdentity.lean) establishes the
counting identity for arbitrary sequences of length at least four, preservation
under canonicalization, the equal-boundary implication, and an integer-score
reparameterization. It uses Lean's standard logical foundations and declares no
new axioms. The frequency formula above follows by dividing the count identity;
the proof does not establish the behavior of floating-point evaluation.

## Fixed-panel check

The [outcome-independent computation](analyze.py) checks the identity on all
55,276 existing 200-base panel sequences using canonical trimer counts. The
boundary statistic `B` takes values 0 through 6 with counts 3,702; 7,038;
11,202; 12,956; 10,602; 6,506; and 3,270, respectively. GC count and the
GC-weighted trimer count have Pearson correlation approximately **0.999576**
on these records. This is a descriptive property of the fixed panel, including
its duplicate or related sequences; it is not a population or biological
effect estimate. Exact sufficient statistics and the squared correlation as a
rational number are recorded in [summary.json](summary.json).
An [independent raw-window implementation](../gc-trimer-independent-v1/README.md)
reproduces the full summary and checks the rational-score transformation.

The panel is an existing explored source; its JSON parser instantiates the
previously used observation fields, but the script analyzes only DNA sequence.
It reads no new outcome data and makes no claim of genome-wide novelty or
endogenous mechanism. Reproduction from the repository root:

```bash
python3 -m unittest discover -s research/variant-locality/gc-trimer-identity-v1 -p 'test_*.py'
python3 research/variant-locality/gc-trimer-identity-v1/analyze.py --verify
source scripts/env.sh
lean proofs/GCTrimerIdentity.lean
```

For a new output, pass `--out` followed by a nonexistent directory instead of
`--verify`. Tests exhaust all sequences of lengths four through six, check
reverse-complement invariance and edge versus interior edits, exhibit the
same-trimer/different-GC pair, and compare both linear-score parameterizations
using exact rational arithmetic. The result concerns the exact mathematical
features; floating-point evaluation can introduce rounding differences.
