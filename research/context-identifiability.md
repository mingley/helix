# What centered three-cell specificity can identify

Ideal across-cell centering loses one scalar per element. It can preserve
differences between cells while discarding their shared level. This is a
mathematical identifiability result, not evidence that the frozen source
workbook implemented a particular normalization formula. Its unequal weighted
constraint and failed unweighted row-sum check leave that historical formula
unresolved. None of the results below corrects, projects, refits, or reinterprets
the published numerical targets as raw transcription.

## The observable equivalence class

Let `x_i=(x_i1,x_i2,x_i3)` be the three coordinates for element `i` **at the input
of a centering operation**, and let

```text
P = I - 11ᵀ/3
specificity_i = P x_i = x_i - mean_cell(x_i) · 1.
```

For any scalar `h_i`,

```text
P(x_i + h_i·1) = P x_i.
```

The shifts can vary arbitrarily across elements. For example, `(2,4,3)` and
`(7,9,8)` both center to `(-1,1,0)`, although their shared level differs by five.
The units in this example are arbitrary input-stage units, not measured source
activities. Every constant row `(h,h,h)` centers to zero; zero specificity does
not establish inactivity. The same observation has arbitrarily distant input
witnesses if no additional bounds or upstream constraints are imposed.

The converse is also exact: two rows have equal centered observations **if and
only if** they differ by a common shift. Thus the observable object is an
equivalence class `{x_i+h·1}`, rather than a unique row of activities. No
function taking only the centered row can reconstruct every original row.

Two contrasts suffice to specify the centered observation. Set
`u=x_1-x_3` and `v=x_2-x_3`. Then

```text
3 P x = (2u-v, 2v-u, -u-v).
```

Conversely, differences between coordinates of `3Px` recover `3u` and `3v`.
Over the rationals or reals, this is the familiar rank-two projection with
kernel `span(1,1,1)`. The formal file proves the factorization, contrast
identities, and complete common-shift equivalence directly; it does not invoke
an unimplemented matrix-rank library.

## Normalization order changes the hidden direction

“Common activity” above refers to the input stage of `P`. It need not be a
literal common shift in raw activity units. If a **fixed** cell normalizer is

```text
z_ic = (r_ic - μ_c) / σ_c,     σ_c > 0,
E_i = P z_i,
```

then replacing `r_ic` by `r_ic+h_i·σ_c` adds the same `h_i` to all normalized
coordinates and leaves `E_i` unchanged. With unequal cell scales, adding the
same raw `h_i` to every cell generally does not leave `E_i` unchanged. The Lean
counterexample uses preprocessing gains `(1,2,3)`: a common input shift of one
changes the centered numerator from `(0,0,0)` to `(-3,0,3)`.

The more general checked affine statement is conditional and explicit. For
fixed gains `g_c`, offsets `b_c`, and a proposed input shift `δ_c`, if
`g_1δ_1=g_2δ_2=g_3δ_3=h`, then

```text
P [g ⊙ (x+δ) + b] = P [g ⊙ x + b].
```

Equal gains therefore preserve literal common-input-shift invariance; offsets
may differ. This theorem does not assume a normalization order for the source.

If column means and SDs are **recomputed from a whole raw panel**, arbitrary
element-specific shifts can change those fitted normalizers and violate the
normalized columns' mean/variance constraints. The row-centering theorem does
not prove that every shifted normalized panel is itself obtainable from that
full preprocessing pipeline. In particular, it does not establish arbitrarily
large hidden variation while also fixing unit column variances. These
constraints nevertheless do **not** restore uniqueness: the following exact
pair respects them.

A separate algebraic ambiguity remains for ordinary recomputed z-scores:
positive per-cell affine changes `r'_ic=a_c r_ic+b_c` give
`μ'_c=a_cμ_c+b_c` and `σ'_c=a_cσ_c`, hence the same z-scores. This conditional
calculation explains why z-scores alone do not fix raw cell-specific origins
and scales. The Lean file does not formalize estimation of means/SDs or this
whole-panel z-score calculation.

## A complete-panel counterexample with identical column moments

Consider two hypothetical four-element, three-cell input panels:

| Element | X, cell 1 | X, cell 2 | X, cell 3 | Y, cell 1 | Y, cell 2 | Y, cell 3 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 2 | 1 | 0 | 0 | −1 | −2 |
| 2 | 0 | 2 | 1 | −2 | 0 | −1 |
| 3 | 1 | 0 | 2 | −1 | −2 | 0 |
| 4 | −3 | −3 | −3 | 3 | 3 | 3 |

Every column in both panels has sum zero and sum of squares 14. Consequently
both panels have the same recomputed zero column means and the same positive
cell SDs. For population variance, that SD is `σ=√(14/4)`; for sample variance,
it is `σ=√(14/3)`. The denominator is inside the square root. Divide either
panel by the corresponding common `σ`: every column then has mean zero and
unit variance under that same convention.

Yet `Y-X` is the element-specific common shift vector `(-2,-2,-2,6)`, applied
to all three cells. Before column standardization, both panels center to

```text
P X = P Y = [( 1, 0,-1),
             (-1, 1, 0),
             ( 0,-1, 1),
             ( 0, 0, 0)].
```

The column-standardized panels therefore also have identical specificity:
`P(X/σ)=P(Y/σ)=(PX)/σ`. Their common activities differ: the row means of X are
`(1,1,1,-3)`, whereas those of Y are `(-1,-1,-1,3)`; after normalization, divide
these vectors by `σ`. The observed specificity is nonconstant, so this is not
only an all-zero-output example.

This is a counterexample to exact recovery even when both entire candidate
panels obey the same column means and variances, and even when the normalizers
are recomputed from each panel. It establishes **nonuniqueness**, not arbitrarily
large normalized shifts while those moments remain fixed. It also does not
assert that X, Y, or this preprocessing order describes the source workbook.

Lean checks the four-element lengths, exact column sums and square sums,
element-common shifts, distinct row totals, identical nonconstant `center3`
outputs, and impossibility of an exact decoder within this fixed moment class.
The real-number square-root normalization and its cancellation above are the
explicit algebraic interpretation of those checked integer facts; Lean does
not implement real variance or square roots in this file.

## Rescaled scores and an observed rank-two plane

Any fixed downstream function of centered specificity receives the same input
for equivalent rows, so it cannot restore the discarded scalar. This result
does not require the downstream function to be linear or invertible.

For a fixed invertible diagonal cell rescaling `D` after centering,
`E=D P z`, the scores satisfy a weighted relation

```text
(D⁻ᵀ 1)ᵀ E = 0.
```

Their unweighted sum can be nonzero. Under the additional fixed pre-centering
normalization above, the raw hidden direction is `diag(σ)·1`. These equations
describe one possible model family. A fitted weighted plane in a workbook
does not identify `D`, establish that this was its historical computation,
recover `μ` or `σ`, or recover absolute transcription. Multiple processing
orders or constructions can yield a rank-two observation. Good prediction of
such observations establishes a relationship to those observations under the
evaluation design; it does not invert an unverified source transform.

## Machine-checked scope and reproduction

[ContextIdentifiability.lean](../proofs/ContextIdentifiability.lean) uses Lean
4.19 and `Std`, with no external algebra package. It defines integer rows and
the division-free map `center3(x)=3x-sum(x)·1`. This is exactly three times ideal
centering; **no truncated integer division** is used. Every finite rational
panel and rational shift has a common-denominator encoding, so the integer
identities supply exact rational witnesses. Arbitrary real-number operations,
floating-point programs, biological truth, and empirical source provenance are
outside the checked model.

The file proves:

- Common-shift invariance for a row and an arbitrary per-element shift vector.
- Exact characterization of all equal-specificity rows.
- The two-contrast factorization and preserved coordinate differences.
- Persistence of equivalence after arbitrary downstream postprocessing.
- The conditional fixed cell-affine result and an unequal-gain counterexample.
- Zero specificity with nonzero input, arbitrarily distant equal-observation
  witnesses, and impossibility of an exact input decoder.
- The exact numerical illustration above.
- Two distinct complete four-element panels with identical zero column sums
  and square sums 14, nonconstant equal specificity, and no exact decoder even
  when restricted to that moment class.

Run from the repository root:

```sh
source scripts/env.sh
lean proofs/ContextIdentifiability.lean
```

All proofs are checked by Lean's kernel using explicit definitions, rewriting,
integer arithmetic, or `decide`. There are no `axiom`, `sorry`, or
`native_decide` declarations. This establishes the stated mathematical
implications under their assumptions; it does not verify the Rust/Python
pipeline, choose a biological mechanism, or resolve the source's formula.
