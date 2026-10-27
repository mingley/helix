# Exact locality of the composition predictor

An independent, outcome-free enumeration of all 1,024 oriented DNA 5mers and
three alternative central bases gives **3,072 edits, 1,536 reverse-complement
pairs and 1,405 distinct exact feature changes**. The representation consists of
32 reverse-complement-canonical 3mer frequencies and GC fraction on a fixed
200nt sequence. No model coefficients or experimental outcomes enter this audit.

Exactly **64 edits** change none of these features. In this finite enumeration,
that occurs exactly when the mutated local 5mer is the reverse complement of its
original, for example `AAATT → AATTT`. There are 16 complementary pairs of outer
flanks and four central bases, giving 64 such edits. Their 1/48 fraction is a
property of the complete synthetic edit space, not an estimated frequency in
human variation. The rest of the 200nt sequence can remain unchanged: only the
three overlapping windows containing an interior edited base can change.

The minimum nonzero 3mer-count L1 change is 2. Including the GC count gives a
minimum joint integer L1 change of 3; in the actual normalized features the
minimum is `2/198 + 1/200 = 299/19800`, for example `AAACC → AACCC`. These norms
do not impose a minimum predicted activity change: coefficients can cancel or
ignore a nonzero feature change.

Every pair of substitutions separated by at least three nucleotides has exactly
zero mixed feature difference,
`F(double) − F(first) − F(second) + F(reference)`. No 3mer window can contain both
edited sites, and GC is additive. Therefore a fixed linear predictor using these
features cannot represent an interaction at those distances. The exhaustive
distance-three check covered 589,824 double edits. Distance-one and distance-two
feature differences were all nonzero in the respective complete synthetic
spaces, but their coefficient-weighted prediction differences can still vanish.

The complete artifact bundle is [exact-audit-v1](exact-audit-v1/summary.json).
It retains the original implementation, all edit vectors, all equivalence
classes, zero-change examples and a SHA-256/byte manifest. It also records
786,432 exhaustive two-base outer-flank checks and 3,072 full 200nt checks.
The enumeration applies to interior sites at zero-based positions 2 through
197. Edits in the first or last two bases require truncating the affected
windows; the general distance argument still applies.

Replay from the repository root with standard-library Python, choosing a new
output directory:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 research/variant-locality/exact-audit-v1/implementation.py \
  --out runs/variant-locality-exact-replay
```

The snapshot is copied byte-for-byte from the independent computation; its
implementation accepts any fresh `--out` path and contains no hardcoded output
root. The replay preserves exact integer and rational arithmetic.

This is a limitation of a feature representation, not evidence of biological
absence of interactions, measured variant effects or clinical consequences.
No new experimental outcome payload has been acquired or inspected for this
audit. Review of candidate external experiments remains limited to metadata and
schema until a separate prediction protocol is frozen. This finite calculation
does not establish global scientific novelty or independently validate the
composition predictor on another experiment.
