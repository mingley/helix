# Independent checks of pair geometry and edit collisions

From the repository root, run:

```bash
python3 research/variant-locality/pair-geometry-review-v1/verify.py
python3 research/variant-locality/pair-geometry-review-v1/verify_collisions.py
```

The first program independently generates the complete set of 2,880 nearby
two-edit classes, forms all four full 200nt sequences with a different flank,
and computes their 95,040 feature-count interaction entries. It compares
11,520 ensemble predictions using exact rational arithmetic on the saved five
binary64 coefficient vectors; all eight published exact-rational extrema also
match. The [receipt](receipt.json) records the input hashes and results.

The second recomputes each of 3,072 single-edit vectors from full 200nt
sequences, reconstructs reverse-complement orbits and groups equal deltas.
Its [receipt](collision-receipt.json) verifies the 93 collision groups and
their panel-slot totals. Slot weights are inherited from the earlier panel
scan; this check does not independently reconstruct the panel. Both scripts
read no experimental outcomes and overwrite only their own receipt files.
They establish numerical agreement about the representation, not biology.
