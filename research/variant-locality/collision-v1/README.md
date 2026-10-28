# Non-reverse-complement collisions in the fixed representation

Among 3,072 directed interior five-base SNV classes in the existing 200nt assay
panel, **448 classes** share their exact 33-component feature change with at
least one class outside their own reverse-complement orbit. These classes account
for **5,687,716 of 32,502,288 panel-slot edit opportunities (17.4994%)**. The
zero-change group contributes 64 classes and 745,052 slots. Excluding it, **384
classes and 4,942,664 slots (15.2071%)** still have a nonzero feature change
that coincides with that of a different reverse-complement orbit.

Each slot is a possible alternate base at one interior position in an assay
record: 55,276 records × 196 positions × three alternate bases. This is not a
count of unique genomic alleles, an observed human mutation frequency, or a
measurement of regulatory effects. Group membership means a matching *class*
exists somewhere in the enumerated panel; it does not mean a paired experiment
was conducted, or that both members have biological effects.

For example, changing the center of `AAAAC` from A to C gives `AACAC`; changing
the center of `CAAAA` from A to C gives `CACAA`. These edits are not reverse
complements of each other, yet their integer counts of all 32 canonical
three-base words change identically, as does their GC count. Since fixed-length
200nt models divide these counts by fixed denominators, every fixed linear
predictor using these features must predict the same signed effect for both
edits, regardless of surrounding sequence and fitted weights. The conclusion
concerns model representation, not measured biological equivalence. A nonlinear
predictor can react to the different baseline features even if delta vectors
are identical.

## Definitions and complete distributions

A directed class consists of the original five-base neighborhood and an
alternate center base. We identify it with its reverse-complemented directed
class to form an RC orbit; 3,072 classes form 1,536 orbits. We group by exact
signed feature delta and count groups containing more than one orbit. Allele
swapping reverses the sign of the delta; it is **not** an equivalence rule in
this calculation. Zero-delta swaps happen to enter the zero group because their
vectors are equal. Actual inserted flanks cannot affect an interior three-mer
feature delta, but can affect biology.

Of the 1,405 exact vectors, 1,312 contain one RC orbit, 88 contain two, four
contain four and one (the zero vector) contains 32. Thus **93 vectors** contain
multiple RC orbits; 92 have nonzero deltas. Counts sum over all original panel
records without deduplication. Group counts weight each class by the number of
panel slots in the pinned source, not equally by class.

## Evidence and replay

[analyze.py](analyze.py) authenticates the published panel class table by SHA-256,
recomputes each of its 3,072 delta vectors from five-base strings, checks the
RC and signed allele-swap identities, checks totals, and groups the exact integer
vectors. It consumes neither model coefficients nor outcomes. The source's
panel-slot counts come from the earlier [panel scan](../panel-scan-v1/README.md);
this computation does not independently rebuild the original panel scan.
An [independent full-sequence review](../pair-geometry-review-v1/README.md)
reconstructs all vectors and collision groups from the pinned class table.

[summary.json](summary.json) contains totals and definitions;
[vectors.json](vectors.json) lists all 1,405 vectors, RC-orbit representatives,
and weighted counts; [collision-classes.csv](collision-classes.csv) lists all
448 directed classes and their vectors. [manifest.json](manifest.json) records
the exact input, code and output hashes. Replay into a new directory:

```bash
python3 research/variant-locality/collision-v1/analyze.py --out runs/collision-v1-replay
```

The result is an exact property of this input and model feature representation.
It does not establish a biological mechanism, generalize to other features or
lengths without analysis, or claim global scientific novelty.
