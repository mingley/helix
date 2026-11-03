# Independent check of S18 feature collisions

This outcome-free check reconstructs the
[four-construct feature audit](../epistasis-feature-collisions-v1/README.md)
directly from the original S18 design workbook. It reads S18 as ZIP/XML using
the Python standard library and independently reclassifies the full 35,616
source rows into 7,542 eligible and 1,362 ineligible pair/window groups. It
compares every eligible group's labels, edit positions, row numbers and source
identifiers to the earlier design audit, and every rejected group's key and
row numbers. It does not open the S20 workbook or read reporter outcomes.

For each eligible group, the verifier computes the complete feature vector
of all four 200-base oligos: it counts each oriented 3-mer separately, adds the
counts of its reverse complement, and counts GC bases. This differs from the
primary implementation's canonicalization of each window before counting.
It checks the separated-edit feature parallelogram, the two 33-integer feature
deltas, all within-square feature-equality pairs, global reverse-complement
status, and every field of the 7,542-record primary ledger. It classifies
collisions from equality patterns, independently of the primary code's
delta-based classification, and recomputes every numerical summary category.

Both implementations find **282 windows across 78 pair identifiers with feature collisions**. Their
506 colliding construct-pair occurrences represent **480 distinct unordered
exact DNA-sequence pairs**; none are global reverse complements. The verifier
also confirms the class breakdown of 104 first-edit-blind windows, 120
second-edit-blind windows, 46 equal nonzero single-edit deltas, 12 opposite
nonzero deltas, and 7,260 windows with no collision. These are design
opportunities, not observed functional effects or independent biological
replicates.

The hg19 and hg38 statuses in the feature ledger are compared to the
hash-bound, previously independently checked reference ledgers; this checker
does not independently repeat their genome mapping. The original source,
prior design report, reference ledgers, feature ledger, summary, primary
manifest and verification code are bound by SHA-256 in the
[receipt](receipt.json). It also checks that the primary manifest binds its
ledger, summary, and implementation. A hash establishes identity within these
checks, not authenticity or biological validity.

Run from the repository root:

```bash
python3 -m unittest discover -s research/epistasis-feature-collisions-independent -p 'test_*.py'
python3 research/epistasis-feature-collisions-independent/verify.py \
  --receipt /tmp/epistasis-feature-collisions-independent.json
cmp research/epistasis-feature-collisions-independent/receipt.json \
  /tmp/epistasis-feature-collisions-independent.json
```

The tests include all six possible collision classes on synthetic full-length
sequences, an adjacent-edit counterexample, DNA reverse-complement invariance,
and source and ledger tampering under `python -O`. Errors occur before an
existing receipt is overwritten. These computational checks do not establish
outcome availability, reproducible activity differences, a causal regulatory
mechanism, or scientific priority.
