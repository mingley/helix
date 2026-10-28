# Independent check of the S18 sequence-design audit

`verify.py` reads the original S18 XLSX directly as ZIP/XML using only the Python
standard library. It independently checks design eligibility by comparing all
four sequences at every nucleotide and requiring the prescribed four genotype
patterns at exactly two varying positions. It reads **no S20 response values**.

The check verifies all 35,616 data rows, 8,904 pair/window group identities,
7,542 accepted groups and 1,362 rejected groups against the primary audit. For
every accepted group it compares variant and window identities, source row
numbers, signed allele labels, positions in the 200nt sequences and distance;
for rejected groups it compares group identities and source row numbers. It
checks label counts, the eligible distance distribution and pair-identifier
counts. It does not independently classify the detailed reason for every
rejected group, independently verify genomic assembly against a reference, or
determine whether measured responses or usable replicates exist. The two
implementations share the declared eligibility criteria and the original source.

The [receipt](receipt.json) binds the source and primary report SHA-256 hashes
and records the checks. To repeat from repository root:

```bash
python3 research/epistasis-s18-independent/verify.py research/epistasis-source-audit/siraj2026-S18.xlsx
python3 -m unittest discover -s research/epistasis-s18-independent -p 'test_*.py'
```

The source SHA-256 is
`92929b4cef82797237c4c2be770feb744daa97255ac6dbfe3d53d5a34c996fbb`.
Tests exercise valid and malformed synthetic squares; they do not establish
any biological effect. Results are design eligibility, not outcome availability.
