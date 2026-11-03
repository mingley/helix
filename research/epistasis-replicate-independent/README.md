# Independent check of the exploratory OL27 name inventory

This verification independently reproduces the
[OL27 name inventory](../epistasis-replicate-metadata/README.md) for the
previously audited S18 design windows. The 7,542 designs and their feature
collision classes are read from the hash-bound design and feature ledgers;
this check does not rerun their DNA or genomic-reference audits. The lookup
file is the separately acquired, hash-pinned `ENCFF823YWM.tsv.gz`, external to
Git. Its source and terms are described in the primary inventory.

The checker generates names by selecting the centered and other variant indices
for each genotype, independently of the primary name-construction routine. It
then uses Python's CSV TSV parser and a separate streaming checksum wrapper to
read **12,054,308** lookup rows, and recomputes every numeric field in the
primary report. It checks both occurrence under any mapping flag and occurrence
on a row with flag 0 and exactly one name. The decompressed byte count and MD5,
compressed SHA-256, prior ledgers, primary report, primary implementation and
this verification code are bound in the [receipt](receipt.json).

The two flag rules give the same overlap counts in this lookup: **256 windows
across 105 variant-pair identifiers** have all four candidate strings present;
**12 windows across five pair identifiers** also carry a previously computed
within-square feature collision. These numbers count strings, not measured
replicates or biological effects. They depend on a tentative name syntax: the
strings are **syntactic candidates**, and this check does not establish that
any string is the intended S18 construct, that its barcode is usable, or that
RNA measurements exist or share a valid normalization.

The parser instantiates the lookup's mapping/QC columns and barcode field as
text, but it interprets only the name list and overall mapping flag. No RNA,
reporter-activity or outcome table is read. Reference-ledger and feature-ledger
checks are computational provenance, not independent biological replication.

With the lookup obtained as described in the primary inventory, run:

```bash
source scripts/research-env.sh
python3 -m unittest discover -s research/epistasis-replicate-independent -p 'test_*.py'
python3 research/epistasis-replicate-independent/verify.py \
  --lookup "$HELIX_TOOLS_DIR/epistasis-replicate-metadata/ENCFF823YWM.tsv.gz" \
  --receipt /tmp/ol27-name-inventory-independent.json
cmp research/epistasis-replicate-independent/receipt.json /tmp/ol27-name-inventory-independent.json
```

Choose a new receipt path: the verifier uses exclusive file creation and does
not overwrite an existing receipt. Synthetic tests check center-specific allele
order, malformed rows, the joint flag-0/single-name requirement, partial versus
complete squares, and both an existing output and a race after the initial
output check. No full-sized external input is needed for those tests.
