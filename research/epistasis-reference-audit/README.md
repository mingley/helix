# S18 design windows: genomic-reference audit

This is an **outcome-free design audit**, not a frozen interaction cohort or a
measurement of biological epistasis. It checks the 7,542 windows previously
accepted by the [S18 four-sequence audit](../epistasis-source-audit/README.md)
against a pinned genome reference and assembly mapping. S18 is the authors'
designed insert table, not a table of measured activities. This audit does not
read an S20 response workbook. The article, its methods, and prior Siraj
single-variant outcomes have already been inspected; see the previous audit and
[endpoint semantics review](../epistasis/source/siraj-s20-semantics.md) for that
exposure and unresolved measurement questions.

## Rule and inputs

The previous audit establishes the literal RR, AR, RA and AA sequence square,
alleles and unsigned separation. Here, each variant ID is interpreted as a
one-based **hg19** coordinate, consistent with the authors' source processing
and the [earlier assembly review](../variant-transfer/source/README.md). For
zero-based insert edit indices `i1` and `i2`, both alleles must give the same
nonnegative hg19 insert start:

```text
start0 = variant1_position1 - 1 - i1 = variant2_position1 - 1 - i2
```

An equal absolute distance alone would also allow reversed positions; this
signed equality rules that out. Each source interval is `[start0, start0+200)`.
The entire interval must map through one aligned block of the pinned UCSC
hg19-to-hg38 chain. All coordinate/strand mappings are deduplicated before
requiring exactly one. That mapping must be to a primary hg38 chromosome, and
the complete RR sequence must match exactly the unmasked hg38 reference
sequence fetched in the mapped strand. We do not search nearby coordinates,
choose a reference-matching mapping among ambiguous mappings, clip a block,
impute a base, or rescue a failed insert.

Inputs are the unchanged publisher S18 XLSX and prior eligibility report, the
repository's UCSC `hg19ToHg38.over.chain.gz`, and the frozen UCSC `hg38.2bit`
reference installed outside Git. SHA-256 hashes for every input, the reused
mapping/reader implementations, this audit implementation, and the complete
ledger are bound in [summary.json](summary.json). The reference is the same
835,393,456-byte file used in the earlier [single-variant transfer source
preparation](../variant-transfer/source/README.md). This audit reuses that
mapping code and twoBit reader, so agreement with earlier work using those
helpers is **not** independent verification.

The parser opens the S18-only design workbook, verifies the workbook, worksheet
and headers, and checks all four source row identities and literal sequences for
each prior-eligible square before mapping. The XLSX parser may instantiate other
S18 design rows/shared strings; S18 contains oligonucleotide designs, not S20
values. All earlier rejected groups remain in the previous audit and are not
reclassified here. Source or helper hash mismatches and any disagreement between
the previous report and source rows raise errors rather than silently changing
the population. Production gates use explicit exceptions, including under
Python's optimized mode.

## Results

| Strict reference status | S18 design windows | Distinct pair identifiers represented |
| --- | ---: | ---: |
| Full 200-base match | 7,480 | 2,140 |
| Full sequence mismatch | 58 | 18 |
| No unique whole-block mapping available: zero candidates | 4 | 1 |
| Signed start disagreement, multiple mappings, or nonprimary mapping | 0 | 0 |
| Total prior-eligible windows | 7,542 | 2,156 |

The pair counts in individual rows are not additive: a pair can have both
matching and nonmatching windows. Of the 2,156 pair identifiers, 2,140 have at
least one exact-match window and 2,137 have all their previously eligible
windows matching. All 7,542 windows passed signed start equality. Among the 58
full-sequence mismatches, 34 differ at one base and 24 at two bases; 44 windows
have a mismatch at one or both edited positions. A separate
[hg19 reference check](../epistasis-hg19-audit/README.md) finds that all 7,542
reference inserts match hg19 at their declared intervals, including these 58 and
the four without a whole-block mapping. This is consistent with differences
between the reference builds; it does not establish the source-generation
history or prove that a publisher sequence or variant identity is erroneous.

[ledger.jsonl](ledger.jsonl) preserves one record for every prior-eligible
window, in prior report order. It contains source row numbers, both edit indices,
the two computed hg19 starts, every distinct whole-block mapping candidate, the
RR sequence hash, a status, and every base difference when a unique primary
mapping is compared. A mismatch's `hg38` value is the target base *oriented to
the source insert*; it may be reverse-complemented relative to forward hg38.
A null mismatch list means no full-sequence comparison was made. No failed
window was silently deleted. The previous audit remains the authoritative
record of design-only eligibility; this additional check is not a replacement
for a source-specific protocol or an outcome-filtered cohort.

## Reproduction and limits

From the repository root, with the pinned research environment and reference:

```bash
source scripts/research-env.sh
python -m unittest discover -s research/epistasis-reference-audit -p 'test_*.py'
python research/epistasis-reference-audit/audit.py \
  --reference "$HELIX_TOOLS_DIR/reference/hg38.2bit" \
  --out /tmp/helix-s18-reference-replay
cmp research/epistasis-reference-audit/ledger.jsonl /tmp/helix-s18-reference-replay/ledger.jsonl
cmp research/epistasis-reference-audit/summary.json /tmp/helix-s18-reference-replay/summary.json
```

The destination must not exist. `scripts/setup_reference.sh` can acquire the
reference and checks its size, MD5, and SHA-256. Seven synthetic tests cover
forward/reverse mapping, signed disagreement, negative starts, partial blocks,
duplicate versus genuinely distinct mappings, nonprimary and full-reference
mismatches, and malformed squares. The `hg38.2bit` file is large and is not
stored in this repository.

A single-block chain rule is deliberately conservative: an otherwise valid
sequence can fail when its mapping crosses a chain gap or is not represented in
this chain. Exact agreement is assembly consistency at the declared coordinate,
not genome-wide uniqueness, construct measurement, assay replication, or proof
of endogenous function. It also does not resolve S20's workbook-to-upstream
linkage, normalization, sample support or replicate covariance. Those gates and
outcome-free cohort, dependence and inference rules remain to be settled before
a biological interaction test.

An [independent implementation](../epistasis-reference-independent/README.md)
reads the original XLSX, chain and twoBit files without the primary helpers. It
reproduces all 7,542 complete ledger records and all status counts; it shares
the same source data and declared criteria and does not establish measurement
or biological replication.
