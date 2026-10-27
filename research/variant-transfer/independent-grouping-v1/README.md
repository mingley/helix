# Independent outcome-free component reconstruction

A separate C++ implementation reconstructs the complete external DNA partition
and prior-exposure exclusions from the 332,719 reference-validated Siraj variant
pairs. It does not import the primary grouper, its labels or its union edges.
No reporter activity or allelic-effect response source is accepted.

The completed independent run yields **1,171 components**. Prior-exposure rules
leave **2,430 variants in 174 components**. Its retained-ID order SHA-256 is
`7f613f2d41172ef55b33804af0ff851cabb41215a52eede7348d951bc77d8d30`.
Independent comparison passes **4,689 categories with zero failures**: all
1,171 complete partitions and every one of the 332,719 variant assignments,
retention/exclusion flags and reason classes match the primary implementation.
Both complete artifact trees and source bindings are authenticated. These remain
DNA-only eligibility results, separate from source-cohort and outcome gates.

## Independent computation

`run.py` directly reads the checksum-bound full mapped-pair JSONL and every
prior-inventory member. It validates literal 200 nt alleles, the centered SNV,
source IDs, GRCh38 interval types and old DNA/RC hashes. The original full
pre-exclusion input is used; the earlier coordinate-filtered subset is excluded
as an input. All 55,276 old Agarwal cores remain in DNA screening, including
valid alternate-contig cores. Only the new Siraj records require primary
chromosomes.

`group.cpp` stores DNA in seven little-endian two-bit words and creates a sorted,
flat disjoint-piece index. Its representation and data structures differ from
the primary Python big-integer/dictionary implementation. It independently
computes source-parent links, **floor(start/1 Mb)** bins keyed by chromosome,
separate strict-positive interval overlaps, exact/RC identity and full-length
at-least-90% ungapped homology between any alleles on either strand.

A crossing insert does not join every touched bin. An adjacent or separated
next-bin insert stays separate unless another declared relation connects it.
Actual positive overlap across the boundary joins the pair. The implementation
and executable boundary fixtures preserve this distinction from the superseded
metadata attempt that used touched bins.

Twenty-one disjoint pieces guarantee a shared exact seed for at most twenty
mismatches. Every cross-component candidate is checked by packed Hamming distance.
Skipping a bucket already inside the current connected component removes only
redundant edges. All new-to-new components are complete before old exposure is
applied; no old locus receives a bin or unions otherwise separate new components.

For genomic guards, old intervals are expanded by 10,000 bp and merged separately
by chromosome and source collection. Closed-intersection tests on the new
interval endpoints implement the inclusive gap convention, including adjacent
intervals. This differs from the primary prefix-maximum witness search.

All explicit old 200 nt windows use a 21-piece index. The factorized old SNV
families use **22 pieces**, preserving every candidate within twenty-one
reference mismatches. At twenty-one, a recorded alternate must repair an actual
mismatch; at twenty or fewer the old reference already qualifies. Every member,
allele and strand of every retained component is cleared. One valid old hit is
sufficient to exclude a whole component, so excluded-component hit catalogs are
witnesses rather than exhaustive enumerations.

The independent candidate budget is one billion; the completed run checked
20,964,505 distinct candidate encounters under its own documented counting
scheme. The count need not equal the primary implementation's intermediate
work count. It records 1,656 new homology union edges. No candidate was discarded
to meet the budget.

S30's 41 exact/RC technical-control variant IDs mark three components, with no
otherwise retained variant affected. These flags are kept separate from the
prior-exposure result until the final source-cohort binding. This computation
does not infer a control class from variant-name patterns.

## Artifacts and reproduction

[run-v1/bindings.json](run-v1/bindings.json) records all input hashes, protocol,
algorithm, prepared TSV fingerprints, compiler version/flags and executable SHA.
[run-v1/manifest.json](run-v1/manifest.json) hashes complete components, every
variant assignment/exclusion reason, observed union/old-hit witnesses, compiler
output, execution log and implementation snapshots. The original source JSONL
and prior inventory remain the complete input data; temporary TSVs are only a
bound transport representation.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s research/variant-transfer/independent-grouping-v1 -p test_grouping.py -v
PYTHONDONTWRITEBYTECODE=1 python3 \
  research/variant-transfer/independent-grouping-v1/run.py run \
  --out /tmp/independent-grouping-review
```

Choose a fresh output directory. A C++20 compiler (`g++`) and Python standard
library suffice; no new numerical dependency is added. Fourteen executable
counterexamples cover the bin/overlap boundary, chromosome namespace, parent
links, 20/21 mismatch boundary, reverse strands, inclusive genomic guard,
whole-component propagation, factorized SNV repair, technical controls and old
nodes remaining separate.

The [completed comparison](comparison-v1/comparison.json) authenticates both
complete artifact manifests and compares
all canonical component memberships, every variant's exclusion status/reason
classes and source bindings. Component IDs are SHA-256 of newline-joined sorted
variant IDs, with no trailing newline. Independent numerical reproduction,
source-cohort binding and biological independence remain different claims.

The independent run manifest SHA-256 is
`e25231628756da32f961f6553999a58132353b9e223e6e3df7f467f817ede33b`.
The final comparison SHA-256 is
`93f2b9d82c5b45d9dbb6344a350e48f24ac99ad0ce496d8993b813fb8d1faa60`.
[verification.json](verification.json) records the source and validation checks.
[witness-history.json](witness-history.json) retains exact union-edge differences
and both old-DNA witness choices, with the coordinate-witness representation
difference explained. Witness traversal order is not selected to force agreement.

The immutable run summary still says its subsequent comparison was pending;
the separate final comparison receipt supplies the completed result. Before that
first comparison, the verifier added full archive authentication to comparison
mode. Every execution/helper function is unchanged by AST comparison, and the
original execution snapshot and final comparison implementation are both retained.
