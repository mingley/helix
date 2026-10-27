#!/usr/bin/env python3
"""Necessary-condition support bound using frozen start bins plus true overlaps.

This omits mandatory homology and parent relations, so surviving support is an
upper bound rather than a final cohort. Additional unions/exposure flags cannot
restore an excluded component. No outcome cell is read.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path


def components(rows, width=1_000_000):
    parents = {}

    def find(x):
        parents.setdefault(x, x)
        while parents[x] != x:
            parents[x] = parents[parents[x]]
            x = parents[x]
        return x

    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parents[b] = a

    by_chrom = defaultdict(list)
    for r in rows:
        m = r["hg38"]
        node = (m["chromosome"], m["start"] // width)
        find(node)
        by_chrom[m["chromosome"]].append((m["start"], m["end"], node))
    for spans in by_chrom.values():
        highest_end = -1
        representative = None
        for start, end, node in sorted(spans):
            if start < highest_end:
                union(node, representative)
                highest_end = max(highest_end, end)
            else:
                representative = node
                highest_end = end
    groups = defaultdict(list)
    for r in rows:
        m = r["hg38"]
        groups[find((m["chromosome"], m["start"] // width))].append(r)
    return groups


def self_test():
    def r(start, end):
        return {"hg38": {"chromosome": "chr1", "start": start, "end": end}}
    # A crossed bin boundary alone must not union both complete bins.
    assert len(components([r(95, 105), r(110, 120)], 100)) == 2
    assert len(components([r(95, 105), r(104, 120)], 100)) == 1
    assert len(components([r(95, 105), r(105, 120)], 100)) == 2
    assert len(components([r(1, 10), r(90, 99)], 100)) == 1
    assert len(components([r(95, 205), r(104, 110), r(204, 210)], 100)) == 1
    print("5 frozen-coordinate grouping controls passed")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--ledger", type=Path)
    p.add_argument("--out", type=Path)
    args = p.parse_args()
    if args.self_test:
        self_test()
        return
    if args.ledger is None or args.out is None:
        p.error("--ledger and --out required")
    invalid = {"no_whole_insert_single_block_mapping", "ambiguous_whole_insert_chain_mapping",
               "mapped_nonprimary_hg38_contig", "full_reference_insert_not_exact_hg38_match",
               "not_primary_chromosome_ACGT_200nt_SNV_pair", "pair_does_not_encode_ID_ref_alt_at_index99"}
    rows = []
    with gzip.open(args.ledger, "rt") as f:
        for line in f:
            row = json.loads(line)
            if "hg38" in row and not set(row["reasons"]) & invalid:
                rows.append(row)
    groups = components(rows)
    bad = {key for key, values in groups.items()
           if any(any(reason.startswith("within_inclusive_10000bp_") for reason in row["reasons"])
                  for row in values)}
    good = [(key, vals) for key, vals in groups.items() if key not in bad]
    result = {
        "schema_version": 1, "phase": "metadata-only necessary-condition upper bound, not final homology cohort",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "coordinate_rule": "floor(start/1000000) bins plus positive interval-overlap unions; never touched-bin unions",
        "input_ledger_sha256": hashlib.sha256(args.ledger.read_bytes()).hexdigest(),
        "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "full_reference_validated_variants": len(rows), "coordinate_only_components": len(groups),
        "coordinate_only_bad_components": len(bad), "upper_bound_eligible_components": len(good),
        "upper_bound_eligible_variants": sum(len(vals) for _, vals in good),
        "omitted_relations": ["source parents", "exact/RC sequence", "full90% homology", "S30 technical controls"],
        "logic": "Additional mandatory unions/exposure flags can merge or remove survivors; they cannot restore excluded components.",
        "good_components": [{"chromosome": key[0], "start_bin": key[1], "size": len(vals)}
                            for key, vals in sorted(good)],
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "good_components"}, indent=2))


if __name__ == "__main__":
    main()
