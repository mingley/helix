#!/usr/bin/env python3
"""Independently compare the S18 genomic-design audit directly with hg19.

This reads only S18 design sequences, source audit reports, and the genome
references. No response or S20 workbook is read.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from verify import (ROOT, SOURCE_SHA, TwoBit, VARIANT, collect_design, require,
                    require_new_receipt, sha256, write_receipt)

HG19_SHA = "797bdb4391cbe23bd30268d41d0d86676fd7b36b7f8e76797b1882d585766fa3"


def expected_record(item, reference):
    first = VARIANT.fullmatch(item["variant1"])
    second = VARIANT.fullmatch(item["variant2"])
    start = int(first[2]) - 1 - item["position1_zero_based"]
    other_start = int(second[2]) - 1 - item["position2_zero_based"]
    require(start >= 0 and start == other_start, "inconsistent source interval")
    end = start + 200
    genomic = reference.fetch(first[1], start, end, "+")
    differences = []
    for offset, (source, genome) in enumerate(zip(item["reference"], genomic)):
        if source != genome:
            if offset == item["position1_zero_based"]:
                site = "variant1"
            elif offset == item["position2_zero_based"]:
                site = "variant2"
            else:
                site = "background"
            differences.append({"offset_zero_based": offset, "position_one_based": start + offset + 1,
                                "source_base": source, "hg19_base": genome, "site": site})
    return {"variant1": item["variant1"], "variant2": item["variant2"],
            "center": item["center"], "window": item["window"],
            "source_rr_row_1based": item["rows"]["RR"], "hg19_chromosome": first[1],
            "hg19_start": start, "hg19_end": end,
            "source_rr_sha256": hashlib.sha256(item["reference"].encode("ascii")).hexdigest(),
            "differences": differences, "exact_match": not differences}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=ROOT / "research/epistasis-source-audit/siraj2026-S18.xlsx")
    p.add_argument("--design-report", type=Path, default=ROOT / "research/epistasis-source-audit/s18-design-audit.json")
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--ledger", type=Path, default=ROOT / "research/epistasis-hg19-audit/results/windows.jsonl")
    p.add_argument("--hg38-ledger", type=Path, default=ROOT / "research/epistasis-reference-audit/ledger.jsonl")
    p.add_argument("--hg38-receipt", type=Path, default=ROOT / "research/epistasis-reference-independent/receipt.json")
    p.add_argument("--receipt", type=Path, required=True)
    args = p.parse_args()
    require_new_receipt(args.receipt)
    require(sha256(args.source) == SOURCE_SHA, "S18 input SHA-256 mismatch")
    require(sha256(args.reference) == HG19_SHA, "hg19 input SHA-256 mismatch")
    count, accepted, rejected = collect_design(args.source)
    report = json.loads(args.design_report.read_text())
    published_accepted = {(x["variant1"], x["variant2"], x["center"], x["window"]): x
                          for x in report["eligible_groups"]}
    published_rejected = {tuple(x["key"]): sorted(x["rows"]) for x in report["ineligible_groups"]}
    require(len(published_accepted) == len(report["eligible_groups"]) and
            len(published_rejected) == len(report["ineligible_groups"]), "duplicate design report keys")
    require(count == report["data_rows"] and len(accepted) + len(rejected) == report["design_groups"],
            "design row or group count differs")
    require(rejected == published_rejected and set(accepted) == set(published_accepted),
            "design eligibility differs")
    for key, item in accepted.items():
        require({k: v for k, v in item.items() if k != "reference"} == published_accepted[key],
                f"design details differ: {key}")
    old_receipt = json.loads(args.hg38_receipt.read_text())
    require(old_receipt.get("all_checks_passed") is True and
            old_receipt.get("source_sha256") == SOURCE_SHA and
            old_receipt.get("design_report_sha256") == sha256(args.design_report) and
            old_receipt.get("implementation_sha256") == sha256(Path(__file__).with_name("verify.py")) and
            old_receipt.get("ledger_sha256") == sha256(args.hg38_ledger),
            "hg38 ledger is not bound to the independently checked receipt")
    hg38_status = {}
    with args.hg38_ledger.open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            key = tuple(record[k] for k in ("variant1", "variant2", "center", "window"))
            require(key in accepted and key not in hg38_status, "unexpected hg38 ledger key")
            hg38_status[key] = record["status"]
    require(set(hg38_status) == set(accepted), "hg38 ledger missing source groups")
    reference = TwoBit(args.reference)
    seen, joint = set(), Counter()
    try:
        with args.ledger.open(encoding="utf-8") as handle:
            for number, line in enumerate(handle, 1):
                actual = json.loads(line)
                key = tuple(actual[k] for k in ("variant1", "variant2", "center", "window"))
                require(key in accepted and key not in seen,
                        f"unexpected or duplicate hg19 ledger key at line {number}")
                seen.add(key)
                expected = expected_record(accepted[key], reference)
                if actual != expected:
                    fields = [k for k in set(actual) | set(expected) if actual.get(k) != expected.get(k)]
                    raise ValueError(f"hg19 ledger line {number} differs: {key}; fields={sorted(fields)}")
                joint[(hg38_status[key], expected["exact_match"])] += 1
    finally:
        reference.close()
    require(seen == set(accepted), "hg19 ledger missing source groups")
    receipt = {"method": "Independent XLSX XML genotype patterns and on-demand hg19 2bit decoding",
               "implementation_sha256": sha256(Path(__file__).resolve()),
               "shared_reader_sha256": sha256(Path(__file__).with_name("verify.py")),
               "source_sha256": sha256(args.source), "design_report_sha256": sha256(args.design_report),
               "hg19_reference_sha256": sha256(args.reference), "ledger_sha256": sha256(args.ledger),
               "hg38_ledger_sha256": sha256(args.hg38_ledger), "hg38_receipt_sha256": sha256(args.hg38_receipt),
               "source_rows": count, "accepted_groups": len(accepted), "rejected_groups": len(rejected),
               "ledger_groups_checked": len(seen),
               "hg19_exact_match_by_hg38_status": {status: {"match": joint[(status, True)],
                                                              "mismatch": joint[(status, False)]}
                                                     for status in sorted(set(hg38_status.values()))},
               "detailed_rejection_reasons_checked": False, "outcome_values_read": False,
               "all_checks_passed": True}
    write_receipt(args.receipt, receipt)
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
