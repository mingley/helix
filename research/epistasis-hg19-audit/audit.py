#!/usr/bin/env python3
"""Compare all sequence-eligible S18 RR inserts with the frozen hg19 reference.

This is a design-only analysis: it never reads the S20 outcome workbook.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys

import openpyxl

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

SOURCE = ROOT / "research/epistasis-source-audit/siraj2026-S18.xlsx"
LEDGER = ROOT / "research/epistasis-source-audit/s18-design-audit.json"
SOURCE_SHA256 = "92929b4cef82797237c4c2be770feb744daa97255ac6dbfe3d53d5a34c996fbb"
LEDGER_SHA256 = "70dc4a6bd1813ea879b9699565cf75e8e9c21bff9ea1af5a46c25c1ad4686b0f"
HG19_SHA256 = "797bdb4391cbe23bd30268d41d0d86676fd7b36b7f8e76797b1882d585766fa3"
HG19_BYTES = 816241703
READER_SHA256 = "ffbb77304eb324b2028cf6fb29ca8e0ea4b6a1e9daaf77ee42ce4a8f629847a2"
HEADERS = ("Variant 1", "Variant 2", "Distance", "Center Variant", "Window", "Diplotype", "Oligo")
VARIANT = re.compile(r"(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)):([1-9][0-9]*):([ACGT]):([ACGT])\Z")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_hash(path, expected):
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"SHA-256 mismatch for {path}: {actual}")


def interval_for_group(group):
    """Require both one-based variant coordinates to imply the same interval."""
    a = VARIANT.fullmatch(group["variant1"])
    b = VARIANT.fullmatch(group["variant2"])
    if not a or not b or a[1] != b[1]:
        raise ValueError("invalid or different-chromosome variant IDs")
    if a[3] == a[4] or b[3] == b[4]:
        raise ValueError("variant is not a substitution")
    i, j = group["position1_zero_based"], group["position2_zero_based"]
    if type(i) is not int or type(j) is not int or not (0 <= i < 200 and 0 <= j < 200):
        raise ValueError("invalid insert edit offsets")
    start = int(a[2]) - 1 - i
    if start < 0 or start != int(b[2]) - 1 - j or abs(i - j) < 3:
        raise ValueError("variant coordinates do not imply one valid forward insert")
    if group["distance"] != abs(int(a[2]) - int(b[2])):
        raise ValueError("ledger distance disagrees with variant coordinates")
    return a[1], start, start + 200, i, j


def check_rr_row(group, row):
    """Bind an RR row to its ledger identity and both stated source alleles."""
    chrom, start, end, i, j = interval_for_group(group)
    if len(row) != 7 or any(isinstance(v, str) and v.startswith("=") for v in row):
        raise ValueError("malformed row or formula in S18")
    v1, v2, distance, center, window, label, oligo = row
    if (v1, v2, distance, center, window, label) != (
        group["variant1"], group["variant2"], group["distance"],
        group["center"], group["window"], "RR"
    ):
        raise ValueError("RR row identity disagrees with the frozen design ledger")
    if not isinstance(oligo, str) or len(oligo) != 200 or set(oligo) - set("ACGT"):
        raise ValueError("RR insert is not 200 A/C/G/T bases")
    if oligo[i] != v1.split(":")[2] or oligo[j] != v2.split(":")[2]:
        raise ValueError("RR insert alleles disagree with the variant IDs")
    return chrom, start, end, i, j, oligo


def differences(source, reference, start, first, second):
    if len(source) != 200 or len(reference) != 200:
        raise ValueError("both sequences must be exactly 200 bases")
    return [
        {"offset_zero_based": offset, "position_one_based": start + offset + 1,
         "source_base": s, "hg19_base": r,
         "site": "variant1" if offset == first else "variant2" if offset == second else "background"}
        for offset, (s, r) in enumerate(zip(source, reference)) if s != r
    ]


def read_rows(groups):
    rows = {}
    for group in groups:
        row = group["rows"]["RR"]
        if type(row) is not int or row < 3 or row in rows:
            raise ValueError("invalid or duplicate RR row in ledger")
        rows[row] = group
    workbook = openpyxl.load_workbook(SOURCE, read_only=True, data_only=False)
    seen = {}
    try:
        if workbook.sheetnames != ["S18.Variant Pair Oligos"]:
            raise ValueError("unexpected S18 worksheet")
        iterator = workbook.active.iter_rows(values_only=True)
        next(iterator)
        if tuple(next(iterator)) != HEADERS:
            raise ValueError("unexpected S18 headers")
        for row_number, row in enumerate(iterator, start=3):
            if row_number in rows:
                seen[row_number] = check_rr_row(rows[row_number], row)
    finally:
        workbook.close()
    if set(seen) != set(rows):
        raise ValueError("one or more ledger RR rows are absent from S18")
    return seen


def audit(reference):
    require_hash(SOURCE, SOURCE_SHA256)
    require_hash(LEDGER, LEDGER_SHA256)
    require_hash(ROOT / "scripts/twobit_reference.py", READER_SHA256)
    from twobit_reference import TwoBitReader

    if reference.stat().st_size != HG19_BYTES:
        raise ValueError("unexpected hg19 twoBit size")
    require_hash(reference, HG19_SHA256)
    ledger = json.loads(LEDGER.read_text())
    if ledger["source_sha256"] != SOURCE_SHA256:
        raise ValueError("design ledger refers to a different source")
    groups = ledger["eligible_groups"]
    if len(groups) != 7542 or ledger["status_counts_not_mutually_exclusive"]["eligible"] != len(groups):
        raise ValueError("unexpected eligible design cohort")
    keys = [(g["variant1"], g["variant2"], g["center"], g["window"]) for g in groups]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate eligible design key")
    rows = read_rows(groups)
    records = []
    with TwoBitReader(reference) as hg19:
        for group in groups:
            number = group["rows"]["RR"]
            chrom, start, end, i, j, oligo = rows[number]
            ref = hg19.fetch(chrom, start, end)
            diff = differences(oligo, ref, start, i, j)
            records.append({
                "variant1": group["variant1"], "variant2": group["variant2"],
                "center": group["center"], "window": group["window"],
                "source_rr_row_1based": number, "hg19_chromosome": chrom,
                "hg19_start": start, "hg19_end": end,
                "source_rr_sha256": hashlib.sha256(oligo.encode("ascii")).hexdigest(),
                "differences": diff, "exact_match": not diff,
            })
    counts = Counter("exact_match" if r["exact_match"] else "mismatch" for r in records)
    summary = {
        "source_sha256": SOURCE_SHA256, "design_ledger_sha256": LEDGER_SHA256,
        "hg19_reference_sha256": HG19_SHA256, "hg19_reference_bytes": HG19_BYTES,
        "windows_checked": len(records),
        "distinct_pair_ids_checked": len({(r["variant1"], r["variant2"]) for r in records}),
        "exact_match_windows": counts["exact_match"], "mismatch_windows": counts["mismatch"],
        "mismatch_pair_ids": len({(r["variant1"], r["variant2"]) for r in records if not r["exact_match"]}),
        "mismatch_base_occurrences": sum(len(r["differences"]) for r in records),
        "coordinates": "0-based half-open insert intervals; 1-based mismatch positions",
        "scope": "S18 RR designs only; no outcome workbook or observations read",
    }
    return summary, records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path, help="hash-bound UCSC hg19.2bit")
    parser.add_argument("--out", required=True, type=Path, help="new output directory")
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("output directory already exists; refusing to overwrite")
    summary, records = audit(args.reference)
    args.out.mkdir(parents=True, exist_ok=False)
    with (args.out / "windows.jsonl").open("w", encoding="utf-8", newline="\n") as f:
        for record in records:
            f.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    manifest = {name: {"bytes": (args.out / name).stat().st_size,
                       "sha256": sha256(args.out / name)} for name in ("summary.json", "windows.jsonl")}
    manifest["implementation_sha256"] = sha256(Path(__file__))
    manifest["twobit_reader_sha256"] = sha256(ROOT / "scripts/twobit_reference.py")
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
