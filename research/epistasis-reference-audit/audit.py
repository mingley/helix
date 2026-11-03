#!/usr/bin/env python3
"""Reference-check prior S18 sequence squares without reading outcome data.

Reproduce from the repository root. Output is written only to a new directory;
see README.md for the source and reference contracts and scope.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys

import openpyxl

ROOT = Path(__file__).resolve().parents[2]
S18 = ROOT / "research/epistasis-source-audit/siraj2026-S18.xlsx"
PRIOR = ROOT / "research/epistasis-source-audit/s18-design-audit.json"
CHAIN = ROOT / "research/variant-transfer/source/hg19ToHg38.over.chain.gz"
MAPPER = ROOT / "research/variant-transfer/source/map_design.py"
TWOBIT = ROOT / "scripts/twobit_reference.py"
EXPECTED_HASHES = {
    "s18_workbook": "92929b4cef82797237c4c2be770feb744daa97255ac6dbfe3d53d5a34c996fbb",
    "prior_design_audit": "70dc4a6bd1813ea879b9699565cf75e8e9c21bff9ea1af5a46c25c1ad4686b0f",
    "hg19_to_hg38_chain": "5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d",
    "hg38_twobit": "1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0",
    "chain_mapper": "4bfb65df9ceb9add57bbcde62f1f6b8f005dbc63991f08ed73ce1bed8bbd0212",
    "twobit_reader": "ffbb77304eb324b2028cf6fb29ca8e0ea4b6a1e9daaf77ee42ce4a8f629847a2",
}
HEADERS = ("Variant 1", "Variant 2", "Distance", "Center Variant", "Window", "Diplotype", "Oligo")
LABELS = ("RR", "AR", "RA", "AA")
VARIANT = re.compile(r"^(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)):([0-9]+):([ACGT]):([ACGT])$")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def mapper():
    # Validate shared code before importing it. This does not claim that the
    # shared implementation independently verifies itself.
    if sha(MAPPER) != EXPECTED_HASHES["chain_mapper"] or sha(TWOBIT) != EXPECTED_HASHES["twobit_reader"]:
        raise ValueError("mapping or twoBit helper hash mismatch")
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("helix_frozen_map_design", MAPPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_sequences(entries, source):
    """Bind all four retained row identities and sequences to the S18 workbook."""
    required = {}
    keys = set()
    for index, entry in enumerate(entries):
        key = tuple(entry[k] for k in ("variant1", "variant2", "center", "window"))
        if key in keys or set(entry["rows"]) != set(LABELS):
            raise ValueError("duplicate prior window or malformed source row list")
        keys.add(key)
        for label in LABELS:
            rownum = entry["rows"][label]
            if type(rownum) is not int or rownum < 3 or rownum in required:
                raise ValueError("duplicate or invalid prior source row number")
            required[rownum] = (index, label)
    seqs = [{} for _ in entries]
    wb = openpyxl.load_workbook(source, read_only=True, data_only=False)
    try:
        if wb.sheetnames != ["S18.Variant Pair Oligos"]:
            raise ValueError("unexpected S18 sheet names")
        rows = wb.active.iter_rows(values_only=True)
        title, headers = next(rows), next(rows)
        if title[0] != "Supplementary Table 18 | Oligo sequences for variant pairs" or tuple(headers) != HEADERS:
            raise ValueError("unexpected S18 title or header")
        for rownum, row in enumerate(rows, 3):
            if rownum not in required:
                continue
            index, label = required.pop(rownum)
            entry = entries[index]
            if len(row) != 7 or row[:6] != (entry["variant1"], entry["variant2"], entry["distance"],
                                           entry["center"], entry["window"], label):
                raise ValueError(f"prior row identity disagrees with workbook at row {rownum}")
            seqs[index][label] = row[6]
    finally:
        wb.close()
    if required:
        raise ValueError("prior report refers to missing source rows")
    return seqs


def validate_square(entry, seqs):
    """Reject a malformed prior audit binding rather than repair an insert."""
    m1 = VARIANT.fullmatch(entry["variant1"])
    m2 = VARIANT.fullmatch(entry["variant2"])
    if not m1 or not m2 or m1.group(1) != m2.group(1):
        raise ValueError("prior eligible record has malformed or cross-chromosome IDs")
    chrom, n1, ref1, alt1 = m1.groups()
    _, n2, ref2, alt2 = m2.groups()
    n1, n2 = int(n1), int(n2)
    p1, p2 = entry["position1_zero_based"], entry["position2_zero_based"]
    if (n1 < 1 or n2 < 1 or type(p1) is not int or type(p2) is not int or
            not 0 <= p1 < 200 or not 0 <= p2 < 200 or p1 == p2 or
            ref1 == alt1 or ref2 == alt2 or abs(n1 - n2) < 3 or
            abs(n1 - n2) != entry["distance"] or abs(p1 - p2) != entry["distance"]):
        raise ValueError("prior eligible record has invalid alleles, positions or distance")
    if (set(seqs) != set(LABELS) or any(not isinstance(s, str) or len(s) != 200 or
                                        set(s) - set("ACGT") for s in seqs.values())):
        raise ValueError("prior eligible record lacks four literal 200nt sequences")
    rr = seqs["RR"]
    if rr[p1] != ref1 or rr[p2] != ref2:
        raise ValueError("prior eligible record disagrees with reference alleles")
    ar = list(rr)
    ar[p1] = alt1
    ra = list(rr)
    ra[p2] = alt2
    aa = list(ar)
    aa[p2] = alt2
    if seqs["AR"] != "".join(ar) or seqs["RA"] != "".join(ra) or seqs["AA"] != "".join(aa):
        raise ValueError("prior eligible record is not the stated exact four-construct square")
    return chrom, n1, n2, p1, p2, rr


def audit_one(entry, seqs, blocks, reference, helpers):
    chrom, n1, n2, p1, p2, rr = validate_square(entry, seqs)
    start1, start2 = n1 - 1 - p1, n2 - 1 - p2
    consistent = start1 == start2 and start1 >= 0
    record = {
        "variant1": entry["variant1"], "variant2": entry["variant2"],
        "center": entry["center"], "window": entry["window"],
        "source_rows": entry["rows"],
        "position1_zero_based": p1, "position2_zero_based": p2,
        "hg19": {"chromosome": chrom, "start_from_variant1": start1,
                 "start_from_variant2": start2, "end": start1 + 200 if consistent else None},
        "mapping_candidates": [],
        "reference_sequence_sha256": hashlib.sha256(rr.encode("ascii")).hexdigest(),
        "mismatches": None,
        "status": None,
    }
    if not consistent:
        record["status"] = "inconsistent_start"
        return record
    block_index = blocks.get(chrom)
    candidates = sorted({helpers.mapped_interval(block, start1, start1 + 200)
                         for block in (() if block_index is None else block_index.enclosing(start1, start1 + 200))})
    record["mapping_candidates"] = [dict(zip(("chromosome", "start", "end", "strand"), candidate))
                                    for candidate in candidates]
    if len(candidates) != 1:
        record["status"] = "no_map" if not candidates else "ambiguous_map"
        return record
    qchrom, qstart, qend, strand = candidates[0]
    if not helpers.PRIMARY.fullmatch(qchrom):
        record["status"] = "nonprimary"
        return record
    genomic = reference.fetch(qchrom, qstart, qend, strand=strand)
    if len(genomic) != 200:
        raise ValueError("reference fetch did not return exactly 200 bases")
    record["mismatches"] = [{"position_zero_based": i, "source": a, "hg38": b}
                            for i, (a, b) in enumerate(zip(rr, genomic)) if a != b]
    record["status"] = "full_reference_mismatch" if record["mismatches"] else "valid"
    return record


def run(reference_path, output):
    inputs = {"s18_workbook": S18, "prior_design_audit": PRIOR,
              "hg19_to_hg38_chain": CHAIN, "hg38_twobit": reference_path,
              "chain_mapper": MAPPER, "twobit_reader": TWOBIT}
    for name, path in inputs.items():
        if sha(path) != EXPECTED_HASHES[name]:
            raise ValueError(f"input hash mismatch: {name}")
    helpers = mapper()
    audit = json.loads(PRIOR.read_text())
    entries = audit["eligible_groups"]
    if len(entries) != audit["status_counts_not_mutually_exclusive"]["eligible"]:
        raise ValueError("prior eligible count disagrees with records")
    seqs = source_sequences(entries, S18)
    blocks = helpers.read_chain(CHAIN)
    # Do not overwrite a prior report or conflate stale outputs with this run.
    output.mkdir(parents=True, exist_ok=False)
    counts = Counter()
    pairs = defaultdict(list)
    ledger_path = output / "ledger.jsonl"
    with helpers.TwoBitReader(reference_path) as reference, ledger_path.open("w", encoding="utf-8", newline="\n") as ledger:
        for entry, sequences in zip(entries, seqs):
            record = audit_one(entry, sequences, blocks, reference, helpers)
            counts[record["status"]] += 1
            pairs[(record["variant1"], record["variant2"])].append(record["status"])
            ledger.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
    summary = {
        "schema_version": 1,
        "scope": "S18 design-only genomic consistency; no S20 outcome workbook read",
        "source_assembly": "hg19; variant IDs 1-based, insert intervals zero-based half-open",
        "target_reference": "UCSC hg38.2bit, unmasked; source-oriented via chain strand",
        "rules": ["same nonnegative hg19 start implied by both signed positions",
                  "unique deduplicated 200nt single-chain-block mapping before sequence validation",
                  "primary hg38 chromosome and exact full 200nt RR sequence match"],
        "counts_by_status": dict(sorted(counts.items())),
        "design_windows": len(entries), "pair_identifiers": len(pairs),
        "pairs_with_at_least_one_valid_window": sum("valid" in x for x in pairs.values()),
        "pairs_with_all_eligible_windows_valid": sum(all(y == "valid" for y in x) for x in pairs.values()),
        "inputs_sha256": EXPECTED_HASHES,
        "implementation_sha256": sha(Path(__file__)),
        "ledger_sha256": sha(ledger_path),
        "independent_biological_replicates_established": False,
        "outcome_availability_established": False,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("counts_by_status", "pair_identifiers", "pairs_with_at_least_one_valid_window", "pairs_with_all_eligible_windows_valid")}, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True, help="frozen hg38.2bit")
    parser.add_argument("--out", type=Path, required=True, help="new output directory")
    args = parser.parse_args()
    run(args.reference, args.out)


if __name__ == "__main__":
    main()
