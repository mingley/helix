#!/usr/bin/env python3
"""Inspect public design DNA against an assembly chain and explored loci.

This is an outcome-free source gate. Whole inserts must map within one chain
block, map uniquely before sequence validation, and match the frozen hg38 DNA
exactly. It never aligns against an assay response or repairs source sequences.
Source class and assay endpoint gates are separate and remain required.
"""
import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from twobit_reference import TwoBitReader

PRIMARY = re.compile(r"chr(?:[1-9]|1\d|2[0-2]|X|Y)$")
DNA_EXPORT_SHA = "c8cbb47713f6020f1c134f10644d1e991694b78379e9f3699a589635401848df"
CHAIN_SHA = "5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d"
REFERENCE_SHA = "1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0"
PANEL_SHA = "fce6b2a1aaf5aa4d530a2e298dccb840bc7c594a05b6831fdb0c6c0f9e2a9be3"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical(s):
    return min(s, rc(s))


class Blocks:
    def __init__(self, rows):
        self.rows = sorted(rows)
        self.starts = [r[0] for r in self.rows]
        self.maxends = []
        top = -1
        for r in self.rows:
            top = max(top, r[1])
            self.maxends.append(top)

    def enclosing(self, start, end):
        i = bisect_right(self.starts, start) - 1
        while i >= 0 and self.maxends[i] >= end:
            r = self.rows[i]
            if r[0] <= start and end <= r[1]:
                yield r
            i -= 1


def read_chain(path):
    blocks = defaultdict(list)
    header = None
    with gzip.open(path, "rt") as f:
        for line in f:
            v = line.split()
            if not v:
                continue
            if v[0] == "chain":
                if len(v) != 13 or v[4] != "+" or v[9] not in {"+", "-"}:
                    raise ValueError("unsupported or malformed chain header")
                header = v
                ts, qs = int(v[5]), int(v[10])
                continue
            if header is None or len(v) not in {1, 3}:
                raise ValueError("chain block without supported header")
            size = int(v[0])
            if size <= 0:
                raise ValueError("nonpositive chain block")
            # tStart/tEnd, qName, qSize, qStrand, qStart, chainId
            blocks[header[2]].append((ts, ts + size, header[7],
                                     int(header[8]), header[9], qs, header[12]))
            ts += size
            qs += size
            if len(v) == 3:
                dt, dq = int(v[1]), int(v[2])
                if dt < 0 or dq < 0:
                    raise ValueError("negative chain gap")
                ts += dt
                qs += dq
            elif ts != int(header[6]) or qs != int(header[11]):
                raise ValueError("chain terminal coordinates disagree")
    return {k: Blocks(v) for k, v in blocks.items()}


def mapped_interval(block, start, end):
    tstart, tend, chrom, qsize, strand, qstart, _ = block
    if not tstart <= start < end <= tend:
        raise ValueError("whole interval is not in one aligned block")
    left, right = qstart + start - tstart, qstart + end - tstart
    if strand == "-":
        left, right = qsize - right, qsize - left
    if not 0 <= left < right <= qsize:
        raise ValueError("mapped coordinates out of bounds")
    return chrom, left, right, strand


def merged_guard_intervals(intervals, guard):
    out = defaultdict(list)
    for chrom, start, end in intervals:
        out[chrom].append((max(0, start - guard), end + guard))
    result = {}
    for chrom, vals in out.items():
        merged = []
        for left, right in sorted(vals):
            if merged and left <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(right, merged[-1][1]))
            else:
                merged.append((left, right))
        result[chrom] = (merged, [x[0] for x in merged])
    return result


def hits_guard(index, chrom, start, end):
    if chrom not in index:
        return False
    vals, starts = index[chrom]
    i = bisect_right(starts, end) - 1
    return i >= 0 and vals[i][1] >= start


def self_test():
    plus = (10, 20, "chr2", 1000, "+", 100, "1")
    minus = (10, 20, "chr2", 1000, "-", 100, "2")
    assert mapped_interval(plus, 12, 17) == ("chr2", 102, 107, "+")
    assert mapped_interval(minus, 12, 17) == ("chr2", 893, 898, "-")
    try:
        mapped_interval(plus, 9, 17)
    except ValueError:
        pass
    else:
        raise AssertionError("clipped insert admitted")
    a = Blocks([(0, 100, "a"), (20, 30, "b"), (40, 50, "c")])
    assert list(a.enclosing(45, 48)) == [(40, 50, "c"), (0, 100, "a")]
    guards = merged_guard_intervals([("chr1", 100, 200)], 10)
    assert hits_guard(guards, "chr1", 210, 230)
    assert not hits_guard(guards, "chr1", 211, 230)
    assert not hits_guard(guards, "chr2", 100, 200)
    assert canonical("ACG") == canonical(rc("ACG"))
    print("8 design-mapping controls passed")


def run(args):
    bindings = {
        "source_DNA_export": (args.dna, DNA_EXPORT_SHA),
        "hg19_to_hg38_chain": (args.chain, CHAIN_SHA),
        "hg38_reference": (args.reference, REFERENCE_SHA),
        "explored_Agarwal_panel": (ROOT / "datasets/context-panel/panel.json", PANEL_SHA),
    }
    for name, (path, expected) in bindings.items():
        if digest(path) != expected:
            raise ValueError(f"input hash mismatch: {name}")
    blocks = read_chain(args.chain)
    panel = json.loads(bindings["explored_Agarwal_panel"][0].read_text())["records"]
    old_sequences = {canonical(x["sequence"]) for x in panel}
    panel_intervals = [(x["chromosome"], x["start"], x["end"]) for x in panel]
    metadata_path = ROOT / "datasets/regulatory-cohort/metadata-selection.json"
    explored = json.loads(metadata_path.read_text())["assays"]
    explored_intervals = []
    for assay in explored:
        m = assay["mapping"]
        if m["assembly"] != "GRCh38":
            raise ValueError("explored mapping uses another genome assembly")
        explored_intervals.append(("chr" + m["chromosome"].removeprefix("chr"), m["start"], m["end"]))
        s = assay["target_sequence"]
        old_sequences.update(canonical(s[i:i + 200]) for i in range(len(s) - 199))
    hbb_path = ROOT / "datasets/hbb-reference/annotation.json"
    hbb = json.loads(hbb_path.read_text())["ncbiRefSeq"]
    explored_intervals.extend((x["chrom"], x["txStart"], x["txEnd"]) for x in hbb)
    panel_guard = merged_guard_intervals(panel_intervals, 10000)
    other_guard = merged_guard_intervals(explored_intervals, 10000)
    args.out.mkdir(parents=True, exist_ok=False)
    counts = Counter()
    eligible_ids = []
    with TwoBitReader(args.reference) as reference, gzip.open(args.dna, "rt") as source, \
            gzip.GzipFile(filename=str(args.out / "design-ledger.jsonl.gz"), mode="wb", mtime=0) as ledger, \
            gzip.GzipFile(filename=str(args.out / "all-reference-validated-pairs.jsonl.gz"), mode="wb", mtime=0) as all_mapped, \
            gzip.GzipFile(filename=str(args.out / "eligible-pairs.jsonl.gz"), mode="wb", mtime=0) as eligible:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames != ["row_1based", "variant", "allele_1_oligo", "allele_2_oligo"]:
            raise ValueError("unbound design-export header")
        for row in reader:
            counts["all_design_rows"] += 1
            vid, a, b = row["variant"], row["allele_1_oligo"], row["allele_2_oligo"]
            record = {"id": vid, "source_row_1based": int(row["row_1based"]), "reasons": []}
            parts = vid.split(":")
            good = (len(parts) == 4 and PRIMARY.fullmatch(parts[0]) and parts[1].isdigit()
                    and len(parts[2]) == len(parts[3]) == 1 and parts[2] in "ACGT"
                    and parts[3] in "ACGT" and parts[2] != parts[3]
                    and len(a) == len(b) == 200 and set(a + b) <= set("ACGT"))
            if not good:
                record["reasons"].append("not_primary_chromosome_ACGT_200nt_SNV_pair")
            else:
                counts["primary_chromosome_200nt_SNV_pairs"] += 1
                diffs = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
                if diffs != [99] or a[99] != parts[2] or b[99] != parts[3]:
                    record["reasons"].append("pair_does_not_encode_ID_ref_alt_at_index99")
                else:
                    chrom, pos = parts[0], int(parts[1])
                    start, end = pos - 100, pos + 100
                    record["source_hg19"] = {"chromosome": chrom, "start": start, "end": end,
                                              "variant_position_1based": pos}
                    candidates = sorted({mapped_interval(v, start, end)
                                         for v in blocks.get(chrom, Blocks([])).enclosing(start, end)})
                    record["whole_insert_chain_candidates"] = [list(x) for x in candidates]
                    if len(candidates) != 1:
                        record["reasons"].append("no_whole_insert_single_block_mapping" if not candidates
                                                 else "ambiguous_whole_insert_chain_mapping")
                    else:
                        qchrom, qstart, qend, strand = candidates[0]
                        record["hg38"] = {"chromosome": qchrom, "start": qstart, "end": qend,
                                          "strand": strand}
                        if not PRIMARY.fullmatch(qchrom):
                            record["reasons"].append("mapped_nonprimary_hg38_contig")
                        elif reference.fetch(qchrom, qstart, qend, strand=strand) != a:
                            record["reasons"].append("full_reference_insert_not_exact_hg38_match")
                        else:
                            counts["exact_hg38_reference_validated_pairs"] += 1
                            # Preserve the complete reference-valid external population.
                            # The frozen transfer protocol forms components on this
                            # population BEFORE applying any prior-exposure guard.
                            full_pair = dict(record, sequence_ref=a, sequence_alt=b,
                                             gc_preserving=(a[99] in "GC") == (b[99] in "GC"))
                            all_mapped.write((json.dumps(full_pair, separators=(",", ":")) + "\n").encode())
                            if hits_guard(panel_guard, qchrom, qstart, qend):
                                record["reasons"].append("within_inclusive_10000bp_of_explored_Agarwal_insert")
                            if hits_guard(other_guard, qchrom, qstart, qend):
                                record["reasons"].append("within_inclusive_10000bp_of_explored_Kircher_or_HBB_locus")
                    if canonical(a) in old_sequences or canonical(b) in old_sequences:
                        record["reasons"].append("exact_or_RC_200nt_match_to_explored_construct")
            record["DNA_gate_passed"] = not record["reasons"]
            for reason in record["reasons"]:
                counts[reason] += 1
            if record["DNA_gate_passed"]:
                counts["DNA_gate_passed_source_class_still_unbound"] += 1
                pair = dict(record, sequence_ref=a, sequence_alt=b,
                            gc_preserving=(a[99] in "GC") == (b[99] in "GC"))
                eligible.write((json.dumps(pair, separators=(",", ":")) + "\n").encode())
                eligible_ids.append(vid)
            ledger.write((json.dumps(record, separators=(",", ":")) + "\n").encode())
    summary = {
        "schema_version": 1,
        "phase": "design_DNA_only_no_new_outcome_payload_acquired",
        "source_class_gate_passed": False,
        "remaining_gates": ["source class excludes technical controls and reused experiments",
                            "full reference-valid population components before whole-component exposure exclusion",
                            "exhaustive old-homology gate", "frozen endpoint/header/value-type binding",
                            "complete prespecified context coverage"],
        "counts": dict(counts),
        "definitions": {"source_coordinates": "hg19, zero-based half-open interval [ID_position_1based-100, ID_position_1based+100)",
                        "mapping": "whole200nt within one chain aligned block; unique coordinate/strand mapping before exact hg38 validation",
                        "old_locus_guard": "inclusive10000bp from every explored Agarwal insert and Kircher target/HBB gene interval",
                        "all_reference_validated_pairs": "all exact hg38 DNA pairs before old-guard/exact-match exclusion; final protocol components must use this file",
                        "eligible_pairs": "rowwise inventory diagnostic only; do not use this prefiltered population to form final components",
                        "not_established": "global genome uniqueness or exhaustive homology; independent experiment status requires source class gate"},
        "eligible_ID_order_sha256": hashlib.sha256("\n".join(eligible_ids).encode()).hexdigest(),
        "inputs": {k: {"path": str(v[0]), "sha256": v[1]} for k, v in bindings.items()},
        "other_explored_inputs": {str(p.relative_to(ROOT)): digest(p) for p in [metadata_path, hbb_path]},
        "implementation_sha256": digest(__file__),
        "files": {p.name: {"sha256": digest(p), "bytes": p.stat().st_size} for p in args.out.iterdir()},
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["counts"], indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--dna", type=Path)
    p.add_argument("--chain", type=Path)
    p.add_argument("--reference", type=Path)
    p.add_argument("--out", type=Path)
    args = p.parse_args()
    if args.self_test:
        self_test()
    else:
        if any(x is None for x in [args.dna, args.chain, args.reference, args.out]):
            p.error("--dna --chain --reference --out are required")
        run(args)


if __name__ == "__main__":
    main()
