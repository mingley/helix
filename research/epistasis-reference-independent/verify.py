#!/usr/bin/env python3
"""Independent, outcome-free reconstruction of the S18 reference ledger.

The XLSX, chain and 2bit readers here use only the Python standard library.
No source outcome workbook is opened. Coordinates are zero-based half-open.
"""
import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import re
import struct
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[2]
SOURCE_SHA = "92929b4cef82797237c4c2be770feb744daa97255ac6dbfe3d53d5a34c996fbb"
CHAIN_SHA = "5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d"
REFERENCE_SHA = "1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
VARIANT = re.compile(r"(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)):([0-9]+):([ACGT]):([ACGT])\Z")
PRIMARY = re.compile(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)\Z")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def require_new_receipt(path):
    """Reject existing files, including symlinks, before doing the audit."""
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"receipt path already exists: {path}")


def write_receipt(path, receipt):
    """Create a receipt exclusively, including if the path appeared after preflight."""
    with path.open("x", encoding="utf-8", newline="\n") as output:
        output.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


def design_rows(path):
    """Read cells directly from the original S18 XLSX ZIP/XML, without openpyxl."""
    with zipfile.ZipFile(path) as z:
        workbook = ET.fromstring(z.read("xl/workbook.xml"))
        sheets = workbook.findall(".//" + NS + "sheet")
        require(len(sheets) == 1 and sheets[0].get("name") == "S18.Variant Pair Oligos",
                "unexpected S18 workbook")
        strings = []
        with z.open("xl/sharedStrings.xml") as raw:
            for _, elem in ET.iterparse(raw, events=("end",)):
                if elem.tag == NS + "si":
                    strings.append("".join(item.text or "" for item in elem.iter(NS + "t")))
                    elem.clear()
        with z.open("xl/worksheets/sheet1.xml") as raw:
            previous = 0
            for _, elem in ET.iterparse(raw, events=("end",)):
                if elem.tag != NS + "row":
                    continue
                n = int(elem.attrib["r"])
                require(n == previous + 1, "nonconsecutive worksheet rows")
                previous = n
                cells = {}
                for cell in elem.findall(NS + "c"):
                    address = re.fullmatch(r"([A-G])([0-9]+)", cell.attrib["r"])
                    require(address is not None and int(address.group(2)) == n,
                            "invalid cell address")
                    require(cell.find(NS + "f") is None, "formula cell in design source")
                    value = cell.find(NS + "v")
                    require(value is not None and value.text is not None, "missing cell value")
                    typ = cell.attrib.get("t", "n")
                    if typ == "s":
                        index = int(value.text)
                        require(0 <= index < len(strings), "invalid shared string index")
                        decoded = strings[index]
                    else:
                        require(typ == "n" and re.fullmatch(r"[0-9]+", value.text) is not None,
                                "unexpected cell type")
                        decoded = int(value.text)
                    col = address.group(1)
                    require(col not in cells, "duplicate cell")
                    cells[col] = decoded
                if n == 1:
                    require(cells == {"A": "Supplementary Table 18 | Oligo sequences for variant pairs"},
                            "unexpected title")
                elif n == 2:
                    require(tuple(cells.get(c) for c in "ABCDEFG") ==
                            ("Variant 1", "Variant 2", "Distance", "Center Variant", "Window",
                             "Diplotype", "Oligo"), "unexpected S18 columns")
                else:
                    require(set(cells) == set("ABCDEFG"), "missing design cell")
                    yield n, tuple(cells[c] for c in "ABCDEFG")
                elem.clear()


def square(key, rows):
    """Independently recognize a square by four-column nucleotide patterns."""
    v1, v2, center, window = key
    a, b = VARIANT.fullmatch(v1), VARIANT.fullmatch(v2)
    if not a or not b or a[1] != b[1] or a[3] == a[4] or b[3] == b[4]:
        return None
    distance = abs(int(a[2]) - int(b[2]))
    if distance < 3 or center not in ("var1", "var2") or window not in ("left", "middle", "right"):
        return None
    if len(rows) != 4 or {r[2] for r in rows} != {"RR", "AR", "RA", "AA"}:
        return None
    if any(r[1] != distance or not isinstance(r[3], str) or len(r[3]) != 200 or
           set(r[3]) - set("ACGT") for r in rows):
        return None
    records = {r[2]: r for r in rows}
    labels = ("RR", "AR", "RA", "AA")
    columns = list(zip(*(records[label][3] for label in labels)))
    changed = [(i, bases) for i, bases in enumerate(columns) if len(set(bases)) > 1]
    first = [i for i, bases in changed if bases == (a[3], a[4], a[3], a[4])]
    second = [i for i, bases in changed if bases == (b[3], b[3], b[4], b[4])]
    if len(changed) != 2 or len(first) != 1 or len(second) != 1 or first[0] == second[0]:
        return None
    if abs(first[0] - second[0]) != distance:
        return None
    return {"variant1": v1, "variant2": v2, "center": center, "window": window,
            "distance": distance, "rows": {label: records[label][0] for label in labels},
            "position1_zero_based": first[0], "position2_zero_based": second[0],
            "reference": records["RR"][3]}


def collect_design(path):
    groups = defaultdict(list)
    row_count = 0
    for n, (v1, v2, distance, center, window, label, seq) in design_rows(path):
        groups[(v1, v2, center, window)].append((n, distance, label, seq))
        row_count += 1
    accepted, rejected = {}, {}
    for key, items in groups.items():
        result = square(key, items)
        if result is None:
            rejected[key] = sorted(item[0] for item in items)
        else:
            accepted[key] = result
    return row_count, accepted, rejected


class ChainIndex:
    def __init__(self, records):
        self.by_chrom = {}
        for chrom, block in records:
            self.by_chrom.setdefault(chrom, []).append(block)
        for chrom, blocks in self.by_chrom.items():
            blocks.sort()
            starts = [b[0] for b in blocks]
            running, ends = -1, []
            for block in blocks:
                running = max(running, block[1])
                ends.append(running)
            self.by_chrom[chrom] = (blocks, starts, ends)

    def map(self, chrom, start, end):
        """Map whole intervals via first and last bases within a single block."""
        require(0 <= start < end, "invalid source interval")
        if chrom not in self.by_chrom:
            return []
        blocks, starts, ends = self.by_chrom[chrom]
        i = bisect_right(starts, start) - 1
        out = set()
        while i >= 0 and ends[i] >= end:
            left, right, dst, size, strand, target = blocks[i]
            if left <= start and end <= right:
                first = target + (start - left)
                last = target + (end - 1 - left)
                if strand == "-":
                    first, last = size - 1 - first, size - 1 - last
                lo, hi = min(first, last), max(first, last) + 1
                require(0 <= lo < hi <= size, "mapped interval out of bounds")
                out.add((dst, lo, hi, strand))
            i -= 1
        return [{"chromosome": c, "start": s, "end": e, "strand": strand}
                for c, s, e, strand in sorted(out)]


def parse_chain(lines):
    records = []
    state = None
    for line in lines:
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "chain":
            require(state is None and len(fields) == 13 and fields[4] == "+" and
                    fields[9] in ("+", "-"), "invalid chain header")
            state = {"src": fields[2], "srcsize": int(fields[3]), "a": int(fields[5]),
                     "aend": int(fields[6]), "dst": fields[7], "dstsize": int(fields[8]),
                     "strand": fields[9], "b": int(fields[10]), "bend": int(fields[11])}
            require(0 <= state["a"] < state["aend"] <= state["srcsize"] and
                    0 <= state["b"] < state["bend"] <= state["dstsize"], "bad chain bounds")
            continue
        require(state is not None and len(fields) in (1, 3), "bad chain block")
        size = int(fields[0])
        require(size > 0 and state["a"] + size <= state["aend"] and
                state["b"] + size <= state["bend"], "bad chain block length")
        records.append((state["src"], (state["a"], state["a"] + size, state["dst"],
                                       state["dstsize"], state["strand"], state["b"])))
        state["a"] += size
        state["b"] += size
        if len(fields) == 3:
            da, db = int(fields[1]), int(fields[2])
            require(da >= 0 and db >= 0, "negative chain gap")
            state["a"] += da
            state["b"] += db
        else:
            require(state["a"] == state["aend"] and state["b"] == state["bend"],
                    "chain endpoint mismatch")
            state = None
    require(state is None, "unterminated chain")
    return ChainIndex(records)


class TwoBit:
    """Independent on-demand twoBit reader; N blocks are restored, masks ignored."""
    def __init__(self, path):
        self.handle = path.open("rb")
        magic = self.handle.read(4)
        self.endian = {b"\x1a\x41\x27\x43": ">", b"\x43\x27\x41\x1a": "<"}.get(magic)
        require(self.endian is not None, "bad 2bit signature")
        version, count, reserved = self.ints(3)
        require(version == 0 and reserved == 0, "unsupported 2bit header")
        self.index = {}
        for _ in range(count):
            length = self.handle.read(1)
            require(len(length) == 1 and length[0] > 0, "invalid 2bit index")
            name = self.handle.read(length[0]).decode("ascii")
            offset, = self.ints(1)
            require(name not in self.index, "duplicate 2bit sequence")
            self.index[name] = offset
        self.records = {}

    def ints(self, count):
        raw = self.handle.read(4 * count)
        require(len(raw) == 4 * count, "truncated 2bit")
        return struct.unpack(self.endian + str(count) + "I", raw)

    def record(self, chrom):
        if chrom not in self.records:
            require(chrom in self.index, "chromosome missing from 2bit")
            self.handle.seek(self.index[chrom])
            size, n_count = self.ints(2)
            starts, lengths = self.ints(n_count), self.ints(n_count)
            mask_count, = self.ints(1)
            self.ints(mask_count)
            self.ints(mask_count)
            reserved, = self.ints(1)
            require(reserved == 0, "bad 2bit record")
            blocks = list(zip(starts, lengths))
            previous = 0
            for start, length in blocks:
                require(length > 0 and previous <= start and start + length <= size,
                        "invalid 2bit N block")
                previous = start + length
            self.records[chrom] = (size, blocks, self.handle.tell())
        return self.records[chrom]

    def fetch(self, chrom, start, end, strand):
        size, nblocks, offset = self.record(chrom)
        require(0 <= start < end <= size and strand in ("+", "-"), "bad 2bit interval")
        self.handle.seek(offset + start // 4)
        count = (end + 3) // 4 - start // 4
        packed = self.handle.read(count)
        require(len(packed) == count, "truncated packed sequence")
        bases = "TCAG"
        decoded = "".join(bases[(byte >> shift) & 3] for byte in packed for shift in (6, 4, 2, 0))
        sequence = list(decoded[start % 4:start % 4 + end - start])
        for left, length in nblocks:
            for i in range(max(left, start), min(left + length, end)):
                sequence[i - start] = "N"
        result = "".join(sequence)
        return result if strand == "+" else result.translate(COMPLEMENT)[::-1]

    def close(self):
        self.handle.close()


def expected_record(item, chains, reference):
    a, b = VARIANT.fullmatch(item["variant1"]), VARIANT.fullmatch(item["variant2"])
    first = int(a[2]) - 1 - item["position1_zero_based"]
    second = int(b[2]) - 1 - item["position2_zero_based"]
    end = first + 200 if first == second and first >= 0 else None
    result = {"variant1": item["variant1"], "variant2": item["variant2"],
              "center": item["center"], "window": item["window"],
              "source_rows": item["rows"], "position1_zero_based": item["position1_zero_based"],
              "position2_zero_based": item["position2_zero_based"],
              "hg19": {"chromosome": a[1], "start_from_variant1": first,
                       "start_from_variant2": second, "end": end},
              "mapping_candidates": [],
              "reference_sequence_sha256": hashlib.sha256(item["reference"].encode("ascii")).hexdigest(),
              "mismatches": None}
    if end is None:
        result["status"] = "inconsistent_start"
        return result
    candidates = chains.map(a[1], first, end)
    result["mapping_candidates"] = candidates
    if len(candidates) == 0:
        result["status"] = "no_map"
    elif len(candidates) > 1:
        result["status"] = "ambiguous_map"
    elif not PRIMARY.fullmatch(candidates[0]["chromosome"]):
        result["status"] = "nonprimary"
    else:
        c = candidates[0]
        genomic = reference.fetch(c["chromosome"], c["start"], c["end"], c["strand"])
        result["mismatches"] = [{"position_zero_based": i, "source": x, "hg38": y}
                                for i, (x, y) in enumerate(zip(item["reference"], genomic)) if x != y]
        result["status"] = "full_reference_mismatch" if result["mismatches"] else "valid"
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=ROOT / "research/epistasis-source-audit/siraj2026-S18.xlsx")
    p.add_argument("--design-report", type=Path, default=ROOT / "research/epistasis-source-audit/s18-design-audit.json")
    p.add_argument("--chain", type=Path, default=ROOT / "research/variant-transfer/source/hg19ToHg38.over.chain.gz")
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--ledger", type=Path, default=ROOT / "research/epistasis-reference-audit/ledger.jsonl")
    p.add_argument("--receipt", type=Path, required=True)
    args = p.parse_args()
    require_new_receipt(args.receipt)
    for path, expected, name in [(args.source, SOURCE_SHA, "S18"), (args.chain, CHAIN_SHA, "chain"),
                                 (args.reference, REFERENCE_SHA, "hg38")]:
        require(sha256(path) == expected, f"{name} input SHA-256 mismatch")
    count, accepted, rejected = collect_design(args.source)
    report = json.loads(args.design_report.read_text())
    published_accepted = {(x["variant1"], x["variant2"], x["center"], x["window"]): x
                          for x in report["eligible_groups"]}
    published_rejected = {tuple(x["key"]): sorted(x["rows"]) for x in report["ineligible_groups"]}
    require(len(published_accepted) == len(report["eligible_groups"]) and
            len(published_rejected) == len(report["ineligible_groups"]), "duplicate design report keys")
    require(count == report["data_rows"] and len(accepted) + len(rejected) == report["design_groups"],
            "design row or group count differs")
    require(rejected == published_rejected, "rejected design identities or row numbers differ")
    require(set(accepted) == set(published_accepted), "accepted design identities differ")
    for key, item in accepted.items():
        without_sequence = {k: v for k, v in item.items() if k != "reference"}
        require(without_sequence == published_accepted[key], f"design details differ: {key}")
    with gzip.open(args.chain, "rt", encoding="ascii") as handle:
        chains = parse_chain(handle)
    reference = TwoBit(args.reference)
    seen, statuses = set(), Counter()
    try:
        with args.ledger.open(encoding="utf-8") as handle:
            for number, line in enumerate(handle, 1):
                actual = json.loads(line)
                key = tuple(actual[k] for k in ("variant1", "variant2", "center", "window"))
                require(key in accepted and key not in seen, f"unexpected or duplicate ledger key at line {number}")
                seen.add(key)
                expected = expected_record(accepted[key], chains, reference)
                if actual != expected:
                    differences = [k for k in set(actual) | set(expected) if actual.get(k) != expected.get(k)]
                    raise ValueError(f"reference ledger line {number} differs: {key}; fields={sorted(differences)}")
                statuses[expected["status"]] += 1
    finally:
        reference.close()
    require(seen == set(accepted), "reference ledger missing accepted design groups")
    receipt = {"method": "Independent XLSX XML genotype patterns, chain block endpoints, and on-demand 2bit decoding",
               "implementation_sha256": sha256(Path(__file__).resolve()),
               "source_sha256": sha256(args.source), "design_report_sha256": sha256(args.design_report),
               "chain_sha256": sha256(args.chain), "reference_sha256": sha256(args.reference),
               "ledger_sha256": sha256(args.ledger), "source_rows": count,
               "accepted_groups": len(accepted), "rejected_groups": len(rejected),
               "ledger_groups_checked": len(seen), "status_counts": dict(sorted(statuses.items())),
               "detailed_rejection_reasons_checked": False, "outcome_values_read": False,
               "all_checks_passed": True}
    write_receipt(args.receipt, receipt)
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
