#!/usr/bin/env python3
"""Independently check S18 four-construct feature collisions.

The original S18 design XLSX is decoded directly as ZIP/XML. Each feature vector
is computed from the complete 200-base insert, without importing the primary
collision code. No workbook containing response outcomes is opened.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "research/epistasis-source-audit/siraj2026-S18.xlsx"
PRIOR = ROOT / "research/epistasis-source-audit/s18-design-audit.json"
HG38 = ROOT / "research/epistasis-reference-audit/ledger.jsonl"
HG19 = ROOT / "research/epistasis-hg19-audit/results/windows.jsonl"
PRIMARY = ROOT / "research/epistasis-feature-collisions-v1/results"
HASHES = {
    "source": "92929b4cef82797237c4c2be770feb744daa97255ac6dbfe3d53d5a34c996fbb",
    "prior": "70dc4a6bd1813ea879b9699565cf75e8e9c21bff9ea1af5a46c25c1ad4686b0f",
    "hg38": "72800cd3748ecff22fc1cd3d2a9593053353274f91aa9479e47d67dc0b4c55e0",
    "hg19": "666b62d2259dc03e8a3e9148dadf4d1a503ed18353fdd0bd5997d1e8d2b53265",
}
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
LABELS = ("RR", "AR", "RA", "AA")
VARIANT = re.compile(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y):[0-9]+:[ACGT]:[ACGT]\Z")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
PAIRS = tuple(combinations(LABELS, 2))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


WORDS = ["".join(x) for x in product("ACGT", repeat=3)]
CANONICAL = sorted({min(w, reverse_complement(w)) for w in WORDS})
require(len(CANONICAL) == 32, "unexpected canonical 3-mer classes")


def feature(seq):
    """Count full insert 3-mers by RC class and GC; no floating point."""
    require(len(seq) == 200 and not (set(seq) - set("ACGT")), "not a literal 200-base insert")
    # Generate classes from all 64 words; count both orientations separately.
    counts = Counter(seq[i:i + 3] for i in range(198))
    features = [counts[word] + counts[reverse_complement(word)] for word in CANONICAL]
    require(sum(features) == 198, "3-mer counts do not cover full insert")
    return features + [seq.count("C") + seq.count("G")]


def classify_equalities(equal):
    """Classify equality patterns directly, independent of delta-vector tests."""
    equal = frozenset(equal)
    patterns = {
        frozenset(PAIRS): "both_blind",
        frozenset((("RR", "AR"), ("RA", "AA"))): "first_blind_only",
        frozenset((("RR", "RA"), ("AR", "AA"))): "second_blind_only",
        frozenset((("AR", "RA"),)): "same_nonzero",
        frozenset((("RR", "AA"),)): "opposite_nonzero",
        frozenset(): "no_collision",
    }
    require(equal in patterns, f"unsupported feature equality pattern: {sorted(equal)}")
    return patterns[equal]


def expected_features(seqs):
    require(set(seqs) == set(LABELS), "expected four different diplotype labels")
    vectors = {label: feature(seqs[label]) for label in LABELS}
    equal = [(a, b) for a, b in PAIRS if vectors[a] == vectors[b]]
    delta1 = [a - b for a, b in zip(vectors["AR"], vectors["RR"])]
    delta2 = [a - b for a, b in zip(vectors["RA"], vectors["RR"])]
    # This is independently verified on every whole-sequence vector; it is
    # required for the mutually exclusive square classification above.
    require(all(x + y == u + v for x, y, u, v in
                zip(vectors["AA"], vectors["RR"], vectors["AR"], vectors["RA"])),
            "distant-edit feature additivity failed")
    return {
        "first_feature_delta": delta1,
        "second_feature_delta": delta2,
        "classification": classify_equalities(equal),
        "feature_equal_pairs": [
            {"construct1": a, "construct2": b,
             "reverse_complement": seqs[a] == reverse_complement(seqs[b])}
            for a, b in equal
        ],
    }


def xml_rows(source):
    """Read the S18 worksheet directly from OpenXML, including all S18 rows."""
    require(sha(source) == HASHES["source"], "unexpected S18 source hash")
    with zipfile.ZipFile(source) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheets = workbook.findall(".//" + NS + "sheet")
        require(len(sheets) == 1 and sheets[0].get("name") == "S18.Variant Pair Oligos",
                "unexpected workbook sheets")
        strings = []
        with archive.open("xl/sharedStrings.xml") as raw:
            for _, elem in ET.iterparse(raw, events=("end",)):
                if elem.tag == NS + "si":
                    strings.append("".join(x.text or "" for x in elem.iter(NS + "t")))
                    elem.clear()
        last_row = 0
        with archive.open("xl/worksheets/sheet1.xml") as raw:
            for _, elem in ET.iterparse(raw, events=("end",)):
                if elem.tag != NS + "row":
                    continue
                number = int(elem.attrib["r"])
                require(number == last_row + 1, "nonconsecutive source rows")
                last_row = number
                cells = {}
                for cell in elem.findall(NS + "c"):
                    addr = re.fullmatch(r"([A-G])([0-9]+)", cell.attrib["r"])
                    require(addr is not None and int(addr.group(2)) == number, "bad cell address")
                    require(cell.find(NS + "f") is None, "formula in source cell")
                    value = cell.find(NS + "v")
                    require(value is not None and value.text is not None, "missing source value")
                    typ = cell.attrib.get("t", "n")
                    if typ == "s":
                        index = int(value.text)
                        require(0 <= index < len(strings), "bad shared string index")
                        decoded = strings[index]
                    else:
                        require(typ == "n" and re.fullmatch(r"[0-9]+", value.text) is not None,
                                "unexpected source cell type")
                        decoded = int(value.text)
                    require(addr.group(1) not in cells, "duplicate source cell")
                    cells[addr.group(1)] = decoded
                if number == 1:
                    require(cells == {"A": "Supplementary Table 18 | Oligo sequences for variant pairs"},
                            "unexpected worksheet title")
                elif number == 2:
                    require(tuple(cells.get(c) for c in "ABCDEFG") ==
                            ("Variant 1", "Variant 2", "Distance", "Center Variant", "Window",
                             "Diplotype", "Oligo"), "unexpected worksheet columns")
                else:
                    require(set(cells) == set("ABCDEFG"), "missing or extra source columns")
                    yield number, tuple(cells[c] for c in "ABCDEFG")
                elem.clear()


def square(key, rows):
    """Reconstruct the prior eligibility from the four complete sequences."""
    v1, v2, center, window = key
    if not all(isinstance(x, str) and VARIANT.fullmatch(x) for x in (v1, v2)):
        return None
    c1, n1, r1, a1 = v1.split(":")
    c2, n2, r2, a2 = v2.split(":")
    distance = abs(int(n1) - int(n2))
    if (c1 != c2 or distance < 3 or r1 == a1 or r2 == a2 or
            center not in ("var1", "var2") or window not in ("left", "middle", "right")):
        return None
    if len(rows) != 4 or {x[2] for x in rows} != set(LABELS):
        return None
    if any(d != distance or not isinstance(s, str) or len(s) != 200 or set(s) - set("ACGT")
           for _, d, _, s in rows):
        return None
    rowbylabel = {item[2]: item for item in rows}
    seqs = {label: rowbylabel[label][3] for label in LABELS}
    columns = list(zip(*(seqs[label] for label in LABELS)))
    varied = [(i, bases) for i, bases in enumerate(columns) if len(set(bases)) != 1]
    first = [i for i, bases in varied if bases == (r1, a1, r1, a1)]
    second = [i for i, bases in varied if bases == (r2, r2, a2, a2)]
    if len(varied) != 2 or len(first) != 1 or len(second) != 1:
        return None
    if first[0] == second[0] or abs(first[0] - second[0]) != distance:
        return None
    prior = {"variant1": v1, "variant2": v2, "center": center, "window": window,
             "distance": distance, "rows": {k: rowbylabel[k][0] for k in LABELS},
             "position1_zero_based": first[0], "position2_zero_based": second[0]}
    return prior, seqs


def collect_design(source):
    groups = defaultdict(list)
    count = 0
    for number, (v1, v2, distance, center, window, label, seq) in xml_rows(source):
        groups[(v1, v2, center, window)].append((number, distance, label, seq))
        count += 1
    accepted, rejected = {}, {}
    for key, rows in groups.items():
        result = square(key, rows)
        if result is None:
            rejected[key] = sorted(x[0] for x in rows)
        else:
            accepted[key] = result
    return count, accepted, rejected


def key_of(record):
    return tuple(record[x] for x in ("variant1", "variant2", "center", "window"))


def reference_status(path, hash_key, attr):
    require(sha(path) == HASHES[hash_key], f"{hash_key} input SHA-256 mismatch")
    result = {}
    with path.open(encoding="utf-8") as f:
        for line in f:
            record = json.loads(line)
            key = key_of(record)
            require(key not in result, f"duplicate key in {hash_key} reference audit")
            result[key] = record[attr]
    return result


def verify(args):
    for path, key in ((args.source, "source"), (args.prior, "prior")):
        require(sha(path) == HASHES[key], f"{key} input SHA-256 mismatch")
    count, accepted, rejected = collect_design(args.source)
    prior = json.loads(args.prior.read_text())
    expected_accepted = {key_of(x): x for x in prior["eligible_groups"]}
    expected_rejected = {tuple(x["key"]): sorted(x["rows"]) for x in prior["ineligible_groups"]}
    require(len(expected_accepted) == len(prior["eligible_groups"]) and
            len(expected_rejected) == len(prior["ineligible_groups"]), "duplicate prior design keys")
    require(count == prior["data_rows"] and len(accepted) + len(rejected) == prior["design_groups"],
            "S18 group or row totals differ")
    require(rejected == expected_rejected, "rejected design keys or rows differ")
    require({key: result[0] for key, result in accepted.items()} == expected_accepted,
            "accepted design details differ")
    hg38 = reference_status(args.hg38, "hg38", "status")
    hg19 = reference_status(args.hg19, "hg19", "exact_match")
    require(set(hg38) == set(accepted) == set(hg19), "reference audit group keys differ")
    seen = set()
    classifications, pair_count, rc_pairs = Counter(), 0, 0
    categories_by_hg38 = defaultdict(Counter)
    category_pair_ids = defaultdict(set)
    colliding_pair_ids, distinct_dna_pairs = set(), set()
    with args.ledger.open(encoding="utf-8") as f:
        for number, line in enumerate(f, 1):
            actual = json.loads(line)
            key = key_of(actual)
            require(key in accepted and key not in seen, f"unexpected or duplicate ledger key at {number}")
            seen.add(key)
            prior_record, seqs = accepted[key]
            expected = {k: prior_record[k] for k in ("variant1", "variant2", "center", "window",
                                                       "position1_zero_based", "position2_zero_based")}
            expected.update({"source_rows": prior_record["rows"], "hg38_status": hg38[key],
                             "hg19_exact_match": hg19[key]})
            expected.update(expected_features(seqs))
            if actual != expected:
                different = [k for k in set(actual) | set(expected) if actual.get(k) != expected.get(k)]
                raise ValueError(f"feature ledger line {number} differs: {key}; fields={sorted(different)}")
            classifications[expected["classification"]] += 1
            categories_by_hg38[hg38[key]][expected["classification"]] += 1
            category_pair_ids[expected["classification"]].add(key[:2])
            pair_count += len(expected["feature_equal_pairs"])
            if expected["feature_equal_pairs"]:
                colliding_pair_ids.add(key[:2])
            for item in expected["feature_equal_pairs"]:
                a, b = seqs[item["construct1"]], seqs[item["construct2"]]
                require(a != b, "two different diplotypes contain identical DNA")
                distinct_dna_pairs.add(tuple(sorted((a, b))))
                rc_pairs += item["reverse_complement"]
    require(seen == set(accepted), "primary feature ledger does not cover every eligible group")
    summary = json.loads(args.summary.read_text())
    expected_summary = {
        "feature_order": ["rc3mer_" + word for word in CANONICAL] + ["gc_count"],
        "design_windows": len(seen),
        "distinct_pair_identifiers": len({key[:2] for key in seen}),
        "classification_windows": dict(sorted(classifications.items())),
        "classification_pair_identifiers": {k: len(v) for k, v in sorted(category_pair_ids.items())},
        "classification_windows_by_hg38_status": {k: dict(sorted(v.items()))
                                                    for k, v in sorted(categories_by_hg38.items())},
        "windows_with_any_equal_feature_pair": sum(v for k, v in classifications.items()
                                                    if k != "no_collision"),
        "pair_identifiers_with_any_equal_feature_pair": len(colliding_pair_ids),
        "equal_feature_construct_pairs_over_windows": pair_count,
        "distinct_unordered_exact_dna_pairs": len(distinct_dna_pairs),
        "equal_feature_pairs_that_are_global_reverse_complements": rc_pairs,
    }
    for name, value in expected_summary.items():
        require(summary[name] == value, f"primary summary differs for {name}")
    manifest = json.loads(args.manifest.read_text())
    require(manifest["ledger_sha256"] == sha(args.ledger) and
            manifest["summary_sha256"] == sha(args.summary), "primary manifest output hashes disagree")
    require(manifest["implementation_sha256"] == sha(args.implementation),
            "primary manifest implementation hash disagrees")
    return {
        "source_sha256": sha(args.source), "prior_design_sha256": sha(args.prior),
        "hg38_reference_audit_sha256": sha(args.hg38), "hg19_reference_audit_sha256": sha(args.hg19),
        "primary_ledger_sha256": sha(args.ledger), "primary_summary_sha256": sha(args.summary),
        "primary_manifest_sha256": sha(args.manifest), "implementation_sha256": sha(Path(__file__).resolve()),
        "source_rows": count, "accepted_design_groups": len(accepted), "rejected_design_groups": len(rejected),
        "full_ledger_records_compared": len(seen), "classes": dict(sorted(classifications.items())),
        "windows_with_feature_collisions": expected_summary["windows_with_any_equal_feature_pair"],
        "pair_identifiers_with_feature_collisions": len(colliding_pair_ids),
        "feature_equal_construct_pairs": pair_count, "distinct_dna_pairs": len(distinct_dna_pairs),
        "global_reverse_complement_pairs": rc_pairs,
        "outcome_workbooks_read": False,
        "method": "Independent XML design reconstruction; complete-sequence RC 3-mer and GC counting",
        "all_checks_passed": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--prior", type=Path, default=PRIOR)
    parser.add_argument("--hg38", type=Path, default=HG38)
    parser.add_argument("--hg19", type=Path, default=HG19)
    parser.add_argument("--ledger", type=Path, default=PRIMARY / "ledger.jsonl")
    parser.add_argument("--summary", type=Path, default=PRIMARY / "summary.json")
    parser.add_argument("--manifest", type=Path, default=PRIMARY / "manifest.json")
    parser.add_argument("--implementation", type=Path, default=PRIMARY.parent / "audit.py")
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    receipt = verify(args)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
