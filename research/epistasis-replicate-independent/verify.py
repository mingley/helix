#!/usr/bin/env python3
"""Independently reproduce the exploratory OL27 name-overlap inventory.

Candidate names are only syntactic hypotheses, not proven insert identities.
The gzipped barcode lookup contains mapping/QC fields; they are instantiated as
text by the TSV parser, but only name and overall flag columns are interpreted.
No RNA, reporter-activity or outcome table is read.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "research/epistasis-source-audit/s18-design-audit.json"
FEATURES = ROOT / "research/epistasis-feature-collisions-v1/results/ledger.jsonl"
PRIMARY = ROOT / "research/epistasis-replicate-metadata/ol27-name-overlap.json"
PRIMARY_IMPLEMENTATION = ROOT / "research/epistasis-replicate-metadata/inventory.py"
INPUT_SHA = {
    "design": "70dc4a6bd1813ea879b9699565cf75e8e9c21bff9ea1af5a46c25c1ad4686b0f",
    "feature_ledger": "53c2751cd093da20a2c1bad0cd602e8ed173b7f164456c60f9938a7e1efa9ce7",
    "lookup_gzip": "58f0c48c6d0a21b4a558794f31dcf38ccd6e8c766d56450512a28a27de059a15",
}
RAW_MD5 = "c18ba6775c64fbef1da73cede0d29dea"
RAW_BYTES = 1039687724
LABELS = ("RR", "AR", "RA", "AA")
ALLOWED = {"no_collision", "both_blind", "first_blind_only", "second_blind_only",
           "same_nonzero", "opposite_nonzero"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def key(group):
    return tuple(group[k] for k in ("variant1", "variant2", "center", "window"))


def syntactic_names(group):
    """Construct four tentative names by selecting centered and other indices."""
    ids = [group["variant1"], group["variant2"]]
    variants = []
    for value in ids:
        parts = value.split(":")
        require(len(parts) == 4 and parts[0].startswith("chr") and len(parts[0]) > 3,
                "invalid source variant identifier")
        variants.append((parts[0][3:], parts[1], parts[2], parts[3]))
    center_names = {"var1": 0, "var2": 1}
    window_names = {"left": "L", "middle": "C", "right": "R"}
    require(group["center"] in center_names and group["window"] in window_names,
            "invalid source center or window")
    centered = center_names[group["center"]]
    other = 1 - centered
    chrom, pos, ref, alt = variants[centered]
    other_chrom, other_pos, _, _ = variants[other]
    code = window_names[group["window"]]
    result = []
    for label in LABELS:
        result.append(f"{chrom}:{pos}:{ref}:{alt}:{label[centered]}:w{code}:"
                      f"Alt_{other_chrom}-{other_pos}-{label[other]}")
    require(len(set(result)) == 4, "construct candidate names not distinct within window")
    return tuple(result)


def fields(row):
    """Interpret only the oligo-name list and overall flag in a 10-column row."""
    require(len(row) == 10, "lookup record is not ten-column TSV")
    require(row[4] in ("0", "1", "2"), "unexpected mapping flag")
    names = row[1].split(",")
    require(names and all(names), "empty candidate oligo name")
    return names, row[4]


class HashingStream(io.RawIOBase):
    """Hash uncompressed bytes while the separate CSV reader consumes them."""
    def __init__(self, source):
        super().__init__()
        self.source = source
        self.digest = hashlib.md5()
        self.bytes_read = 0

    def readable(self):
        return True

    def readinto(self, output):
        data = self.source.read(len(output))
        size = len(data)
        output[:size] = data
        self.digest.update(data)
        self.bytes_read += size
        return size


def scan_lookup(path, candidates):
    seen = set()
    strict_seen = set()
    totals = Counter()
    with gzip.open(path, "rb") as compressed:
        hashing = HashingStream(compressed)
        with io.TextIOWrapper(io.BufferedReader(hashing, 1024 * 1024), encoding="utf-8", newline="") as text:
            for row in csv.reader(text, delimiter="\t", quoting=csv.QUOTE_NONE, strict=True):
                names, flag = fields(row)
                totals["rows"] += 1
                totals[f"flag_{flag}_rows"] += 1
                if len(names) > 1:
                    totals["multiple_name_rows"] += 1
                for name in names:
                    if name in candidates:
                        seen.add(name)
                        if flag == "0" and len(names) == 1:
                            strict_seen.add(name)
    return seen, strict_seen, totals, hashing.digest.hexdigest(), hashing.bytes_read


def counts(groups, classifications, seen):
    overall, colliding, collision_classes = Counter(), Counter(), Counter()
    complete_pair_ids, collision_pair_ids = set(), set()
    for group in groups:
        num_seen = sum(candidate in seen for candidate in syntactic_names(group))
        overall[num_seen] += 1
        pair_id = (group["variant1"], group["variant2"])
        classification = classifications[key(group)]
        if num_seen == 4:
            complete_pair_ids.add(pair_id)
        if classification != "no_collision":
            colliding[num_seen] += 1
            if num_seen == 4:
                collision_pair_ids.add(pair_id)
                collision_classes[classification] += 1
    return {
        "windows_by_number_of_four_candidate_strings_present": {str(i): overall[i] for i in range(5)},
        "pair_ids_with_at_least_one_four_string_window": len(complete_pair_ids),
        "collision_windows_by_number_of_four_candidate_strings_present": {str(i): colliding[i] for i in range(5)},
        "collision_pair_ids_with_four_string_window": len(collision_pair_ids),
        "collision_classifications_for_four_string_windows": dict(sorted(collision_classes.items())),
    }


def numeric_paths(value, prefix=()):
    """Collect all numeric report fields, failing on unexpected new fields."""
    if type(value) in (int, float):
        return {prefix: value}
    if isinstance(value, dict):
        result = {}
        for k, v in value.items():
            result.update(numeric_paths(v, prefix + (str(k),)))
        return result
    if isinstance(value, list):
        result = {}
        for i, v in enumerate(value):
            result.update(numeric_paths(v, prefix + (str(i),)))
        return result
    return {}


def run(args):
    for path, name in ((args.design, "design"), (args.features, "feature_ledger"),
                       (args.lookup, "lookup_gzip")):
        require(sha(path) == INPUT_SHA[name], f"{name} input SHA-256 mismatch")
    design = json.loads(args.design.read_text())
    groups = design["eligible_groups"]
    group_keys = [key(g) for g in groups]
    require(len(set(group_keys)) == len(group_keys), "duplicate design group")
    classifications = {}
    with args.features.open(encoding="utf-8") as stream:
        for line in stream:
            entry = json.loads(line)
            k = key(entry)
            require(k not in classifications and entry["classification"] in ALLOWED,
                    "duplicate or unknown feature classification")
            classifications[k] = entry["classification"]
    require(set(group_keys) == set(classifications), "design and feature groups differ")
    candidates = defaultdict(set)
    for group in groups:
        for label, name in zip(LABELS, syntactic_names(group)):
            candidates[name].add((key(group), label))
    seen, strict_seen, totals, raw_md5, raw_bytes = scan_lookup(args.lookup, candidates)
    require(raw_md5 == RAW_MD5 and raw_bytes == RAW_BYTES, "lookup uncompressed checksum or size mismatch")
    expected_numeric = {
        "schema_version": 1,
        "design_windows": len(groups),
        "distinct_candidate_strings": len(candidates),
        "candidate_strings_reused_by_multiple_window_labels": sum(len(x) > 1 for x in candidates.values()),
        "lookup_uncompressed_bytes": raw_bytes,
        "lookup_row_counts": dict(sorted(totals.items())),
        "seen_under_any_flag": counts(groups, classifications, seen),
        "seen_on_a_flag_0_single_name_row": counts(groups, classifications, strict_seen),
    }
    primary = json.loads(args.primary.read_text())
    actual_fields, expected_fields = numeric_paths(primary), numeric_paths(expected_numeric)
    require(actual_fields == expected_fields,
            f"primary numerical fields differ: {sorted(set(actual_fields) ^ set(expected_fields)) or 'values'}")
    require(primary["design_sha256"] == INPUT_SHA["design"] and
            primary["feature_ledger_sha256"] == INPUT_SHA["feature_ledger"] and
            primary["lookup_gzip_sha256"] == INPUT_SHA["lookup_gzip"] and
            primary["lookup_uncompressed_md5"] == raw_md5 and
            primary["implementation_sha256"] == sha(args.primary_implementation),
            "primary provenance fields disagree")
    return {
        "method": "Independent csv TSV parser, streaming decompressed checksum, and centered-index name construction",
        "design_sha256": INPUT_SHA["design"], "feature_ledger_sha256": INPUT_SHA["feature_ledger"],
        "lookup_gzip_sha256": INPUT_SHA["lookup_gzip"], "lookup_uncompressed_md5": raw_md5,
        "primary_report_sha256": sha(args.primary), "primary_implementation_sha256": sha(args.primary_implementation),
        "implementation_sha256": sha(Path(__file__).resolve()),
        "numerical_fields_compared": len(expected_fields), "numerical_summary": expected_numeric,
        "source_scope": "Lookup names are syntactic candidates; mapping/QC columns instantiated but not interpreted except overall flag",
        "rna_or_reporter_outcome_table_read": False, "all_checks_passed": True,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--lookup", type=Path, required=True)
    p.add_argument("--design", type=Path, default=DESIGN)
    p.add_argument("--features", type=Path, default=FEATURES)
    p.add_argument("--primary", type=Path, default=PRIMARY)
    p.add_argument("--primary-implementation", type=Path, default=PRIMARY_IMPLEMENTATION)
    p.add_argument("--receipt", type=Path, required=True)
    args = p.parse_args()
    require(not args.receipt.exists(), "receipt already exists")
    receipt = run(args)
    with args.receipt.open("x", encoding="utf-8", newline="\n") as out:
        out.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
