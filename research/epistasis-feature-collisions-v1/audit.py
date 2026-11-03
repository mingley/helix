#!/usr/bin/env python3
"""Inspect feature collisions within S18 design squares, without outcome values."""
import argparse
from collections import Counter, defaultdict
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "research/epistasis-source-audit/siraj2026-S18.xlsx"
DESIGN = ROOT / "research/epistasis-source-audit/s18-design-audit.json"
HG38 = ROOT / "research/epistasis-reference-audit/ledger.jsonl"
HG19 = ROOT / "research/epistasis-hg19-audit/results/windows.jsonl"
DESIGN_READER = ROOT / "research/epistasis-reference-audit/audit.py"
INPUT_HASHES = {
    "s18_workbook": "92929b4cef82797237c4c2be770feb744daa97255ac6dbfe3d53d5a34c996fbb",
    "s18_design_report": "70dc4a6bd1813ea879b9699565cf75e8e9c21bff9ea1af5a46c25c1ad4686b0f",
    "hg38_reference_ledger": "72800cd3748ecff22fc1cd3d2a9593053353274f91aa9479e47d67dc0b4c55e0",
    "hg19_reference_ledger": "666b62d2259dc03e8a3e9148dadf4d1a503ed18353fdd0bd5997d1e8d2b53265",
    "design_reader": "66413f09afe3ad4ecb71ff8667cf8a95403ea9e4aa2a308585b3f73888ac8c4f",
}
LABELS = ("RR", "AR", "RA", "AA")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


WORDS = sorted({min(word, rc(word)) for word in map("".join, itertools.product("ACGT", repeat=3))})


def features(sequence):
    if len(sequence) != 200 or set(sequence) - set("ACGT"):
        raise ValueError("expected a literal 200-base A/C/G/T sequence")
    counts = Counter(min(sequence[i:i + 3], rc(sequence[i:i + 3])) for i in range(198))
    return tuple(counts[word] for word in WORDS) + (sum(base in "GC" for base in sequence),)


def subtract(first, second):
    return tuple(a - b for a, b in zip(first, second))


def classify(d1, d2):
    """The two separated edit deltas classify all possible square collisions."""
    if len(d1) != 33 or len(d2) != 33:
        raise ValueError("expected 33-dimensional feature deltas")
    zero = (0,) * 33
    if d1 == zero and d2 == zero:
        return "both_blind"
    if d1 == zero:
        return "first_blind_only"
    if d2 == zero:
        return "second_blind_only"
    if d1 == d2:
        return "same_nonzero"
    if d1 == tuple(-x for x in d2):
        return "opposite_nonzero"
    return "no_collision"


def square_features(sequences):
    if set(sequences) != set(LABELS):
        raise ValueError("expected all four construct labels")
    values = {label: features(sequences[label]) for label in LABELS}
    first = subtract(values["AR"], values["RR"])
    second = subtract(values["RA"], values["RR"])
    if subtract(values["AA"], values["RR"]) != tuple(a + b for a, b in zip(first, second)):
        raise ValueError("separated-edit feature parallelogram failed")
    equal = [{"construct1": a, "construct2": b,
              "reverse_complement": sequences[a] == rc(sequences[b])}
             for a, b in itertools.combinations(LABELS, 2) if values[a] == values[b]]
    if any(sequences[p["construct1"]] == sequences[p["construct2"]] for p in equal):
        raise ValueError("two distinct genotype labels have the same DNA sequence")
    kind = classify(first, second)
    expected_pairs = {
        "both_blind": {tuple(x) for x in itertools.combinations(LABELS, 2)},
        "first_blind_only": {("RR", "AR"), ("RA", "AA")},
        "second_blind_only": {("RR", "RA"), ("AR", "AA")},
        "same_nonzero": {("AR", "RA")},
        "opposite_nonzero": {("RR", "AA")},
        "no_collision": set(),
    }
    if {(x["construct1"], x["construct2"]) for x in equal} != expected_pairs[kind]:
        raise ValueError("feature collision classification disagrees with direct comparison")
    return first, second, kind, equal


def key(item):
    return tuple(item[k] for k in ("variant1", "variant2", "center", "window"))


def load_ledger(path):
    result = {}
    with path.open(encoding="utf-8") as f:
        for line in f:
            record = json.loads(line)
            k = key(record)
            if k in result:
                raise ValueError("duplicate reference ledger key")
            result[k] = record
    return result


def checked_inputs():
    paths = {"s18_workbook": SOURCE, "s18_design_report": DESIGN,
             "hg38_reference_ledger": HG38, "hg19_reference_ledger": HG19,
             "design_reader": DESIGN_READER}
    for name, path in paths.items():
        if digest(path) != INPUT_HASHES[name]:
            raise ValueError(f"input hash mismatch: {name}")
    spec = importlib.util.spec_from_file_location("frozen_s18_reference_audit", DESIGN_READER)
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    return reader


def audit(out):
    reader = checked_inputs()
    prior = json.loads(DESIGN.read_text())
    entries = prior["eligible_groups"]
    if len(entries) != 7542 or len(entries) != prior["status_counts_not_mutually_exclusive"]["eligible"]:
        raise ValueError("unexpected previously eligible design population")
    keys = [key(item) for item in entries]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate previously eligible design key")
    sequences = reader.source_sequences(entries, SOURCE)
    hg38, hg19 = load_ledger(HG38), load_ledger(HG19)
    if set(hg38) != set(keys) or set(hg19) != set(keys):
        raise ValueError("reference ledgers do not cover the same design keys")

    out.mkdir(parents=True, exist_ok=False)
    ledger = out / "ledger.jsonl"
    categories, pair_ids, reference_categories = Counter(), defaultdict(set), defaultdict(Counter)
    collision_windows, collisions, rc_collisions = 0, 0, 0
    distinct_dna_pairs = set()
    with ledger.open("w", encoding="utf-8", newline="\n") as f:
        for entry, seqs in zip(entries, sequences):
            chrom, n1, n2, p1, p2, rr = reader.validate_square(entry, seqs)
            if n1 - 1 - p1 != n2 - 1 - p2 or abs(p1 - p2) < 3:
                raise ValueError("design violates separated signed-coordinate contract")
            k = key(entry)
            r38, r19 = hg38[k], hg19[k]
            rr_hash = hashlib.sha256(rr.encode("ascii")).hexdigest()
            if (r38["source_rows"] != entry["rows"] or
                    r38["position1_zero_based"] != p1 or r38["position2_zero_based"] != p2 or
                    r38["reference_sequence_sha256"] != rr_hash or
                    r19["source_rr_row_1based"] != entry["rows"]["RR"] or
                    r19["source_rr_sha256"] != rr_hash or
                    r19["hg19_chromosome"] != chrom or r19["hg19_start"] != n1 - 1 - p1 or
                    not isinstance(r19["exact_match"], bool)):
                raise ValueError("reference ledger identity or RR source binding disagrees")
            first, second, category, equal = square_features(seqs)
            record = {"variant1": entry["variant1"], "variant2": entry["variant2"],
                      "center": entry["center"], "window": entry["window"],
                      "source_rows": entry["rows"], "position1_zero_based": p1,
                      "position2_zero_based": p2, "hg38_status": r38["status"],
                      "hg19_exact_match": r19["exact_match"],
                      "first_feature_delta": first, "second_feature_delta": second,
                      "classification": category, "feature_equal_pairs": equal}
            f.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
            categories[category] += 1
            pair_ids[category].add(k[:2])
            reference_categories[r38["status"]][category] += 1
            if equal:
                collision_windows += 1
            for pair in equal:
                a, b = seqs[pair["construct1"]], seqs[pair["construct2"]]
                distinct_dna_pairs.add(tuple(sorted((a, b))))
                collisions += 1
                rc_collisions += pair["reverse_complement"]
    summary = {
        "schema_version": 1,
        "scope": "S18 design sequences only; no outcomes or fitted coefficients read",
        "feature_order": ["rc3mer_" + word for word in WORDS] + ["gc_count"],
        "features": "32 canonical overlapping 3mer integer counts followed by GC count, on 200-base DNA",
        "first_delta": "AR minus RR", "second_delta": "RA minus RR",
        "rules": "compare all four complete feature vectors; require AA - RR = first delta + second delta",
        "design_windows": len(entries),
        "distinct_pair_identifiers": len({k[:2] for k in keys}),
        "classification_windows": dict(sorted(categories.items())),
        "classification_pair_identifiers": {k: len(v) for k, v in sorted(pair_ids.items())},
        "classification_windows_by_hg38_status": {k: dict(sorted(v.items())) for k, v in sorted(reference_categories.items())},
        "windows_with_any_equal_feature_pair": collision_windows,
        "pair_identifiers_with_any_equal_feature_pair": len(set().union(*(v for k, v in pair_ids.items() if k != "no_collision"))),
        "equal_feature_construct_pairs_over_windows": collisions,
        "distinct_unordered_exact_dna_pairs": len(distinct_dna_pairs),
        "equal_feature_pairs_that_are_global_reverse_complements": rc_collisions,
        "inputs_sha256": INPUT_HASHES,
    }
    summary_path = out / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {"implementation_sha256": digest(Path(__file__)),
                "ledger_sha256": digest(ledger), "summary_sha256": digest(summary_path)}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("classification_windows", "windows_with_any_equal_feature_pair",
                      "pair_identifiers_with_any_equal_feature_pair", "equal_feature_construct_pairs_over_windows",
                      "distinct_unordered_exact_dna_pairs", "equal_feature_pairs_that_are_global_reverse_complements")},
                     indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="new output directory")
    args = parser.parse_args()
    audit(args.out)


if __name__ == "__main__":
    main()
