#!/usr/bin/env python3
"""Post-hoc full-panel feature collisions and v3 observed-label loss floors."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "datasets/context-panel/panel.json"
FEATURES = ROOT / "research/results/context-activity-v3/core/features.json.gz"
ACTIVITY = ROOT / "research/results/context-activity-v3/core/activity.json.gz"
GROUPING = ROOT / "research/results/context-activity-v3/core/grouping.json"
INPUTS = ROOT / "research/results/context-activity-v3/core/inputs.json"
PUBLICATION = ROOT / "research/results/context-activity-v3/publication.json"
MANIFEST = ROOT / "datasets/context-activity/manifest.json"
PROTOCOL = ROOT / "experiments/context-activity-v3.json"
GROUPING_PROTOCOL = ROOT / "experiments/context-specificity-v2.json"
SOURCE = ROOT / "datasets/context-endpoint-recovery/raw/jointdata.xlsx"
EXPECTED = {
    "panel": "fce6b2a1aaf5aa4d530a2e298dccb840bc7c594a05b6831fdb0c6c0f9e2a9be3",
    "features": "1a9b093ca7d0ae56b19294a53c427f0d7d00ee4821ff09e03d8384d601b9417a",
    "activity": "83b2df30128f3ea8326d6d5d1de73b36e32ced113cccebebc915d4fc06e8054b",
    "grouping": "333069c4ab41ed8ea6976b0bc50cd08de9a33118ab60d683533b7ee32a6f517f",
    "inputs": "fe78220d9f94b5022582e79e02ed8d1dbd927c8d292825655ac3d018d04ef6d1",
    "publication": "39cd6a1902ed3b7a4709294b293590309af5f38d2d8c578f82922da74a63a335",
    "source_manifest": "6a3371bed4061f9ca9fa1a95bfc0947ffee91daa22fe509063708e549d453730",
    "protocol": "64e14037fb1f42ddc1ecaca1d23bb8096fabac1f1e892d105eaa20f525fbbb20",
    "grouping_protocol": "8f0576845fecb5cfcf9bfa354fab446d849244448a3489905a6ce33bdeb21cec",
    "raw_workbook": "9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e",
}
CONTEXTS = ("HepG2", "K562", "WTC11")
ENDPOINT = "published_joint_log2_rna_dna_activity"
DNA_TRANSLATE = str.maketrans("ACGT", "TGCA")
JSON_STRING = r'"(?:[^"\\]|\\.)*"'
PANEL_PREFIX = re.compile(
    r'^\{"id":(' + JSON_STRING + r'),"sequence":"([ACGT]+)","chromosome":('
    + JSON_STRING + r'),"start":(\d+),"end":(\d+),"parent_id":(' + JSON_STRING + r')'
)
PANEL_SOURCE = re.compile(
    r'"reference_match_strand":"([+-])","source_declared_strand":"([+-])",'
    r'"source_design_metadata_record":' + JSON_STRING + r',"core_sha256":"([0-9a-f]{64})"'
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def reverse_complement(sequence):
    return sequence.translate(DNA_TRANSLATE)[::-1]


ALL_WORDS = tuple(map("".join, itertools.product("ACGT", repeat=3)))
CANONICAL = {word: min(word, reverse_complement(word)) for word in ALL_WORDS}
WORDS = tuple(sorted(set(CANONICAL.values())))
FEATURE_NAMES = ["rc3mer_" + word for word in WORDS] + ["gc_fraction"]


def features(sequence):
    require(isinstance(sequence, str) and len(sequence) == 200 and not set(sequence) - set("ACGT"),
            "sequence is not literal 200-base A/C/G/T")
    counts = Counter(CANONICAL[sequence[i:i + 3]] for i in range(198))
    return tuple(counts[word] for word in WORDS) + (sum(base in "GC" for base in sequence),)


def feature_fingerprint(counts):
    return hashlib.sha256(json.dumps(list(counts), separators=(",", ":")).encode("ascii")).hexdigest()


def rational(value):
    # JSON consumers commonly decode numbers as binary64; decimal strings keep
    # arbitrarily large exact rational integers portable across those readers.
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def squared_error_sum(values):
    require(bool(values), "empty feature class")
    mean = sum(values, Fraction()) / len(values)
    return sum(((v - mean) ** 2 for v in values), Fraction())


def records_by_id(records, what):
    result = {}
    for record in records:
        identifier = record["id"]
        require(isinstance(identifier, str) and identifier and identifier not in result,
                f"invalid or duplicate {what} ID")
        result[identifier] = record
    return result


def read_json(path):
    with (gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz"
          else path.open(encoding="utf-8")) as f:
        return json.load(f, parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def read_panel():
    """Select only ordered DNA and metadata fields; do not interpret old ESS labels."""
    result = {}
    with PANEL.open(encoding="utf-8") as f:
        header = f.readline().rstrip("\n")
        require(header == '{"schema_version":1,"endpoint":"published_element_specificity_score",'
                          '"contexts":["HepG2","K562","WTC11"],"records":[',
                "unexpected panel header")
        for line in f:
            if line == "]}\n" or line == "]}":
                break
            match, source = PANEL_PREFIX.match(line), PANEL_SOURCE.search(line)
            require(match is not None and source is not None, "malformed panel record")
            identifier, sequence, chromosome, start, end, parent = match.groups()
            identifier, chromosome, parent = map(json.loads, (identifier, chromosome, parent))
            start, end = int(start), int(end)
            require(identifier not in result and end - start == 200, "duplicate or invalid panel record")
            require(len(sequence) == 200 and hashlib.sha256(sequence.encode("ascii")).hexdigest() == source[3],
                    "panel source core hash or length differs")
            result[identifier] = {"id": identifier, "sequence": sequence, "chromosome": chromosome,
                                  "start": start, "end": end, "parent_id": parent,
                                  "reference_match_strand": source[1], "source_declared_strand": source[2]}
        else:
            raise ValueError("panel closing delimiter missing")
        require(not f.read().strip(), "unexpected trailing panel content")
    return result


def provenance():
    paths = {"panel": PANEL, "features": FEATURES, "activity": ACTIVITY, "grouping": GROUPING,
             "inputs": INPUTS, "publication": PUBLICATION, "source_manifest": MANIFEST,
             "protocol": PROTOCOL, "grouping_protocol": GROUPING_PROTOCOL, "raw_workbook": SOURCE}
    for name, path in paths.items():
        require(sha(path) == EXPECTED[name], f"input SHA-256 mismatch: {name}")
    inputs, publication, manifest, protocol = map(read_json, (INPUTS, PUBLICATION, MANIFEST, PROTOCOL))
    require(inputs["activity_sha256"] == EXPECTED["activity"] and
            inputs["dna_panel_sha256"] == EXPECTED["panel"] and
            inputs["grouping_sha256"] == EXPECTED["grouping"] and
            inputs["protocol_sha256"] == EXPECTED["protocol"] and
            inputs["raw_source_sha256"] == EXPECTED["raw_workbook"] and
            inputs["source_manifest_sha256"] == EXPECTED["source_manifest"] and
            inputs["endpoint"] == ENDPOINT, "v3 input provenance disagrees")
    require(publication["published_activity_sha256"] == EXPECTED["activity"] and
            publication["protocol_sha256"] == EXPECTED["protocol"] and
            publication["source"]["sha256"] == EXPECTED["raw_workbook"] and
            publication["source_manifest_sha256"] == EXPECTED["source_manifest"] and
            publication["immutable_previous_protocols"]["context_specificity_v2"] == EXPECTED["grouping_protocol"] and
            publication["literal_activity_binding_verified"] is True and
            publication["numeric_reproduction"] is True and publication["numerical_checks"] == 258,
            "v3 publication provenance disagrees")
    for name in ("activity", "features", "grouping", "inputs"):
        suffix = {"activity": "activity.json.gz", "features": "features.json.gz",
                  "grouping": "grouping.json", "inputs": "inputs.json"}[name]
        require(publication["files"]["core/" + suffix]["sha256"] == EXPECTED[name],
                f"publication file commitment disagrees: {name}")
    require(manifest["files"]["activity.json.gz"] == EXPECTED["activity"] and
            manifest["source"]["sha256"] == EXPECTED["raw_workbook"] and
            manifest["coverage"]["complete_required_ids"] == 55276 and
            manifest["coverage"]["unchanged_published_values"] is True and
            manifest["coverage"]["no_activity_based_selection"] is True and
            manifest["coverage"]["no_imputation"] is True and
            all(value is True for value in manifest["gates"].values()),
            "activity source manifest binding or gates disagree")
    require(protocol["protocol_id"] == "context-activity-v3" and
            protocol["source"]["jointdata_expected_sha256"] == EXPECTED["raw_workbook"] and
            protocol["source"]["worksheet"] == "all_cell_types_summary" and
            protocol["source"]["contexts"] == list(CONTEXTS), "v3 protocol source disagrees")
    require(SOURCE.stat().st_size == 53557805 and publication["source"]["bytes"] == 53557805,
            "raw source byte count disagrees")


def groups_and_folds(ids):
    data = read_json(GROUPING)
    records = records_by_id(data["records"], "group")
    require(set(records) == ids and data["protocol_sha256"] == EXPECTED["grouping_protocol"] and
            data["panel_sha256"] == EXPECTED["panel"], "group coverage or input differs")
    folds = {}
    for fold in data["folds"]:
        number = fold["fold"]
        for identifier in fold["test_ids"]:
            require(identifier in ids and identifier not in folds, "invalid or duplicate test-fold assignment")
            folds[identifier] = number
    require(set(folds) == ids, "test folds do not cover panel")
    group_folds = defaultdict(set)
    for identifier in ids:
        group_folds[records[identifier]["group"]].add(folds[identifier])
    require(all(len(values) == 1 for values in group_folds.values()), "a source group crosses test folds")
    return {identifier: (records[identifier]["group"], folds[identifier]) for identifier in ids}


def compute_design(panel):
    """Discover classes solely from source DNA before any activity join."""
    snapshot = read_json(FEATURES)
    require(snapshot["feature_names"] == FEATURE_NAMES, "unexpected v3 feature order")
    saved = records_by_id(snapshot["records"], "feature")
    require(set(saved) == set(panel), "feature snapshot ID coverage differs")
    classes = defaultdict(list)
    counts_by_id = {}
    fingerprints = {}
    for identifier, record in panel.items():
        counts = features(record["sequence"])
        old = saved[identifier]
        require(old["length"] == 200 and old["canonical_trimer_counts"] == list(counts[:32]) and
                old["gc_count"] == counts[32] and
                old["vector"] == [c / 198 for c in counts[:32]] + [counts[32] / 200],
                f"recomputed feature differs from v3 snapshot for {identifier}")
        fingerprint = feature_fingerprint(counts)
        require(fingerprint not in fingerprints or fingerprints[fingerprint] == counts,
                "different exact feature tuples share a fingerprint")
        fingerprints[fingerprint] = counts
        classes[counts].append(identifier)
        counts_by_id[identifier] = counts
    return classes, counts_by_id


def activity_values(ids):
    """Load only the already-explored, published v3 activity endpoint."""
    data = read_json(ACTIVITY)
    require(data["endpoint"] == ENDPOINT and data["protocol_sha256"] == EXPECTED["protocol"] and
            data["contexts"] == list(CONTEXTS), "activity endpoint or context differs")
    records = records_by_id(data["records"], "activity")
    require(set(records) == ids, "v3 activity coverage differs from panel")
    result = {}
    for identifier, record in records.items():
        obs = record["observations"]
        require(set(obs) == set(CONTEXTS), "activity context coverage differs")
        values = {}
        for context in CONTEXTS:
            value = obs[context]
            require(type(value) in (int, float) and math.isfinite(value), "invalid activity value")
            values[context] = float(value)
        result[identifier] = values
    return result


def audit(out):
    provenance()
    panel = read_panel()
    require(len(panel) == 55276, "unexpected panel size")
    metadata = groups_and_folds(set(panel))
    classes, counts_by_id = compute_design(panel)
    collision_classes = [(counts, sorted(ids)) for counts, ids in classes.items() if len(ids) > 1]
    collision_classes.sort(key=lambda item: item[1])
    # Do not use labels for grouping, eligibility or source-feature validation.
    activity = activity_values(set(panel))
    sse_by_cell = {cell: Fraction() for cell in CONTEXTS}
    pair_counts = Counter()
    members_with_collisions = sum(len(ids) for _, ids in collision_classes)
    collision_data = []
    for counts, ids in collision_classes:
        members = []
        for identifier in ids:
            source = panel[identifier]
            group, fold = metadata[identifier]
            members.append({k: source[k] for k in ("id", "chromosome", "start", "end", "parent_id",
                                                    "reference_match_strand", "source_declared_strand")} |
                           {"sequence_sha256": hashlib.sha256(source["sequence"].encode("ascii")).hexdigest(),
                            "group": group, "test_fold": fold,
                            "observations": {c: {"repr": repr(activity[identifier][c]),
                                                  "binary64_hex": activity[identifier][c].hex()}
                                             for c in CONTEXTS}})
        relations = []
        for a, b in itertools.combinations(ids, 2):
            s1, s2 = panel[a]["sequence"], panel[b]["sequence"]
            equal = s1 == s2
            is_rc = s1 == reverse_complement(s2)
            pair_counts["same_sequence"] += equal
            pair_counts["distinct_reverse_complement"] += not equal and is_rc
            pair_counts["distinct_non_reverse_complement"] += not equal and not is_rc
            same_interval = all(panel[a][x] == panel[b][x] for x in ("chromosome", "start", "end"))
            same_group, same_fold = metadata[a][0] == metadata[b][0], metadata[a][1] == metadata[b][1]
            pair_counts["same_interval"] += same_interval
            pair_counts["same_group"] += same_group
            pair_counts["same_test_fold"] += same_fold
            relations.append({"id1": a, "id2": b, "same_sequence": equal,
                              "reverse_complement": is_rc, "same_interval": same_interval,
                              "same_group": same_group, "same_test_fold": same_fold})
        cell_results = {}
        for cell in CONTEXTS:
            y = [Fraction.from_float(activity[identifier][cell]) for identifier in ids]
            sse = squared_error_sum(y)
            sse_by_cell[cell] += sse
            cell_results[cell] = {"squared_error_sum": rational(sse),
                                  "class_mean_mse_floor": rational(sse / len(ids)),
                                  "pair_differences": [{"id1": ids[i], "id2": ids[j],
                                                        "first_minus_second": rational(y[i] - y[j])}
                                                       for i, j in itertools.combinations(range(len(ids)), 2)]}
        collision_data.append({"feature_counts": list(counts), "members": members,
                               "pair_relations": relations, "per_cell": cell_results})
    n = len(panel)
    per_cell = {c: {"collision_label_sse": rational(sse_by_cell[c]),
                    "equal_element_panel_mse_floor": rational(sse_by_cell[c] / n),
                    "collision_only_mse_floor": rational(sse_by_cell[c] / members_with_collisions)
                    if members_with_collisions else None} for c in CONTEXTS}
    pooled_sse = sum(sse_by_cell.values(), Fraction())
    pooled = {"equal_element_panel_mse_floor": rational(pooled_sse / (n * len(CONTEXTS))),
              "collision_only_mse_floor": rational(pooled_sse / (members_with_collisions * len(CONTEXTS)))
              if members_with_collisions else None}
    sequences = [item["sequence"] for item in panel.values()]
    summary = {
        "schema_version": 1, "endpoint": ENDPOINT, "contexts": list(CONTEXTS),
        "scope": "post-hoc published v3 activity labels; no new assay outcome source decoded",
        "input_sha256": EXPECTED, "source_workbook_values_reparsed": False,
        "records": n, "features_recomputed_and_snapshot_exact": n,
        "distinct_exact_sequences": len(set(sequences)),
        "distinct_reverse_complement_sequences": len({min(s, reverse_complement(s)) for s in sequences}),
        "feature_classes": len(classes), "multiple_member_feature_classes": len(collision_classes),
        "members_in_multiple_member_classes": members_with_collisions,
        "collision_class_size_distribution": dict(sorted(Counter(str(len(ids)) for _, ids in collision_classes).items())),
        "equal_feature_pair_relations": dict(sorted(pair_counts.items())),
        "per_cell": per_cell, "pooled": pooled,
        "arithmetic": "exact rational arithmetic on JSON-decoded binary64 activity labels",
        "weighting": "equal elements and equal cells, as in published context-activity-v3",
    }
    out.mkdir(parents=True, exist_ok=False)
    assignments_path = out / "assignments.jsonl.gz"
    with assignments_path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="\n") as f:
                for identifier in sorted(panel):
                    counts = counts_by_id[identifier]
                    group, fold = metadata[identifier]
                    item = {"id": identifier, "feature_fingerprint": feature_fingerprint(counts),
                            "class_size": len(classes[counts]), "group": group, "test_fold": fold}
                    f.write(json.dumps(item, sort_keys=True, separators=(",", ":")) + "\n")
    collisions_path, summary_path = out / "collisions.json", out / "summary.json"
    collisions_path.write_text(json.dumps({"schema_version": 1, "classes": collision_data},
                                          indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {"implementation_sha256": sha(Path(__file__)),
                "assignments_sha256": sha(assignments_path), "collisions_sha256": sha(collisions_path),
                "summary_sha256": sha(summary_path)}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("records", "feature_classes", "multiple_member_feature_classes",
                     "members_in_multiple_member_classes", "equal_feature_pair_relations", "pooled")}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="new output directory")
    args = parser.parse_args()
    audit(args.out)


if __name__ == "__main__":
    main()
