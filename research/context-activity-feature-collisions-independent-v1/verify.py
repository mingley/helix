#!/usr/bin/env python3
"""Independently verify a post-hoc collision in the explored v3 activity panel.

The labels checked here are already-published v3 values from the same
experiment. They are not fresh validation or independent laboratory evidence.
"""
import argparse
from collections import Counter, defaultdict
from decimal import Decimal
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile


ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "datasets/context-panel/panel.json"
CORE = ROOT / "research/results/context-activity-v3/core"
S7 = ROOT / "datasets/context-panel/raw/Table-S7-joint-design.xlsx"
S11 = ROOT / "datasets/context-endpoint-recovery/raw/jointdata.xlsx"
PRIMARY = ROOT / "research/context-activity-feature-collisions-v1/results"
HASHES = {
    "panel": "fce6b2a1aaf5aa4d530a2e298dccb840bc7c594a05b6831fdb0c6c0f9e2a9be3",
    "features": "1a9b093ca7d0ae56b19294a53c427f0d7d00ee4821ff09e03d8384d601b9417a",
    "activity": "83b2df30128f3ea8326d6d5d1de73b36e32ced113cccebebc915d4fc06e8054b",
    "grouping": "333069c4ab41ed8ea6976b0bc50cd08de9a33118ab60d683533b7ee32a6f517f",
    "study": "ac45559e8f16a094564250cc31a424f060791193f95cbdd9d4634fd37bf317d7",
    "design_xlsx": "f9ec4c6c2d230b0ef507a39e5d5c73b3b36f90ccc5758c33b442c7561e20e7e5",
    "activity_xlsx": "9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e",
    "inputs": "fe78220d9f94b5022582e79e02ed8d1dbd927c8d292825655ac3d018d04ef6d1",
    "publication": "39cd6a1902ed3b7a4709294b293590309af5f38d2d8c578f82922da74a63a335",
    "source_manifest": "6a3371bed4061f9ca9fa1a95bfc0947ffee91daa22fe509063708e549d453730",
    "protocol": "64e14037fb1f42ddc1ecaca1d23bb8096fabac1f1e892d105eaa20f525fbbb20",
    "grouping_protocol": "8f0576845fecb5cfcf9bfa354fab446d849244448a3489905a6ce33bdeb21cec",
}
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
DNA = "ACGT"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
CONTEXTS = ("HepG2", "K562", "WTC11")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def reverse_complement(sequence):
    return "".join(COMPLEMENT[base] for base in reversed(sequence))


def canonical(word):
    return min(word, reverse_complement(word))


CLASSES = sorted({canonical("".join(word)) for word in itertools.product(DNA, repeat=3)})
CLASS_INDEX = {word: i for i, word in enumerate(CLASSES)}
INDEX_FOR_RAW = {"".join(word): CLASS_INDEX[canonical("".join(word))]
                 for word in itertools.product(DNA, repeat=3)}


def vector(sequence):
    """Compute integer counts from complete DNA, with no published feature reader."""
    require(isinstance(sequence, str) and len(sequence) >= 3 and not set(sequence) - set(DNA),
            "invalid literal DNA sequence")
    counts = [0] * len(CLASSES)
    for i in range(len(sequence) - 2):
        counts[INDEX_FOR_RAW[sequence[i:i + 3]]] += 1
    return tuple(counts) + (sequence.count("C") + sequence.count("G"),)


def floor_for_equal_predictions(values):
    """Smallest within-class SSE, treating each observed binary64 value exactly."""
    require(len(values) >= 2 and all(type(x) is float for x in values),
            "floor requires at least two binary64 observations")
    exact = [Fraction.from_float(value) for value in values]
    mean = sum(exact, Fraction()) / len(exact)
    return sum(((value - mean) ** 2 for value in exact), Fraction())


def fraction_object(value):
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def shared_strings(archive):
    result = []
    with archive.open("xl/sharedStrings.xml") as stream:
        for _, elem in ET.iterparse(stream, events=("end",)):
            if elem.tag == NS + "si":
                result.append("".join(node.text or "" for node in elem.iter(NS + "t")))
                elem.clear()
    return result


def cells_for_row(row, strings, wanted=None):
    result = {}
    rownum = int(row.attrib["r"])
    for cell in row.findall(NS + "c"):
        address = re.fullmatch(r"([A-Z]+)([0-9]+)", cell.attrib.get("r", ""))
        require(address is not None and int(address.group(2)) == rownum,
                "bad source cell address")
        column = address.group(1)
        if wanted is not None and column not in wanted:
            continue
        require(column not in result and cell.find(NS + "f") is None,
                "duplicate or formula source cell")
        value = cell.find(NS + "v")
        require(value is not None and value.text is not None, "missing source cell value")
        typ = cell.attrib.get("t", "n")
        if typ == "s":
            index = int(value.text)
            require(0 <= index < len(strings), "bad source shared-string index")
            result[column] = ("s", strings[index])
        else:
            require(typ == "n", "unexpected source cell type")
            result[column] = ("n", value.text)
    return result


def workbook_rows(path, sheet_name, selected_rows=None, selected_ids=None, header_row=1):
    """Read only needed sheet cells; XML parsing can instantiate other row content."""
    require((selected_rows is None) != (selected_ids is None), "choose row numbers or IDs")
    found = {}
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheets = workbook.findall(".//" + NS + "sheet")
        require(sheets and sheets[0].attrib.get("name") == sheet_name,
                "unexpected first source worksheet")
        relationship = sheets[0].get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        matching = [r for r in rels if r.get("Id") == relationship]
        require(len(matching) == 1 and matching[0].get("Target") == "worksheets/sheet1.xml",
                "source worksheet relationship differs")
        strings = shared_strings(archive)
        with archive.open("xl/worksheets/sheet1.xml") as source:
            for _, row in ET.iterparse(source, events=("end",)):
                if row.tag != NS + "row":
                    continue
                number = int(row.attrib["r"])
                if number == header_row:
                    found["header"] = cells_for_row(row, strings)
                elif number > header_row:
                    if selected_rows is not None:
                        selected = number in selected_rows
                    else:
                        one = cells_for_row(row, strings, {"A"})
                        selected = one.get("A", (None, None))[1] in selected_ids
                    if selected:
                        cells = cells_for_row(row, strings)
                        identity = cells["A"][1]
                        require(identity not in found, "duplicate selected workbook ID")
                        found[identity] = (number, cells)
                row.clear()
    require("header" in found, "source header missing")
    return found


def source_activity(path, ids):
    found = workbook_rows(path, "all_cell_types_summary", selected_ids=set(ids))
    expected_headers = {"A": "name", "B": "HepG2 [log2(rna/dna)]",
                        "C": "K562 [log2(rna/dna)]", "D": "WTC11 [log2(rna/dna)]"}
    require({k: v[1] for k, v in found.pop("header").items()} == expected_headers,
            "unexpected original activity headers")
    require(set(found) == set(ids), "original activity rows missing")
    output = {}
    for identity, (row, cells) in found.items():
        values = {}
        for col, context in zip("BCD", CONTEXTS):
            require(col in cells and cells[col][0] == "n", "source activity is not numeric")
            raw = cells[col][1]
            number = Decimal(raw)
            require(number.is_finite(), "nonfinite original activity")
            values[context] = raw
        output[identity] = {"row": row, "lexemes": values}
    return output


def source_design(path, records):
    found = workbook_rows(path, "joint library", selected_ids={row["id"] for row in records}, header_row=2)
    headers = found.pop("header")
    for col, header in {"A": "name", "G": "chr.hg38", "H": "start.hg38", "I": "stop.hg38",
                        "J": "str.hg38", "K": "230nt sequence (15nt 5' adaptor - 200nt element - 15nt 3' adaptor)"}.items():
        require(headers[col][1] == header, "unexpected original design header")
    require(set(found) == {row["id"] for row in records}, "original design rows missing")
    result = {}
    for row in records:
        number, cells = found[row["id"]]
        oligo = cells["K"][1]
        require(number == row["source_metadata"]["s7_row_1based"] and len(oligo) == 230 and
                oligo[15:215] == row["sequence"], "original design oligo differs from panel DNA")
        require("chr" + cells["G"][1] == row["chromosome"] and
                int(cells["H"][1]) == row["start"] and int(cells["I"][1]) == row["end"] and
                cells["J"][1] == row["source_metadata"]["source_declared_strand"],
                "original design coordinates or strand differ")
        result[row["id"]] = {"row": number, "oligo": oligo, "strand": cells["J"][1]}
    return result


def load_json(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return json.load(stream)


def keyed(rows, name):
    output = {}
    for row in rows:
        identity = row["id"]
        require(isinstance(identity, str) and identity and identity not in output,
                f"duplicate or invalid {name} ID")
        output[identity] = row
    return output


def check_equal(actual, expected, name):
    if actual != expected:
        raise ValueError(f"{name} differs")


def fingerprint(counts):
    payload = json.dumps(list(counts), separators=(",", ":")).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def source_bindings():
    paths = {
        "panel": PANEL, "features": CORE / "features.json.gz", "activity": CORE / "activity.json.gz",
        "grouping": CORE / "grouping.json", "study": CORE / "study.json.gz",
        "design_xlsx": S7, "activity_xlsx": S11, "inputs": CORE / "inputs.json",
        "publication": ROOT / "research/results/context-activity-v3/publication.json",
        "source_manifest": ROOT / "datasets/context-activity/manifest.json",
        "protocol": ROOT / "experiments/context-activity-v3.json",
        "grouping_protocol": ROOT / "experiments/context-specificity-v2.json",
    }
    for name, path in paths.items():
        require(sha256(path) == HASHES[name], f"{name} source hash mismatch")
    inputs = load_json(paths["inputs"])
    publication = load_json(paths["publication"])
    manifest = load_json(paths["source_manifest"])
    protocol = load_json(paths["protocol"])
    require(inputs["activity_sha256"] == HASHES["activity"] and
            inputs["dna_panel_sha256"] == HASHES["panel"] and
            inputs["grouping_sha256"] == HASHES["grouping"] and
            inputs["raw_source_sha256"] == HASHES["activity_xlsx"] and
            inputs["source_manifest_sha256"] == HASHES["source_manifest"] and
            inputs["protocol_sha256"] == HASHES["protocol"], "v3 input provenance mismatch")
    require(publication["published_activity_sha256"] == HASHES["activity"] and
            publication["source"]["sha256"] == HASHES["activity_xlsx"] and
            publication["source_manifest_sha256"] == HASHES["source_manifest"] and
            publication["protocol_sha256"] == HASHES["protocol"] and
            publication["literal_activity_binding_verified"] is True,
            "v3 publication provenance mismatch")
    require(manifest["source"]["sha256"] == HASHES["activity_xlsx"] and
            manifest["files"]["activity.json.gz"] == HASHES["activity"] and
            protocol["source"]["jointdata_expected_sha256"] == HASHES["activity_xlsx"] and
            protocol["source"]["worksheet"] == "all_cell_types_summary",
            "source workbook binding mismatch")
    return paths


def groups_and_folds(panel):
    data = load_json(CORE / "grouping.json")
    groups = keyed(data["records"], "group")
    require(set(groups) == set(panel), "group membership does not cover panel")
    by_group = defaultdict(set)
    for identity, group in groups.items():
        by_group[group["group"]].add(identity)
    fold_for_group = {}
    for fold in data["folds"]:
        for group in fold["test_groups"]:
            require(group in by_group and group not in fold_for_group, "duplicate/unknown fold group")
            fold_for_group[group] = fold["fold"]
        expected = set().union(*(by_group[group] for group in fold["test_groups"]))
        require(expected == set(fold["test_ids"]) and len(expected) == len(fold["test_ids"]),
                "fold ID and group membership differ")
    require(set(fold_for_group) == set(by_group), "fold groups do not cover panel")
    return {identity: (row["group"], fold_for_group[row["group"]]) for identity, row in groups.items()}


def observed_values(panel):
    source = load_json(CORE / "activity.json.gz")
    require(source["contexts"] == list(CONTEXTS) and
            source["endpoint"] == "published_joint_log2_rna_dna_activity",
            "unexpected activity endpoint")
    data = keyed(source["records"], "activity")
    require(set(data) == set(panel), "activity does not cover panel")
    result = {}
    for identity, row in data.items():
        require(set(row["observations"]) == set(CONTEXTS), "activity cells differ")
        values = {}
        for context, value in row["observations"].items():
            require(type(value) in (float, int) and math.isfinite(value), "invalid activity value")
            values[context] = float(value)
        result[identity] = values
    return result


def recompute_features(panel):
    saved = load_json(CORE / "features.json.gz")
    require(saved["feature_names"] == ["rc3mer_" + c for c in CLASSES] + ["gc_fraction"],
            "unexpected feature order")
    old = keyed(saved["records"], "feature")
    require(set(old) == set(panel), "published feature coverage differs")
    values, classes, digest_to_vector = {}, defaultdict(list), {}
    for identity, row in panel.items():
        sequence = row["sequence"]
        require(len(sequence) == 200, "panel is not fixed length")
        counts = vector(sequence)
        previous = old[identity]
        require(previous["length"] == 200 and previous["canonical_trimer_counts"] == list(counts[:-1]) and
                previous["gc_count"] == counts[-1] and
                previous["vector"] == [x / 198 for x in counts[:-1]] + [counts[-1] / 200],
                f"source features disagree for {identity}")
        digest = fingerprint(counts)
        require(digest not in digest_to_vector or digest_to_vector[digest] == counts,
                "different exact vectors share a feature fingerprint")
        digest_to_vector[digest] = counts
        values[identity] = counts
        classes[counts].append(identity)
    return values, classes


def expected_collision_records(panel, classes, membership, activity):
    pair_counts = Counter()
    totals = {c: Fraction() for c in CONTEXTS}
    output = []
    for counts, identities in sorted(((k, sorted(v)) for k, v in classes.items() if len(v) > 1),
                                    key=lambda entry: entry[1]):
        members = []
        for identity in identities:
            row = panel[identity]
            meta = row["source_metadata"]
            group, fold = membership[identity]
            members.append({"id": identity, "chromosome": row["chromosome"], "start": row["start"],
                            "end": row["end"], "parent_id": row["parent_id"],
                            "reference_match_strand": meta["reference_match_strand"],
                            "source_declared_strand": meta["source_declared_strand"],
                            "sequence_sha256": hashlib.sha256(row["sequence"].encode("ascii")).hexdigest(),
                            "group": group, "test_fold": fold,
                            "observations": {cell: {"repr": repr(activity[identity][cell]),
                                                     "binary64_hex": activity[identity][cell].hex()}
                                             for cell in CONTEXTS}})
        relations = []
        for a, b in itertools.combinations(identities, 2):
            first, second = panel[a], panel[b]
            same = first["sequence"] == second["sequence"]
            rc = first["sequence"] == reverse_complement(second["sequence"])
            interval = all(first[k] == second[k] for k in ("chromosome", "start", "end"))
            same_group = membership[a][0] == membership[b][0]
            same_fold = membership[a][1] == membership[b][1]
            pair_counts["same_sequence"] += same
            pair_counts["distinct_reverse_complement"] += not same and rc
            pair_counts["distinct_non_reverse_complement"] += not same and not rc
            pair_counts["same_interval"] += interval
            pair_counts["same_group"] += same_group
            pair_counts["same_test_fold"] += same_fold
            relations.append({"id1": a, "id2": b, "same_sequence": same,
                              "reverse_complement": rc, "same_interval": interval,
                              "same_group": same_group, "same_test_fold": same_fold})
        per_cell = {}
        for cell in CONTEXTS:
            y = [activity[identity][cell] for identity in identities]
            sse = floor_for_equal_predictions(y)
            totals[cell] += sse
            differences = [{"id1": identities[i], "id2": identities[j],
                            "first_minus_second": fraction_object(Fraction.from_float(y[i]) -
                                                                    Fraction.from_float(y[j]))}
                           for i, j in itertools.combinations(range(len(identities)), 2)]
            per_cell[cell] = {"squared_error_sum": fraction_object(sse),
                              "class_mean_mse_floor": fraction_object(sse / len(identities)),
                              "pair_differences": differences}
        output.append({"feature_counts": list(counts), "members": members,
                       "pair_relations": relations, "per_cell": per_cell})
    return output, pair_counts, totals


def verify_predictions(ids, membership, activity):
    study = load_json(CORE / "study.json.gz")
    matches = {}
    for record in study["predictions"]:
        identity, cell = record["sequence_id"], record["context"]
        if identity not in ids:
            continue
        key = (identity, cell)
        require(key not in matches and cell in CONTEXTS and
                record["observed"] == activity[identity][cell] and
                (record["group"], record["fold_index_0based"]) == membership[identity],
                "published prediction record binding differs")
        matches[key] = record["predictions"]
    require(set(matches) == set(itertools.product(ids, CONTEXTS)), "published predictions missing")
    for cell in CONTEXTS:
        predictions = [matches[(identity, cell)] for identity in sorted(ids)]
        require(set(predictions[0]) == {"global_mean", "context_mean", "shared_dna",
                                        "additive_context_dna", "independent_context_dna"},
                "unexpected saved model names")
        require(all(item == predictions[0] for item in predictions[1:]),
                f"saved model predictions differ within equal-feature class in {cell}")
    return matches


def verify_assignments(path, panel, features, classes, membership):
    seen = set()
    with gzip.open(path, "rt", encoding="utf-8") as source:
        for line in source:
            actual = json.loads(line)
            identity = actual["id"]
            require(identity in panel and identity not in seen, "unknown or duplicate assignment")
            seen.add(identity)
            group, fold = membership[identity]
            expected = {"id": identity, "feature_fingerprint": fingerprint(features[identity]),
                        "class_size": len(classes[features[identity]]), "group": group, "test_fold": fold}
            check_equal(actual, expected, f"assignment {identity}")
    require(seen == set(panel), "assignment records do not cover panel")
    return len(seen)


def verify_primary(primary, panel, membership, features, classes, activity):
    output, pair_counts, totals = expected_collision_records(panel, classes, membership, activity)
    check_equal(load_json(primary / "collisions.json"), {"schema_version": 1, "classes": output},
                "full collision records")
    colliding = sum(len(group) for group in classes.values() if len(group) > 1)
    n = len(panel)
    per_cell = {cell: {"collision_label_sse": fraction_object(totals[cell]),
                       "equal_element_panel_mse_floor": fraction_object(totals[cell] / n),
                       "collision_only_mse_floor": fraction_object(totals[cell] / colliding)
                       if colliding else None} for cell in CONTEXTS}
    total = sum(totals.values(), Fraction())
    pooled = {"equal_element_panel_mse_floor": fraction_object(total / (n * len(CONTEXTS))),
              "collision_only_mse_floor": fraction_object(total / (colliding * len(CONTEXTS)))
              if colliding else None}
    inputs = {key: value for key, value in HASHES.items()
              if key not in ("design_xlsx", "activity_xlsx", "study")}
    inputs["raw_workbook"] = HASHES["activity_xlsx"]
    sequences = [row["sequence"] for row in panel.values()]
    expected = {
        "schema_version": 1, "endpoint": "published_joint_log2_rna_dna_activity",
        "contexts": list(CONTEXTS),
        "scope": "post-hoc published v3 activity labels; no new assay outcome source decoded",
        "input_sha256": inputs, "source_workbook_values_reparsed": False,
        "records": n, "features_recomputed_and_snapshot_exact": n,
        "distinct_exact_sequences": len(set(sequences)),
        "distinct_reverse_complement_sequences": len({min(seq, reverse_complement(seq)) for seq in sequences}),
        "feature_classes": len(classes), "multiple_member_feature_classes": len(output),
        "members_in_multiple_member_classes": colliding,
        "collision_class_size_distribution": dict(sorted(Counter(str(len(v)) for v in classes.values()
                                                                    if len(v) > 1).items())),
        "equal_feature_pair_relations": dict(sorted(pair_counts.items())),
        "per_cell": per_cell, "pooled": pooled,
        "arithmetic": "exact rational arithmetic on JSON-decoded binary64 activity labels",
        "weighting": "equal elements and equal cells, as in published context-activity-v3",
    }
    check_equal(load_json(primary / "summary.json"), expected, "complete primary summary")
    manifest = load_json(primary / "manifest.json")
    for name in ("assignments", "collisions", "summary"):
        filename = "assignments.jsonl.gz" if name == "assignments" else name + ".json"
        require(manifest[name + "_sha256"] == sha256(primary / filename), f"primary {name} hash mismatch")
    require(manifest["implementation_sha256"] == sha256(primary.parent / "audit.py"),
            "primary implementation hash mismatch")
    assignments = verify_assignments(primary / "assignments.jsonl.gz", panel, features, classes, membership)
    return output, expected, assignments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, default=PRIMARY)
    parser.add_argument("--receipt", type=Path, required=True, help="new receipt path")
    args = parser.parse_args()
    if args.receipt.exists() or args.receipt.is_symlink():
        raise FileExistsError("receipt path already exists")
    source_bindings()
    raw_panel = load_json(PANEL)
    panel = keyed(raw_panel["records"], "panel")
    require(len(panel) == 55276, "unexpected panel size")
    membership = groups_and_folds(panel)
    features, classes = recompute_features(panel)
    activity = observed_values(panel)
    collision, summary, assignments = verify_primary(args.primary, panel, membership, features, classes, activity)
    require(len(collision) == 1 and len(collision[0]["members"]) == 2,
            "this checker expects the one two-member collision it reports")
    identities = [row["id"] for row in collision[0]["members"]]
    pair = [panel[identity] for identity in identities]
    require(pair[0]["sequence"] != pair[1]["sequence"] and
            pair[0]["sequence"] == reverse_complement(pair[1]["sequence"]),
            "collision sequences are not distinct reverse complements")
    require(all(pair[0][k] == pair[1][k] for k in ("chromosome", "start", "end")),
            "collision coordinates differ")
    design = source_design(S7, pair)
    require({design[identity]["strand"] for identity in identities} == {"+", "-"},
            "original source strands are not opposite")
    left = {design[identity]["oligo"][:15] for identity in identities}
    right = {design[identity]["oligo"][215:] for identity in identities}
    require(len(left) == len(right) == 1, "design adapters differ")
    require(design[identities[0]]["oligo"] != reverse_complement(design[identities[1]]["oligo"]),
            "unexpected whole-oligo reverse-complement relationship")
    original = source_activity(S11, identities)
    for identity in identities:
        for cell in CONTEXTS:
            require(float(original[identity]["lexemes"][cell]).hex() == activity[identity][cell].hex(),
                    f"original source and v3 binary64 differ: {identity}/{cell}")
    verify_predictions(set(identities), membership, activity)
    receipt = {
        "method": "Independent raw-DNA 3-mer counting, exact binary64 rational floors, direct source XLSX ZIP/XML reads, and saved-fold prediction check",
        "scope": "post-hoc reuse of already-explored v3 labels; no independent assay or biological validation",
        "independent_implementation_sha256": sha256(Path(__file__).resolve()),
        "primary_implementation_sha256": sha256(args.primary.parent / "audit.py"),
        "primary_summary_sha256": sha256(args.primary / "summary.json"),
        "primary_collisions_sha256": sha256(args.primary / "collisions.json"),
        "primary_assignments_sha256": sha256(args.primary / "assignments.jsonl.gz"),
        "inputs_sha256": HASHES,
        "panel_records_checked": len(panel), "published_feature_records_checked": len(features),
        "primary_assignment_records_checked": assignments,
        "complete_collision_and_summary_match": True,
        "collision_member_ids": identities,
        "same_interval": True, "distinct_reverse_complement_cores": True,
        "source_design_rows_and_strands": {identity: {"row": design[identity]["row"],
                                                       "strand": design[identity]["strand"]}
                                            for identity in identities},
        "identical_source_15nt_adapters": {"left": next(iter(left)), "right": next(iter(right))},
        "complete_230nt_designs_are_reverse_complements": False,
        "source_activity_rows_and_raw_lexemes": original,
        "source_lexemes_decode_to_exact_retained_binary64": True,
        "source_numeric_precision": "Original XML numeric lexemes may differ from the v3 JSON decimal spellings; each selected lexeme decodes to the same binary64, and floors use those binary64 values exactly as rationals.",
        "source_XML_exposure": "S11 streaming XML instantiates other row contents and the workbook-wide shared-string table; only selected IDs' six aggregate numeric cells are interpreted. All v3 activity labels were previously explored and are loaded by this verifier.",
        "all_five_saved_predictions_match_pairwise_in_each_cell": True,
        "exact_pooled_panel_floor": summary["pooled"]["equal_element_panel_mse_floor"],
        "exact_collision_only_floor": summary["pooled"]["collision_only_mse_floor"],
        "all_checks_passed": True,
    }
    with args.receipt.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
