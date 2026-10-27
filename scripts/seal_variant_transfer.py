#!/usr/bin/env python3
"""Seal an outcome-free source/component cohort before reading Siraj S4 effects.

Global components are immutable. Technical-control matches exclude complete
components; S4 ID membership then removes individual unavailable measurements
without regrouping. No workbook or reporter-effect value is decoded here.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_SHA = "811930c2d3a5ee763014ea611acda4d2c84404b7d5ae399c7843c62d0a828149"
SOURCE_MANIFEST_SHA = "61b9756522156b2406e692113cf1268df525125441a6d348cfa99f3795875117"
MODEL_SHA = "2ea476bc5413dc13486dec9f3d7e9162a50a0e47f1d83ce29efabc5f5cc5af9d"
RAW = {"S2": ("siraj2026-S2-acquisition.json", "source_S2_raw_cache", 40795798,
               "4b6e20d1789021ee8f11f57af366a65680d25496de9351a9e781e9424716d824"),
       "S4": ("siraj2026-S4-acquisition.json", "source_S4_raw_cache", 224405074,
               "723c4a0e1cfdb2c9c6ae55b3fa9206855dd95e0241319b9fc8d0831d3b42f9c2"),
       "S30": ("S30-bundle-acquisition.json", "source_S30_raw_cache", 19928224,
                "3e62a3876b70e9c75951df26ed09ba552c3ba6f9f3b318934c1b01a5a42403b1")}
HEX = re.compile(r"[0-9a-f]{64}\Z")
HEADERS = {"A": "Variant", "AN": "log2Skew in HEPG2", "AO": "log2Skew in K562"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def strict_object(items):
    result = {}
    for key, value in items:
        require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def decode(text):
    return json.loads(text, object_pairs_hook=strict_object,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def read_json(path):
    return decode(Path(path).read_text())


def read_jsonl(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            yield decode(line)


def encode(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def write_jsonl(path, rows):
    with Path(path).open("xb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as f:
            for row in rows:
                f.write((json.dumps(row, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode())


def member(root, name):
    require(isinstance(name, str), "manifest member must be text")
    rel = PurePosixPath(name)
    require(not rel.is_absolute() and rel.parts and ".." not in rel.parts
            and str(rel) == name and "\\" not in name, "unsafe manifest member")
    target = Path(root) / name
    require(not Path(root).is_symlink(), "symlink inventory root")
    cursor = target
    while cursor != Path(root):
        require(not cursor.is_symlink(), "symlink inventory member")
        cursor = cursor.parent
    require(target.is_file() and target.resolve().is_relative_to(Path(root).resolve()), "missing inventory member")
    return target


def verify_inventory(root, manifest_name, expected, *, exact=False):
    require(isinstance(expected, str) and HEX.fullmatch(expected), "invalid expected inventory SHA256")
    root = Path(root)
    manifest_path = member(root, manifest_name)
    require(sha256(manifest_path) == expected, "inventory manifest SHA256 mismatch")
    manifest = read_json(manifest_path)
    require(type(manifest.get("schema_version")) is int and manifest["schema_version"] in {1, 5}
            and isinstance(manifest.get("files"), dict), "invalid inventory schema")
    for name, receipt in manifest["files"].items():
        path = member(root, name)
        digest = receipt if isinstance(receipt, str) else receipt.get("sha256")
        require(isinstance(digest, str) and HEX.fullmatch(digest) and sha256(path) == digest,
                "inventory member SHA256 mismatch: " + name)
        if not isinstance(receipt, str):
            require(type(receipt.get("bytes")) is int and path.stat().st_size == receipt["bytes"],
                    "inventory member byte count differs: " + name)
    if exact:
        actual = {str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()}
        require(actual == set(manifest["files"]) | {manifest_name}, "uncommitted extra or missing inventory files")
    return manifest


def component_identity(ids):
    return hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()


def select_cohort(components, pairs, source_ids, controls):
    """Pure metadata selection; original component labels never change."""
    assignments, component_rows = {}, []
    for row in components:
        ids = row.get("variant_ids")
        require(isinstance(ids, list) and ids and all(isinstance(v, str) and v for v in ids)
                and len(set(ids)) == len(ids), "invalid or duplicate component membership")
        require(row.get("component") == component_identity(ids) and row.get("size") == len(ids),
                "component identity or size differs")
        require(type(row.get("retained_after_prior_exposure")) is bool
                and isinstance(row.get("exclusion_reasons"), list), "invalid primary component gate")
        require(row["retained_after_prior_exposure"] == (not row["exclusion_reasons"]),
                "contradictory primary component gate")
        reasons = list(row["exclusion_reasons"])
        matched_controls = sorted(set(ids) & controls)
        if matched_controls:
            reasons.append("entire_component_contains_reused_S30_control")
        admitted = sorted(set(ids) & source_ids) if not reasons else []
        component_rows.append({"component": row["component"], "original_size": len(ids),
            "retained_after_prior_exposure": row["retained_after_prior_exposure"],
            "S30_control_variant_ids": matched_controls, "component_exclusion_reasons": sorted(set(reasons)),
            "S4_member_count": len(set(ids) & source_ids), "final_cohort_member_count": len(admitted)})
        for identity in ids:
            require(identity not in assignments, "variant belongs to multiple original components")
            assignments[identity] = component_rows[-1]
    require(controls <= assignments.keys(), "technical-control ID outside complete external population")
    pair_map = {}
    for row in pairs:
        identity = row["id"]
        require(identity not in pair_map and identity in assignments, "duplicate or unknown retained DNA ID")
        component = assignments[identity]
        require(component["retained_after_prior_exposure"] and row.get("component") == component["component"],
                "retained DNA changes original component or prior gate")
        a, b = row.get("sequence_ref"), row.get("sequence_alt")
        require(isinstance(a, str) and isinstance(b, str) and len(a) == len(b) == 200
                and set(a + b) <= set("ACGT") and [i for i in range(200) if a[i] != b[i]] == [99],
                "invalid centered200nt source DNA")
        gc = (a[99] in "GC") == (b[99] in "GC")
        require(type(row.get("gc_preserving")) is bool and row["gc_preserving"] == gc
                and isinstance(row.get("locus"), str) and row["locus"], "invalid source DNA metadata")
        pair_map[identity] = row
    expected_retained = {identity for identity, row in assignments.items() if row["retained_after_prior_exposure"]}
    require(pair_map.keys() == expected_retained, "primary retained DNA does not cover complete retained components")
    cohort = []
    for identity in sorted(pair_map):
        source = pair_map[identity]
        if assignments[identity]["component_exclusion_reasons"] or identity not in source_ids:
            continue
        cohort.append({"id": identity, "component": source["component"], "locus": source["locus"],
                       "sequence_ref": source["sequence_ref"], "sequence_alt": source["sequence_alt"],
                       "position_0based": 99, "gc_preserving": source["gc_preserving"]})
    return cohort, assignments, sorted(component_rows, key=lambda row: row["component"])


def check_group_gate(summary):
    require(summary.get("grouping_and_prior_exposure_gates_passed") is True
            and summary.get("phase") == "DNA_components_and_prior_exposure_complete_no_outcomes",
            "primary grouping/prior-exposure gate failed or incomplete")


def check_comparison(comparison, primary_manifest, independent_manifest):
    checks = comparison.get("checks")
    require(comparison.get("matched") is True and comparison.get("failed_checks") == []
            and comparison.get("new_numeric_outcome_access") is False
            and isinstance(checks, list) and checks
            and all(c.get("passed") is True for c in checks)
            and comparison.get("checks_passed") == len(checks), "independent grouping comparison failed or incomplete")
    require(comparison.get("bindings") == {"primary_manifest_sha256": primary_manifest,
            "independent_manifest_sha256": independent_manifest}, "independent comparison manifest bindings differ")


def finish_manifest(out):
    payload = {"schema_version": 1, "files": {str(p.relative_to(out)): {
        "bytes": p.stat().st_size, "sha256": sha256(p)} for p in sorted(out.rglob("*"))
        if p.is_file() and p.name != "manifest.json"}}
    (out / "manifest.json").write_bytes(encode(payload))


def seal(args):
    require(not args.out.exists(), "sealing output already exists")
    source_manifest = verify_inventory(args.source_root, "manifest.json", SOURCE_MANIFEST_SHA)
    group_manifest = verify_inventory(args.groups, "manifest.json", args.expected_group_manifest_sha256, exact=True)
    required = {"components.json", "eligible-pairs.jsonl.gz", "variant-ledger.jsonl.gz", "summary.json",
                "inputs.json", "parameters.json", "preparation.json", "implementation.py"}
    require(required <= group_manifest["files"].keys() and "failure.json" not in group_manifest["files"],
            "grouping outputs missing, incomplete or failed")
    summary = read_json(args.groups / "summary.json")
    check_group_gate(summary)
    inputs, parameters = read_json(args.groups / "inputs.json"), read_json(args.groups / "parameters.json")
    require(parameters.get("protocol_id") == "variant-transfer-v1" and parameters.get("bin_width") == 1000000
            and parameters.get("guard_bp_inclusive") == 10000 and parameters.get("max_hamming_mismatches") == 20,
            "frozen grouping parameters differ")
    prior = inputs["prior_inventory"]
    verify_inventory(Path(prior["path"]), "manifest.json", prior["manifest_sha256"])
    preparation = read_json(args.groups / "preparation.json")
    require(preparation.get("phase") == "prepared_not_executed", "invalid grouping preparation")
    for name, digest in preparation["files"].items():
        require(sha256(member(args.groups, name)) == digest, "group preparation file changed")
    require(summary["preparation_sha256"] == sha256(args.groups / "preparation.json"), "group preparation cross-binding differs")
    require(sha256(args.protocol) == PROTOCOL_SHA, "frozen semantic protocol bytes changed")
    verify_inventory(args.model_root, "artifacts.json", MODEL_SHA)
    for label, (receipt_name, cache_key, size, digest) in RAW.items():
        receipt = read_json(args.source_root / receipt_name)
        raw_path = Path(source_manifest[cache_key])
        require(receipt.get("sha256") == digest and receipt.get("bytes") == size
                and not raw_path.is_symlink() and raw_path.stat().st_size == size and sha256(raw_path) == digest,
                "raw source cache commitment differs: " + label)
    header = read_json(args.source_root / "S4-header-inspection-v1/headers-before-values.json")
    require(header.get("archive_sha256") == RAW["S4"][3] and len(header.get("sheets", [])) == 1,
            "S4 header source differs")
    sheet = header["sheets"][0]
    actual_headers = {re.sub(r"\d+$", "", c["coordinate"]): c.get("text")
                      for r in sheet["schema_rows"] if r["row_1based"] == 2 for c in r["cells"]}
    require(sheet["sheet"] == "S4.MPRA results" and all(actual_headers.get(k) == v for k, v in HEADERS.items()),
            "exact S4 text headers differ")
    missing_path = args.source_root / "missing-token-evidence.json"
    missing = read_json(missing_path)
    require(missing["missing_policy"]["allowed_tokens"] == [], "unfrozen missing-token policy")
    with gzip.open(args.source_root / "S4-ID-metadata.json.gz", "rt") as f:
        id_metadata = decode(f.read())
    id_rows = id_metadata["rows"]
    source_ids = {r["id"] for r in id_rows}
    require(id_metadata["source_sha256"] == RAW["S4"][3]
            and id_metadata["numeric_or_outcome_column_values_read"] == 0
            and len(source_ids) == len(id_rows) == id_metadata["unique_ID_count"]
            and id_metadata["duplicate_IDs"] == {}, "invalid authoritative S4 ID metadata")
    control_ledger = read_json(args.source_root / "S30-technical-control-matches.json")
    controls = {r["id"] for r in control_ledger["matched_external_variants"]}
    require(len(controls) == control_ledger["matched_external_variant_count"] == 41, "technical-control inventory differs")
    all_pairs_path = args.source_root / "design-mapping-v2/all-reference-validated-pairs.jsonl.gz"
    require(inputs["new_pairs"]["sha256"] == sha256(all_pairs_path)
            == control_ledger["all_reference_validated_pairs_sha256"], "complete new DNA source binding differs")
    components = read_json(args.groups / "components.json")["components"]
    pairs = list(read_jsonl(args.groups / "eligible-pairs.jsonl.gz"))
    cohort, assignments, component_ledger = select_cohort(components, pairs, source_ids, controls)
    pair_map = {r["id"]: r for r in pairs}
    seen = set()
    for row in read_jsonl(all_pairs_path):
        identity = row["id"]
        require(identity not in seen and identity in assignments, "original mapped DNA membership differs")
        seen.add(identity)
        if identity in pair_map:
            require(all(row[key] == pair_map[identity][key] for key in ["sequence_ref", "sequence_alt", "gc_preserving", "hg38"]),
                    "primary retained DNA differs from original design")
    require(seen == assignments.keys() and len(seen) == parameters["new_pair_count_expected"],
            "primary components do not cover complete mapped population")
    require(bool(cohort), "metadata cohort is empty; cannot seal a usable observation binding")
    independence = {"status": "not_attached_at_seal_creation; separate comparison required"}
    if args.independent_comparison is not None:
        require(args.independent_groups is not None and args.expected_independent_comparison_sha256 is not None,
                "independent comparison requires its exact SHA256 and run directory")
        require(sha256(args.independent_comparison) == args.expected_independent_comparison_sha256,
                "independent comparison SHA256 mismatch")
        comparison = read_json(args.independent_comparison)
        independent_manifest_sha = comparison["bindings"]["independent_manifest_sha256"]
        verify_inventory(args.independent_groups, "manifest.json", independent_manifest_sha, exact=True)
        check_comparison(comparison, args.expected_group_manifest_sha256, independent_manifest_sha)
        independent_inputs = read_json(args.independent_groups / "bindings.json")
        require(independent_inputs["new_DNA"]["sha256"] == sha256(all_pairs_path)
                and independent_inputs["prior_inventory"]["manifest_sha256"] == prior["manifest_sha256"]
                and independent_inputs["protocol"]["sha256"] == PROTOCOL_SHA
                and independent_inputs["technical_controls"]["sha256"] == sha256(args.source_root / "S30-technical-control-matches.json")
                and set(independent_inputs["technical_controls"]["ids"]) == controls,
                "independent grouping source/protocol/control bindings differ")
        independence = {"status": "matched", "comparison_path": str(args.independent_comparison),
            "comparison_sha256": args.expected_independent_comparison_sha256,
            "primary_manifest_sha256": args.expected_group_manifest_sha256,
            "independent_manifest_sha256": independent_manifest_sha,
            "checks_passed": comparison["checks_passed"], "failed_checks": [],
            "scope": comparison["scope"], "new_numeric_outcome_access": False}
    else:
        require(args.independent_groups is None and args.expected_independent_comparison_sha256 is None,
                "orphan independent-comparison arguments")
    args.out.mkdir(parents=True)
    shutil.copyfile(__file__, args.out / "implementation.py")
    write_jsonl(args.out / "cohort.jsonl.gz", cohort)
    write_jsonl(args.out / "component-exclusion-ledger.jsonl.gz", component_ledger)
    ledger_seen = set()

    def full_ledger():
        for row in read_jsonl(args.groups / "variant-ledger.jsonl.gz"):
            identity = row["id"]
            require(identity not in ledger_seen and identity in assignments, "duplicate/unknown full-group ledger ID")
            ledger_seen.add(identity)
            component = assignments[identity]
            require(row["component"] == component["component"]
                    and row["component_size"] == component["original_size"]
                    and row["retained_after_prior_exposure"] is component["retained_after_prior_exposure"],
                    "full-group ledger changes component or gate")
            reasons = list(component["component_exclusion_reasons"])
            if identity not in source_ids:
                reasons.append("absent_from_authoritative_S4_ID_membership")
            yield {"id": identity, "component": row["component"], "locus": row["locus"],
                   "original_component_size": row["component_size"], "S4_ID_present": identity in source_ids,
                   "admitted": not reasons, "exclusion_reasons": sorted(set(reasons))}

    write_jsonl(args.out / "variant-exclusion-ledger.jsonl.gz", full_ledger())
    require(ledger_seen == assignments.keys(), "incomplete full-group ledger")
    context_components = {r["component"] for r in cohort}
    binding = {"schema_version": 1, "protocol_id": "variant-transfer-v1", "protocol_sha256": PROTOCOL_SHA,
        "sealed_at_utc": datetime.now(timezone.utc).isoformat(), "model_artifacts_sha256": MODEL_SHA,
        "source_cohort_gate": {"passed": True, "scope": "primary metadata-only source/cohort gate",
            "whole_component_S30_exclusion": True, "S4_ID_membership_verified": True,
            "independent_grouping": independence,
            "numeric_effect_type_coverage_status": "not_evaluated_before_second_binding_commit"},
        "source": {"path": source_manifest["source_S4_raw_cache"], "sha256": RAW["S4"][3], "bytes": RAW["S4"][2]},
        "sheet": {"name": "S4.MPRA results", "header_row_1based": 2, "headers": HEADERS},
        "effects": {"HepG2": "AN", "K562": "AO"},
        "missing_policy": {"allowed_tokens": [], "evidence": str(missing_path), "evidence_sha256": sha256(missing_path)},
        "cohort": {"path": str(args.out / "cohort.jsonl.gz"), "format": "jsonl.gz", "variant_count": len(cohort),
                   "component_count": len(context_components), "sha256": sha256(args.out / "cohort.jsonl.gz")},
        "upstream_source_QC": {"review": str(args.source_root / "source-semantics-review.json"),
            "review_sha256": sha256(args.source_root / "source-semantics-review.json"),
            "policy": "Use unchanged authoritative source-preprocessed unique S4 results; no guessed count rethreshold or alternate library selection."},
        "inputs": {"source_manifest_sha256": SOURCE_MANIFEST_SHA,
            "primary_group_manifest_sha256": args.expected_group_manifest_sha256,
            "prior_inventory_manifest_sha256": prior["manifest_sha256"],
            "S30_match_ledger_sha256": sha256(args.source_root / "S30-technical-control-matches.json"),
            "S4_ID_metadata_sha256": sha256(args.source_root / "S4-ID-metadata.json.gz")},
        "metadata_counts": {"complete_external_population": len(assignments), "global_components": len(components),
            "S30_control_variant_ids": len(controls), "S30_control_components": sum(bool(r["S30_control_variant_ids"]) for r in component_ledger),
            "final_variants": len(cohort), "final_components": len(context_components)},
        "metadata_minimum_support": {"at_least100_variants_and10_components": len(cohort) >= 100 and len(context_components) >= 10,
            "limit": "This is a pre-outcome upper bound; final ordinary-finite numeric coverage and inferential support remain unverified."},
        "no_regrouping_after_exclusion": True, "S4_numeric_effect_values_read": False}
    (args.out / "binding.json").write_bytes(encode(binding))
    finish_manifest(args.out)
    print(json.dumps({"binding_sha256": sha256(args.out / "binding.json"), "manifest_sha256": sha256(args.out / "manifest.json"),
                      "variants": len(cohort), "components": len(context_components)}, sort_keys=True))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-root", type=Path, default=ROOT / "research/variant-transfer/source")
    p.add_argument("--groups", type=Path, required=True)
    p.add_argument("--expected-group-manifest-sha256", required=True)
    p.add_argument("--model-root", type=Path, default=ROOT / "research/results/context-activity-v3/core")
    p.add_argument("--protocol", type=Path, default=ROOT / "experiments/variant-transfer-v1.json")
    p.add_argument("--independent-comparison", type=Path)
    p.add_argument("--expected-independent-comparison-sha256")
    p.add_argument("--independent-groups", type=Path)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        p.error("sealing output already exists; choose a new attempt directory")
    try:
        seal(args)
    except Exception as error:
        if not args.out.exists():
            args.out.mkdir(parents=True)
        else:
            marker = args.out / "implementation.py"
            require(marker.is_file() and sha256(marker) == sha256(__file__),
                    "refusing to alter an output directory not created by this sealer")
        (args.out / "failure.json").write_bytes(encode({"schema_version": 1, "error_type": type(error).__name__,
            "error": str(error), "source_cohort_gate": {"passed": False}, "S4_numeric_effect_values_read": False}))
        finish_manifest(args.out)
        raise


if __name__ == "__main__":
    main()
