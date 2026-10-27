#!/usr/bin/env python3
"""Bind unchanged published activity to the pre-existing DNA-only v3 cohort."""
from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
from io import BytesIO
import json
import math
from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_SHA = "64e14037fb1f42ddc1ecaca1d23bb8096fabac1f1e892d105eaa20f525fbbb20"
RECOVERY_SHA = "00175e41cf9fa0cbf177ce8fc7ef143eb508aec8b11d058f7acb4458284dc868"
RAW_SHA = "9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e"
ENDPOINT = "published_joint_log2_rna_dna_activity"
CONTEXTS = ["HepG2", "K562", "WTC11"]
COLUMNS = {cell: f"{cell} [log2(rna/dna)]" for cell in CONTEXTS}


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            result.update(block)
    return result.hexdigest()


def load(path):
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as stream:
            return json.load(stream)
    return json.loads(path.read_bytes())


def packed(value):
    output = BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=output, mtime=0) as stream:
        stream.write((json.dumps(value, separators=(",", ":"), allow_nan=False) + "\n").encode())
    return output.getvalue()


def finite_source(value):
    if not isinstance(value, str):
        raise ValueError("required activity lacks a raw numeric XML literal")
    try:
        decimal = Decimal(value)
        converted = float(decimal)
    except (InvalidOperation, ValueError, OverflowError) as error:
        raise ValueError("invalid numeric activity source literal") from error
    if not decimal.is_finite() or not math.isfinite(converted):
        raise ValueError("required activity is missing or nonfinite")
    return converted


def activity_records(export, dna_ids):
    """No score selection, aliases, interpolation or outcomes from the old panel."""
    table = export["tables"]["all_cell_types_summary"]
    if table["id_header"] != "name":
        raise ValueError("the frozen literal source ID header differs")
    if not set(COLUMNS.values()).issubset(table["headers"].values()):
        raise ValueError("the frozen empirical activity headers differ")
    seen, records, all_ids = set(), {}, set()
    blank_ids = 0
    for row in table["rows"]:
        identity = row["id"]
        if identity is None or identity == "":
            blank_ids += 1
            continue
        if not isinstance(identity, str) or identity in all_ids:
            raise ValueError("non-string or duplicate source element ID")
        all_ids.add(identity)
        if identity not in dna_ids:
            continue
        seen.add(identity)
        records[identity] = {"id": identity, "observations": {
            cell: finite_source(row["values"].get(header))
            for cell, header in COLUMNS.items()}}
    if seen != dna_ids:
        raise ValueError(f"frozen activity cohort incomplete: {len(dna_ids-seen)} missing IDs")
    return [records[key] for key in sorted(records)], {
        "source_unique_nonempty_ids": len(all_ids), "source_blank_id_rows": blank_ids,
        "fixed_required_ids": len(dna_ids), "complete_required_ids": len(seen),
        "source_ids_outside_fixed_cohort": len(all_ids - dna_ids),
        "no_activity_based_selection": True, "no_imputation": True,
        "unchanged_published_values": True}


def build(args):
    protocol_path = ROOT / "experiments/context-activity-v3.json"
    if digest(protocol_path) != PROTOCOL_SHA:
        raise ValueError("frozen v3 protocol differs")
    protocol = load(protocol_path)
    export, binding = load(args.source_export), load(args.source_binding)
    source_receipt = load(args.source_receipt)
    if (source_receipt.get("protocol_sha256") != RECOVERY_SHA
            or source_receipt.get("schema_binding_sha256") != digest(args.source_binding)
            or source_receipt["export"]["sha256"] != digest(args.source_export)
            or source_receipt["ledger"]["sha256"] != digest(args.source_ledger)):
        raise ValueError("raw source export and numeric-type ledger are not jointly bound")
    if (export["protocol_sha256"] != RECOVERY_SHA
            or export["schema_binding_sha256"] != digest(args.source_binding)
            or export["source"]["sha256"] != RAW_SHA):
        raise ValueError("source export is not bound to frozen recovery/schema/raw source")
    if (not binding["semantic_binding_complete"]
            or binding["status"] != "bound_before_quantitative_inspection"
            or binding["protocol_sha256"] != RECOVERY_SHA
            or binding["source_workbook"]["sha256"] != RAW_SHA):
        raise ValueError("empirical activity schema is not bound")
    raw = ROOT / "datasets/context-endpoint-recovery/raw/jointdata.xlsx"
    if digest(raw) != RAW_SHA or raw.stat().st_size != 53557805:
        raise ValueError("official joint activity workbook differs")
    dna = ROOT / protocol["panel"]["path"]
    groups = ROOT / protocol["grouping"]["path"]
    if digest(dna) != protocol["panel"]["sha256"] or digest(groups) != protocol["grouping"]["sha256"]:
        raise ValueError("frozen DNA cohort or outcome-independent grouping differs")
    design = ROOT / "datasets/context-panel/raw/Table-S7-joint-design.xlsx"
    if digest(design) != "f9ec4c6c2d230b0ef507a39e5d5c73b3b36f90ccc5758c33b442c7561e20e7e5":
        raise ValueError("previously reference-validated design source differs")
    # The official archive validation/extraction receipts are source evidence;
    # their hashes are copied into the manifest, rather than inventing a source
    # signature from an ETag or a downstream project's digest.
    archive_validation = load(args.archive_validation)
    if (archive_validation.get("complete") is not True
            or archive_validation.get("status") != "archive_integrity_verified"
            or archive_validation.get("all_members_crc_verified") is not True
            or archive_validation.get("protocol_sha256") != RECOVERY_SHA):
        raise ValueError("complete official archive integrity validation is required")
    extraction = load(args.extraction_receipt)
    extracted = {item["literal_member_name"]: item for item in extraction["selected_members"]}
    design_member = "2023-03-03628-s5/2023-03-03628C-Table S10 - joint lib design.xlsx"
    activity_member = "2023-03-03628-s5/2023-03-03628C-Table S11 - jointdata.xlsx"
    if (extraction["protocol_sha256"] != RECOVERY_SHA
            or extraction["archive_sha256"] != archive_validation["archive_sha256"]
            or extracted[design_member]["sha256"] != digest(design)
            or extracted[activity_member]["sha256"] != RAW_SHA
            or not all(item["crc32_verified_by_zipfile"] for item in extracted.values())):
        raise ValueError("official design/activity extraction identities differ")
    dna_ids = {record["id"] for record in load(dna)["records"]}
    if len(dna_ids) != 55276:
        raise ValueError("fixed v3 DNA cohort has changed")
    records, coverage = activity_records(export, dna_ids)
    # A selected numeric string, formula cache, boolean or Excel error must not
    # become a measured float just because Decimal happens to parse its text.
    with gzip.open(args.source_ledger, "rt") as ledger:
        for line in ledger:
            item = json.loads(line)
            if (item.get("record_type") == "cell_exception"
                    and item.get("sheet") == "all_cell_types_summary"
                    and item.get("id") in dna_ids
                    and item.get("header") in COLUMNS.values()):
                raise ValueError("required activity is not an ordinary finite numeric source cell")
    recovery = load(args.recovery_summary)
    if (recovery.get("protocol_sha256") != RECOVERY_SHA
            or recovery.get("variant_count") != 36
            or recovery.get("flags", {}).get("predictive_models_fitted") is not False
            or recovery.get("flags", {}).get("empirical_raw_endpoint_binding") is not True
            or not recovery.get("aggregate_comparisons")):
        raise ValueError("record replicate aggregation and normalization review before fitting")
    recovery_root = args.recovery_summary.resolve().parent
    recovery_inputs = load(recovery_root / "inputs.json")
    recovery_artifacts = load(recovery_root / "artifacts.json")
    if (recovery_inputs.get("protocol_sha256") != RECOVERY_SHA
            or recovery_artifacts.get("protocol_sha256") != RECOVERY_SHA):
        raise ValueError("recovery input/result protocol binding differs")
    for name, expected in recovery_artifacts["files"].items():
        artifact = (recovery_root / name).resolve()
        if recovery_root not in artifact.parents or digest(artifact) != expected["sha256"]:
            raise ValueError("recovery artifact inventory is stale or incomplete")
    for path in (args.source_export, args.source_binding, args.source_ledger,
                 args.source_receipt, raw, dna, args.archive_validation):
        relative = path.resolve().relative_to(ROOT).as_posix()
        if recovery_inputs["files"][relative]["sha256"] != digest(path):
            raise ValueError("recovery arithmetic was run on different source inputs")
    activity = {"schema_version": 1, "protocol_sha256": PROTOCOL_SHA,
                "endpoint": ENDPOINT, "contexts": CONTEXTS, "records": records}
    out = Path(args.out)
    if args.verify:
        manifest = load(out / "manifest.json")
        for name, expected in manifest["files"].items():
            if Path(name).name != name or digest(out / name) != expected:
                raise ValueError("activity dataset artifact mismatch")
        if (out / "activity.json.gz").read_bytes() != packed(activity):
            raise ValueError("saved activity values differ from unchanged source export")
        print(json.dumps({"verified": True, "coverage": coverage}))
        return
    out.mkdir(parents=True, exist_ok=False)
    (out / "activity.json.gz").write_bytes(packed(activity))
    paths = {"source-binding.json": args.source_binding,
             "source-export.json.gz": args.source_export,
             "source-ledger.jsonl.gz": args.source_ledger,
             "recovery-summary.json": args.recovery_summary,
             "archive-validation.json": args.archive_validation,
             "extraction-receipt.json": args.extraction_receipt,
             "source-export-receipt.json": args.source_receipt,
             "recovery-inputs.json": recovery_root / "inputs.json",
             "recovery-artifacts.json": recovery_root / "artifacts.json"}
    for name, path in paths.items():
        shutil.copyfile(path, out / name)
    manifest = {"schema_version": 1, "protocol_sha256": PROTOCOL_SHA,
                "endpoint": ENDPOINT, "source": {
                    "raw_path": "../context-endpoint-recovery/raw/jointdata.xlsx",
                    "sha256": RAW_SHA, "bytes": raw.stat().st_size},
                "gates": {"archive_integrity_verified": True,
                          "empirical_activity_schema_bound": True,
                          "coverage_complete": True, "design_identity_verified": True,
                          "replicate_aggregation_status_recorded": True},
                "coverage": coverage, "aggregation_status": recovery,
                "schema_binding": binding,
                "files": {path.name: digest(path) for path in sorted(out.iterdir())},
                "builder_sha256": digest(__file__),
                "scope": "Unchanged published activity aggregate on the fixed prior DNA cohort; no ESS targets or fits reused."}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"manifest_sha256": digest(out / "manifest.json"), "coverage": coverage}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-export", required=True, type=Path)
    parser.add_argument("--source-binding", required=True, type=Path)
    parser.add_argument("--source-ledger", required=True, type=Path)
    parser.add_argument("--source-receipt", required=True, type=Path)
    parser.add_argument("--recovery-summary", required=True, type=Path)
    parser.add_argument("--archive-validation", required=True, type=Path)
    parser.add_argument("--extraction-receipt", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--verify", action="store_true")
    build(parser.parse_args())
