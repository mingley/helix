#!/usr/bin/env python3
"""Six frozen post-hoc ESS-cohort routes; never changes empirical targets."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import platform
import shutil
import sys
import unittest
import xml.etree.ElementTree as ET
from zipfile import ZipFile

import numpy as np


PROTOCOL_SHA = "a0e3f521d1086fca0633fc4c473a3655aaf2fa9d4a4818b0429a735f46d42372"
HELPER_SHA = "2472660613a5f63f41394eb91aaf99b8e5606fddeca0b23e2ee30c372c75a4f4"
READER_SHA = "9a9779f80ea3002eab0ba904530cf1ec5c1bc3d7bb2f96155018d0369a3dbf42"
CONTEXTS = ("HepG2", "K562", "WTC11")
S12_IDS_SHA = "b85af55d0b1f80549d7b70b3ae7f22538215e31627efec0dab76c7272efe9075"
ADMITTED_IDS_SHA = "ba9aa3801e051776e3c899f890a58ce2dfeb18fbb372f6b866b63888635266d7"
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def checked(path, expected, size=None):
    path = Path(path)
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected or (size is not None and len(raw) != size):
        raise ValueError(f"Frozen bytes/SHA256 differ: {path}")
    return raw


def load_module(name, path, expected):
    checked(path, expected)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def primary_path(spec, override):
    return Path(override) / Path(spec["path"]).name if override is not None else Path(spec["path"])


def load_base(helper, name, raw):
    value = json.loads(gzip.decompress(raw))
    if value["name"] != name or value["contexts"] != list(CONTEXTS):
        raise ValueError("Frozen activity-base identity differs")
    return helper.Endpoint(name, value["ids"], np.asarray(value["values"], dtype=float),
                           value["units"], blocked=value.get("blocked_reason"))


def normalize_on_target(helper, base, target_ids, ddof):
    """Require every target row before invoking the archived fixed arithmetic."""
    values = base.subset(target_ids)
    missing = [i for i, row in zip(target_ids, values) if not np.isfinite(row).all()]
    report = {"id": f"{base.name}__exact_ESS_cohort__ddof{ddof}", "activity_base": base.name,
              "pre_zscore_ddof": ddof, "estimable": False,
              "coverage": {"required_ids": len(target_ids), "ids_sha256": helper.ids_hash(target_ids),
                           "missing_or_nonfinite_ids": missing,
                           "missing_contexts": {i: [c for c, v in zip(CONTEXTS, row) if not np.isfinite(v)]
                                                for i, row in zip(target_ids, values) if not np.isfinite(row).all()}}}
    if base.blocked or missing:
        report["nonestimable_reason"] = base.blocked or "required entire ESS cohort not covered; no intersection rescue"
        return report, None
    scoped = helper.Endpoint(base.name, list(target_ids), values, base.units)
    inner, result = helper.pipeline(scoped, helper.POPULATIONS[0], ddof, "none")
    report["normalization_record"] = inner
    report["estimable"] = inner["estimable"]
    if result is None:
        report["nonestimable_reason"] = inner["nonestimable_reason"]
    return report, result


def design_categories(path):
    """Read only exact S10 name/category fields, never choose a cohort from them."""
    with ZipFile(path) as archive:
        strings = ["".join(t.text or "" for t in si.findall(".//m:t", NS))
                   for si in ET.fromstring(archive.read("xl/sharedStrings.xml")).findall("m:si", NS)]
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheet = next(s for s in workbook.findall("m:sheets/m:sheet", NS) if s.attrib["name"] == "joint library")
        rid = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target = next(r.attrib["Target"] for r in relationships if r.attrib["Id"] == rid)
        target = target.lstrip("/") if target.startswith("/") else "xl/" + target
        categories, headers = {}, {}
        for row in ET.fromstring(archive.read(target)).findall("m:sheetData/m:row", NS):
            number = int(row.attrib["r"])
            cells = {}
            for cell in row.findall("m:c", NS):
                column = cell.attrib["r"].rstrip("0123456789")
                if column not in ("A", "B"):
                    continue
                v = cell.find("m:v", NS)
                if cell.attrib.get("t") == "s":
                    cells[column] = strings[int(v.text)] if v is not None else None
                elif cell.attrib.get("t") == "inlineStr":
                    cells[column] = "".join(t.text or "" for t in cell.findall(".//m:t", NS))
                else:
                    cells[column] = v.text if v is not None else None
            if number == 2:
                headers = cells
            elif number >= 3 and cells.get("A"):
                identifier = cells["A"]
                if identifier in categories:
                    raise ValueError("Duplicate S10 source category key")
                categories[identifier] = cells.get("B")
    if headers != {"A": "name", "B": "category"}:
        raise ValueError("Exact S10 category headers differ")
    return categories


def metadata_audit(helper, bases, ess_ids, s12_ids, admitted_ids, categories):
    sets = {"S11_ESS": set(ess_ids), "S12": set(s12_ids), "old_admitted": set(admitted_ids)}
    for base in bases:
        sets[base.name + "__all_ids"] = set(base.ids)
        sets[base.name + "__complete_ids"] = {i for i, row in zip(base.ids, base.values) if np.isfinite(row).all()}

    def describe(ids):
        ordered = sorted(ids)
        counts = Counter(categories.get(i) or "__UNMAPPED_CATEGORY__" for i in ordered)
        return {"count": len(ordered), "ids_sha256": helper.ids_hash(ordered), "identifiers": ordered,
                "existing_source_category_counts": dict(sorted(counts.items()))}

    return {"cohorts": {name: describe(ids) for name, ids in sorted(sets.items())},
            "pairwise_sets": [{"left": left, "right": right,
                               "intersection": describe(sets[left] & sets[right]),
                               "left_minus_right": describe(sets[left] - sets[right]),
                               "right_minus_left": describe(sets[right] - sets[left])}
                              for left, right in itertools.combinations(sorted(sets), 2)],
            "policy": "Categories are metadata counts only; they never select or alter a normalization population."}


def run(args, helper, reader, primary_root):
    if args.out.exists():
        raise FileExistsError("Refusing existing outputs; preserve every earlier attempt")
    protocol_raw = checked(args.protocol, PROTOCOL_SHA)
    protocol = json.loads(protocol_raw)
    if protocol["status"] != "frozen":
        raise ValueError("Protocol is not frozen")
    checked("experiments/context-endpoint-recovery-v1.json", protocol["parent_protocol_sha256"])
    inputs, receipts, loaded = protocol["inputs"], {}, {}
    for key in ("primary_summary", "source_values", "schema_binding"):
        spec = inputs[key]
        path = primary_path(spec, args.primary_root)
        loaded[key] = checked(path, spec["sha256"], spec.get("bytes"))
        receipts[key] = {"transport_path": str(path), "bytes": len(loaded[key]), "sha256": spec["sha256"]}
    source = json.loads(gzip.decompress(loaded["source_values"]))
    binding = json.loads(loaded["schema_binding"])
    official = inputs["official_jointdata"]
    workbook = args.source_workbook or Path(official["path"])
    checked(workbook, official["sha256"], official["bytes"])
    if source["source"]["sha256"] != official["sha256"] or source["schema_binding_sha256"] != helper.sha(primary_path(inputs["schema_binding"], args.primary_root)):
        raise ValueError("Source export ancestry differs")
    target = helper.direct_endpoint("unchanged_S11_ESS", source["tables"], binding["roles"]["published_ess"], "dimensionless published ESS")
    if target.blocked or not np.isfinite(target.values).all():
        raise ValueError("Published target ESS must be unique and finite in all contexts")
    population = protocol["normalization_population"]
    if len(target.ids) != population["expected_ids"] or helper.ids_hash(target.ids) != population["expected_ids_sha256"]:
        raise ValueError("Exact frozen S11 ESS ID population differs")
    bases = []
    for spec in inputs["activity_bases"]:
        path = primary_path(spec, args.primary_root)
        raw = checked(path, spec["sha256"], spec["bytes"])
        bases.append(load_base(helper, spec["id"], raw))
        receipts[spec["id"]] = {"transport_path": str(path), "bytes": len(raw), "sha256": spec["sha256"]}
    s12_raw = checked(inputs["S12"]["path"], inputs["S12"]["sha256"])
    rows, s12_ledger = reader.source_rows(s12_raw)
    s12 = helper.Endpoint("unchanged_S12", [r["id"] for r in rows],
                          np.asarray([[float(v) for v in r["raw"]] for r in rows]),
                          "dimensionless published source score", raw={r["id"]: r["raw"] for r in rows})
    if len(s12.ids) != 55338 or helper.ids_hash(s12.ids) != S12_IDS_SHA:
        raise ValueError("Exact frozen S12 ID population differs")
    admitted_spec = inputs["old_admitted_panel"]
    panel = json.loads(checked(admitted_spec["path"], admitted_spec["sha256"]))
    admitted = [r["id"] for r in panel["records"]]
    if len(admitted) != 55276 or helper.ids_hash(admitted) != ADMITTED_IDS_SHA:
        raise ValueError("Fixed admitted population differs")
    category_spec = inputs["design_category_source"]
    checked(category_spec["path"], category_spec["sha256"])
    categories = design_categories(category_spec["path"])
    metadata = metadata_audit(helper, bases, target.ids, s12.ids, admitted, categories)
    args.out.mkdir(parents=True)
    helper.write_json(args.out / "protocol.json", protocol)
    # Preserve exact protocol bytes rather than the display serializer's choices.
    (args.out / "protocol.json").write_bytes(protocol_raw)
    helper.write_gzip(args.out / "id-category-ledger.json.gz", metadata)
    helper.write_json(args.out / "populations-before-arithmetic.json", {
        "S11_ESS": {"count": len(target.ids), "ids_sha256": helper.ids_hash(target.ids)},
        "S12": {"count": len(s12.ids), "ids_sha256": helper.ids_hash(s12.ids)},
        "old_admitted": {"count": len(admitted), "ids_sha256": helper.ids_hash(admitted)},
        "normalization_will_not_use_categories": True})
    source_means = np.mean(target.values, axis=0)
    routes = []
    for base, ddof in itertools.product(bases, (0, 1)):
        record, result = normalize_on_target(helper, base, target.ids, ddof)
        if result is not None:
            record["comparison_to_unchanged_S11"] = helper.compare(result, target, s12.ids, admitted)
            record["structural_column_mean_audit"] = {
                "source_S11_means": dict(zip(CONTEXTS, source_means.tolist())),
                "reconstructed_means": dict(zip(CONTEXTS, np.mean(result.values, axis=0).tolist())),
                "source_values_forced_to_zero": False}
            file = record["id"] + ".json.gz"
            helper.save_endpoint(args.out / file, result)
            record["derived_values_file"] = file
            derivative, derivative_values = helper.code_hypothesis(result, s12, s12.ids, admitted)
            record["fixed_S12_derivative"] = derivative
            if derivative_values is not None:
                file = record["id"] + "__S12_sample_zscore.json.gz"
                helper.save_endpoint(args.out / file, derivative_values)
                derivative["derived_values_file"] = file
        else:
            record["fixed_S12_derivative"] = {"estimable": False, "nonestimable_reason": "parent ESS-cohort route is nonestimable"}
        routes.append(record)
    if len(routes) != 6 or len({(r["activity_base"], r["pre_zscore_ddof"]) for r in routes}) != 6:
        raise ValueError("Frozen six-route Cartesian grid differs")
    summary = {"protocol_id": protocol["protocol_id"], "protocol_sha256": PROTOCOL_SHA,
               "post_hoc": True, "route_count": 6, "estimable_routes": sum(r["estimable"] for r in routes),
               "normalization_population": {"count": len(target.ids), "ids_sha256": helper.ids_hash(target.ids)},
               "routes": routes, "source_S11_moments": helper.moments(target.values),
               "source_S11_row_sums": helper.row_sum_summary(target),
               "metadata_ledger": "id-category-ledger.json.gz", "S12_source_ledger": s12_ledger,
               "prior_36_results_preserved": True, "source_values_changed": False,
               "predictive_models_fitted": False, "historical_pipeline_certified": False,
               "interpretation": "All six population hypotheses retained. Compatibility is not unique producer history; no replacement target or independent biology."}
    helper.write_json(args.out / "summary.json", summary)
    implementation = Path(__file__)
    inputs_report = {"protocol_sha256": PROTOCOL_SHA, "files": receipts,
                     "primary_root_transport_override": str(args.primary_root) if args.primary_root else None,
                     "official_workbook_transport_path": str(workbook),
                     "official_workbook_sha256": official["sha256"],
                     "archived_primary_helper_sha256": HELPER_SHA,
                     "archived_S12_reader_sha256": READER_SHA,
                     "implementation_sha256": helper.sha(implementation),
                     "runtime": {"python": platform.python_version(), "numpy": np.__version__},
                     "completed_at_utc": datetime.now(timezone.utc).isoformat()}
    helper.write_json(args.out / "inputs.json", inputs_report)
    shutil.copyfile(implementation, args.out / "implementation.py")
    shutil.copyfile(primary_root / "implementation.py", args.out / "primary-helper.py")
    shutil.copyfile(primary_root / "old-s9-reader.py", args.out / "S12-reader.py")
    files = {p.relative_to(args.out).as_posix(): {"bytes": p.stat().st_size, "sha256": helper.sha(p)}
             for p in sorted(args.out.rglob("*")) if p.is_file()}
    helper.write_json(args.out / "artifacts.json", {"protocol_sha256": PROTOCOL_SHA, "files": files})
    print(json.dumps({"output": str(args.out), "routes": 6, "estimable": summary["estimable_routes"],
                      "summary_sha256": helper.sha(args.out / "summary.json"), "targets_changed": False}))


def controls(helper):
    class Controls(unittest.TestCase):
        def endpoint(self, values, ids, name="toy"):
            return helper.Endpoint(name, ids, np.asarray(values, dtype=float), "toy")

        def test_target_population_excludes_extra_outlier(self):
            base = self.endpoint([[-1, 0, 1], [0, 1, -1], [1, -1, 0], [100, 100, 100]], list("abcd"))
            for ddof, factor in ((0, np.sqrt(1.5)), (1, 1)):
                report, value = normalize_on_target(helper, base, list("abc"), ddof)
                self.assertTrue(report["estimable"])
                np.testing.assert_allclose(value.values, base.values[:3] * factor, atol=1e-15)
                self.assertTrue(all(p["count"] == 3 for p in report["normalization_record"]["pre_parameters"]))

        def test_missing_target_blocks_both_routes(self):
            base = self.endpoint([[1, 2, 3], [2, 3, 4]], list("ab"))
            for ddof in (0, 1):
                report, value = normalize_on_target(helper, base, list("abc"), ddof)
                self.assertFalse(report["estimable"])
                self.assertIsNone(value)
                self.assertEqual(report["coverage"]["missing_or_nonfinite_ids"], ["c"])

        def test_nonfinite_context_is_not_partial_mean(self):
            base = self.endpoint([[-1, 0, 1], [0, np.nan, -1], [1, -1, 0]], list("abc"))
            report, value = normalize_on_target(helper, base, list("abc"), 0)
            self.assertFalse(report["estimable"])
            self.assertIsNone(value)
            self.assertEqual(report["coverage"]["missing_contexts"], {"b": ["K562"]})

        def test_decimal_constant_is_nonestimable(self):
            base = self.endpoint([[0.1, i, 6-i] for i in range(7)], [str(i) for i in range(7)])
            report, value = normalize_on_target(helper, base, base.ids, 0)
            self.assertFalse(report["estimable"])
            self.assertIsNone(value)

        def test_derivative_uses_entire_S12_subset(self):
            source = self.endpoint([[-2, 0, 2], [0, 3, -3], [2, -3, 1], [100, 100, 100]], list("abcd"))
            selected = source.values[:3]
            target = self.endpoint((selected-selected.mean(0))/selected.std(0, ddof=1), list("abc"))
            report, _ = helper.code_hypothesis(source, target, list("abc"), ["a"])
            self.assertEqual(report["normalization_population"]["count"], 3)
            self.assertTrue(report["fixed_1e-12_test_by_cohort"]["full_s12_source_id_set"]["matched"])

    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Controls))
    if not result.wasSuccessful():
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=Path("experiments/context-ess-cohort-audit-v1.json"))
    parser.add_argument("--primary-root", type=Path,
                        help="Explicit transport root for byte-identical frozen primary files, matched by basename")
    parser.add_argument("--source-workbook", type=Path, help="Explicit byte-identical official jointdata workbook path")
    parser.add_argument("--out", type=Path, default=Path("runs/context-ess-cohort-audit-v1"))
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    primary_root = args.primary_root or Path("runs/recovered-context-primary")
    helper = load_module("cohort_primary_helper", primary_root / "implementation.py", HELPER_SHA)
    reader = load_module("cohort_source_reader", primary_root / "old-s9-reader.py", READER_SHA)
    if args.self_test:
        controls(helper)
    else:
        run(args, helper, reader, primary_root)


if __name__ == "__main__":
    main()
