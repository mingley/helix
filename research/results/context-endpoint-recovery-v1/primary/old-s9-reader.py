#!/usr/bin/env python3
"""Offline, frozen post-hoc arithmetic audit; never changes prediction targets."""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
from io import BytesIO
import json
import math
from pathlib import Path
import platform
import unittest
import xml.etree.ElementTree as ET
from zipfile import ZipFile

import numpy as np


AUDIT_SHA256 = "8f092c7ff075d7f455791e5a72f189781e611e3a6bda69b0fd9d8a7947990eae"
PARENT_SHA256 = "8f0576845fecb5cfcf9bfa354fab446d849244448a3489905a6ce33bdeb21cec"
CONTEXTS = ("HepG2", "K562", "WTC11")
HEADERS = ("Fold", "Sequence ID") + tuple(
    f"{context}, Element Specificity Score" for context in CONTEXTS
)
SHEET = "Fold assignment, joint library"
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
SEMANTIC_STATUS = (
    "unresolved: published S9 columns fail the unweighted zero-sum formula; "
    "a fitted weighted plane cannot establish historical source preprocessing"
)


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def checked_bytes(path: Path, expected: str) -> bytes:
    raw = path.read_bytes()
    if sha256(raw) != expected:
        raise ValueError(f"SHA256 mismatch: {path}")
    return raw


def source_rows(raw: bytes) -> tuple[list[dict], dict]:
    """Read exact frozen S9 columns and retain raw XML decimal lexemes."""
    with ZipFile(BytesIO(raw)) as archive:
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            for item in ET.fromstring(archive.read("xl/sharedStrings.xml")).findall(
                "m:si", NS
            ):
                strings.append("".join(t.text or "" for t in item.findall(".//m:t", NS)))
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheets = workbook.findall("m:sheets/m:sheet", NS)
        matches = [s for s in sheets if s.attrib["name"] == SHEET]
        if len(matches) != 1:
            raise ValueError("Exact source worksheet missing or duplicated")
        relation = matches[0].attrib[f"{{{REL_NS}}}id"]
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target = next(r.attrib["Target"] for r in relationships if r.attrib["Id"] == relation)
        target = target.lstrip("/") if target.startswith("/") else "xl/" + target
        parsed = ET.fromstring(archive.read(target))
        xml_rows = parsed.findall("m:sheetData/m:row", NS)
        rows = []
        blank_ids = 0
        headers = None
        for row in xml_rows:
            number = int(row.attrib["r"])
            cells = {}
            for cell in row.findall("m:c", NS):
                column = cell.attrib["r"].rstrip("0123456789")
                if column not in ("A", "B", "C", "D", "E"):
                    continue
                value = cell.find("m:v", NS)
                cell_type = cell.attrib.get("t", "n")
                if cell_type == "s":
                    text = strings[int(value.text)] if value is not None else None
                elif cell_type == "inlineStr":
                    text = "".join(t.text or "" for t in cell.findall(".//m:t", NS))
                else:
                    text = value.text if value is not None else None
                cells[column] = (text, cell_type, cell.find("m:f", NS) is not None)
            if number == 1:
                headers = tuple(cells.get(c, (None,))[0] for c in "ABCDE")
                continue
            identifier = cells.get("B", (None,))[0]
            if identifier is None or identifier == "":
                blank_ids += 1
                continue
            decimal_values = []
            for column in "CDE":
                text, cell_type, formula = cells.get(column, (None, None, False))
                if text is None or cell_type != "n" or formula:
                    raise ValueError(f"Non-numeric or formula ESS source cell at row {number}")
                value = Decimal(text)
                if not value.is_finite() or not math.isfinite(float(value)):
                    raise ValueError(f"Non-finite ESS source cell at row {number}")
                decimal_values.append(text)
            rows.append({"id": identifier, "source_row": number, "raw": decimal_values})
    if headers != HEADERS:
        raise ValueError("Exact frozen S9 headers differ")
    if len(rows) != 55338 or len({r["id"] for r in rows}) != len(rows):
        raise ValueError("Expected exactly 55338 unique finite source IDs")
    rows.sort(key=lambda r: r["id"])
    return rows, {
        "worksheet": SHEET,
        "headers": list(headers),
        "xml_row_elements": len(xml_rows),
        "blank_id_rows_after_header": blank_ids,
        "finite_unique_source_ids": len(rows),
        "dimension": parsed.find("m:dimension", NS).attrib.get("ref"),
    }


def matrix(rows: list[dict]) -> np.ndarray:
    return np.asarray([[float(Decimal(v)) for v in row["raw"]] for row in rows], dtype=np.float64)


def residual_summary(values: np.ndarray) -> dict:
    return {
        "mean": float(np.mean(values)),
        "rms": float(np.sqrt(np.mean(values * values))),
        "maximum_absolute": float(np.max(np.abs(values))),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
    }


def singular_analysis(values: np.ndarray, centered: bool) -> dict:
    means = np.mean(values, axis=0)
    analyzed = values - means if centered else values
    _, singular_values, vt = np.linalg.svd(analyzed, full_matrices=False)
    threshold = float(singular_values[0] * 1e-10)
    rank = int(np.sum(singular_values > threshold))
    unit = vt[-1].copy()
    total = float(np.sum(unit))
    if total < 0:
        unit *= -1
        total *= -1
    if total != 0:
        weights = unit * (len(CONTEXTS) / total)
        normalization = "mean one; sum positive"
    else:
        first = next((x for x in unit if x != 0), 0)
        if first < 0:
            unit *= -1
        weights = unit.copy()
        normalization = "sum exactly zero; retained unit norm with first nonzero positive"
    fitted_residual = analyzed @ weights
    raw_residual = values @ weights
    return {
        "matrix": "column_centered" if centered else "raw",
        "singular_values_descending": singular_values.tolist(),
        "relative_rank_cutoff": 1e-10,
        "absolute_rank_cutoff": threshold,
        "numerical_rank": rank,
        "smallest_right_singular_vector_unit_norm": unit.tolist(),
        "weight_normalization": normalization,
        "weights_by_context": dict(zip(CONTEXTS, weights.tolist())),
        "all_weights_strictly_positive": bool(np.all(weights > 0)),
        "minimum_weight": float(np.min(weights)),
        "maximum_weight": float(np.max(weights)),
        "analyzed_matrix_weighted_residual": residual_summary(fitted_residual),
        "raw_matrix_weighted_residual": residual_summary(raw_residual),
        "column_mean_dot_weights": float(means @ weights),
        "interpretation": "Forensic fitted linear constraint only; no identified historical transform.",
    }


def exact_witness(rows: list[dict]) -> dict:
    with localcontext() as ctx:
        ctx.prec = 80
        sums = [sum((Decimal(v) for v in row["raw"]), Decimal(0)) for row in rows]
        position = max(range(len(rows)), key=lambda i: abs(sums[i]))
        row = rows[position]
        return {
            "selection": "maximum absolute exact-decimal row sum, ID order breaks ties",
            "id": row["id"],
            "source_row": row["source_row"],
            "raw_xml_decimal_values": dict(zip(CONTEXTS, row["raw"])),
            "exact_decimal_sum": str(sums[position]),
            "exact_decimal_absolute_sum": str(abs(sums[position])),
            "binary_float_fsum": math.fsum(float(v) for v in row["raw"]),
        }


def cohort_report(rows: list[dict]) -> dict:
    if len(rows) < 2:
        raise ValueError("At least two rows required")
    values = matrix(rows)
    means = np.mean(values, axis=0)
    population_sd = np.std(values, axis=0, ddof=0)
    sample_sd = np.std(values, axis=0, ddof=1)
    moments = {
        c: {
            "mean": float(means[i]),
            "population_sd_ddof0": float(population_sd[i]),
            "sample_sd_ddof1": float(sample_sd[i]),
            "minimum": float(np.min(values[:, i])),
            "maximum": float(np.max(values[:, i])),
        }
        for i, c in enumerate(CONTEXTS)
    }
    row_sums = np.sum(values, axis=1)
    abs_sums = np.abs(row_sums)
    raw = singular_analysis(values, False)
    centered = singular_analysis(values, True)
    witness = exact_witness(rows)
    witness_values = np.asarray([float(witness["raw_xml_decimal_values"][c]) for c in CONTEXTS])
    for analysis in (raw, centered):
        weights = np.asarray([analysis["weights_by_context"][c] for c in CONTEXTS])
        analysis["exact_witness_raw_weighted_residual_float"] = float(witness_values @ weights)
    return {
        "element_count": len(rows),
        "moments_by_context": moments,
        "row_sum": {
            "arithmetic": "float64 context-order sum; exact-decimal maximum witness also reported",
            **residual_summary(row_sums),
            "absolute_type7_quantiles": {
                str(q): float(np.quantile(abs_sums, q, method="linear"))
                for q in (0.0, 0.5, 0.95, 1.0)
            },
            "absolute_exceedance_counts": {
                str(t): int(np.sum(abs_sums > t)) for t in (1e-12, 1e-8, 1e-6)
            },
        },
        "raw_svd": raw,
        "column_centered_svd": centered,
        "fixed_source_standardization_checks": {
            "absolute_tolerance": 1e-8,
            "all_column_means_near_zero": bool(np.all(np.abs(means) <= 1e-8)),
            "all_sample_sd_near_one": bool(np.all(np.abs(sample_sd - 1) <= 1e-8)),
            "raw_smallest_vector_weights_all_positive": raw["all_weights_strictly_positive"],
            "centered_smallest_vector_weights_all_positive": centered["all_weights_strictly_positive"],
            "historical_normalization_verified": False,
        },
        "exact_decimal_row_sum_witness": witness,
    }


def run(protocol_path: Path, out: Path) -> dict:
    if out.exists():
        raise ValueError("Refusing existing output directory; preserve previous attempts")
    protocol_raw = checked_bytes(protocol_path, AUDIT_SHA256)
    protocol = json.loads(protocol_raw)
    if protocol["parent_protocol_sha256"] != PARENT_SHA256:
        raise ValueError("Unexpected parent protocol commitment")
    checked_bytes(Path("experiments/context-specificity-v2.json"), PARENT_SHA256)
    inputs = protocol["inputs"]
    workbook_raw = checked_bytes(Path(inputs["source_workbook"]), inputs["source_sha256"])
    panel_raw = checked_bytes(Path(inputs["panel"]), inputs["panel_sha256"])
    panel = json.loads(panel_raw)
    if panel["contexts"] != list(CONTEXTS) or len(panel["records"]) != 55276:
        raise ValueError("Unexpected frozen admitted panel dimensions")
    rows, schema = source_rows(workbook_raw)
    lookup = {r["id"]: r for r in rows}
    ids = [r["id"] for r in panel["records"]]
    if len(set(ids)) != len(ids) or any(i not in lookup for i in ids):
        raise ValueError("Panel IDs are not a unique subset of source IDs")
    maximum_source_panel_difference = 0.0
    for record in panel["records"]:
        for context, text in zip(CONTEXTS, lookup[record["id"]]["raw"]):
            difference = abs(float(record["observations"][context]) - float(Decimal(text)))
            maximum_source_panel_difference = max(maximum_source_panel_difference, difference)
    if maximum_source_panel_difference != 0:
        raise ValueError("Frozen panel targets differ from exact source-cell float values")
    admitted = sorted((lookup[i] for i in ids), key=lambda r: r["id"])
    implementation = Path(__file__).read_bytes()
    summary = {
        "protocol_id": protocol["protocol_id"],
        "diagnostic_status": "completed post-hoc arithmetic diagnostic; no target/model changes",
        "scientific_semantic_status": SEMANTIC_STATUS,
        "historical_normalization_verified": False,
        "hashes": {
            "audit_protocol_sha256": sha256(protocol_raw),
            "parent_protocol_sha256": PARENT_SHA256,
            "source_workbook_sha256": sha256(workbook_raw),
            "panel_sha256": sha256(panel_raw),
            "implementation_sha256": sha256(implementation),
        },
        "runtime": {"python": platform.python_version(), "numpy": np.__version__},
        "source_schema": schema,
        "source_panel_maximum_absolute_difference": maximum_source_panel_difference,
        "all55338finite_exact_source_score_IDs": cohort_report(rows),
        "all55276frozen_reference_valid_panel_IDs": cohort_report(admitted),
        "prohibitions_observed": [
            "No target modification, projection, rescaling, filtering or predictive refitting.",
            "No fitted source weights enter the frozen benchmark.",
            "No historical normalization identification or biological discovery from a fitted plane.",
        ],
    }
    encoded = (json.dumps(summary, indent=2, allow_nan=False) + "\n").encode()
    out.mkdir(parents=True)
    (out / "summary.json").write_bytes(encoded)
    (out / "protocol.json").write_bytes(protocol_raw)
    (out / "implementation.py").write_bytes(implementation)
    artifacts = {
        p.name: {"bytes": p.stat().st_size, "sha256": sha256(p.read_bytes())}
        for p in sorted(out.iterdir())
    }
    (out / "artifacts.json").write_text(json.dumps(artifacts, indent=2) + "\n")
    return summary


class ArithmeticTests(unittest.TestCase):
    def test_decimal_witness_retains_cancellation(self):
        witness = exact_witness([{"id": "toy", "source_row": 2, "raw": ["0.1", "0.2", "-0.3"]}])
        self.assertEqual(Decimal(witness["exact_decimal_sum"]), Decimal(0))

    def test_rescaled_centered_plane_is_not_unweighted_zero_sum(self):
        base = np.array([[-1., 0., 1.], [1., -2., 1.], [2., 1., -3.], [-2., 1., 1.]])
        scaled = base / base.std(axis=0, ddof=1)
        result = singular_analysis(scaled, False)
        self.assertEqual(result["numerical_rank"], 2)
        self.assertTrue(result["all_weights_strictly_positive"])
        self.assertLess(result["analyzed_matrix_weighted_residual"]["maximum_absolute"], 1e-14)
        self.assertGreater(np.max(np.abs(scaled.sum(axis=1))), 0.1)

    def test_full_rank_is_not_a_rank_two_plane(self):
        self.assertEqual(singular_analysis(np.eye(3), False)["numerical_rank"], 3)

    def test_empty_cohort_rejected(self):
        with self.assertRaises(ValueError):
            cohort_report([])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=Path("experiments/context-source-arithmetic-audit-v1.json"))
    parser.add_argument("--out", type=Path, default=Path("runs/context-source-arithmetic-audit-v1"))
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(ArithmeticTests))
        if not result.wasSuccessful():
            raise SystemExit(1)
        return
    summary = run(args.protocol, args.out)
    print(json.dumps({"out": str(args.out), "summary_sha256": sha256((args.out / "summary.json").read_bytes()), "semantic_status": summary["scientific_semantic_status"]}))


if __name__ == "__main__":
    main()
