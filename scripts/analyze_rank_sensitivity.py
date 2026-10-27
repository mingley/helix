#!/usr/bin/env python3
"""Post-hoc numerical rank sensitivity; never a replacement for frozen v1."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import pathlib
import platform
import sys
from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np

PROTOCOL_SHA = "b968f0b3519b3890c36d2756ce97780f3a42effb8ee06b138e4cf395138b89ed"
PARENT_SHA = "c8e97fabef2afe8ccaa5a01a7ef7bb952f9c9f15857e5c540af112a66e78cb58"
LEVELS = ["raw_exact", "1e-12", "1e-10", "1e-8", "1e-6", "0.01"]
TOLERANCE = 1e-8


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def file_hash(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    path = pathlib.Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return json.load(handle)
    return json.loads(path.read_text())


def checked_files(root, inventory):
    root = pathlib.Path(root).resolve()
    for name, expected in inventory.items():
        path = (root / name).resolve()
        if path.parent != root or pathlib.Path(name).name != name:
            raise ValueError("artifact path escapes its directory")
        if file_hash(path) != expected:
            raise ValueError(f"artifact checksum mismatch: {name}")


def average_ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    first = 0
    while first < len(order):
        end = first + 1
        while end < len(order) and values[order[first]] == values[order[end]]:
            end += 1
        mean_rank = (first + end + 1) / 2.0
        for index in order[first:end]:
            ranks[index] = mean_rank
        first = end
    return ranks


def rank_correlation(observed, predicted):
    if len(observed) != len(predicted) or not observed:
        raise ValueError("rank inputs must be paired and nonempty")
    left, right = average_ranks(observed), average_ranks(predicted)
    if len(left) < 2 or min(left) == max(left) or min(right) == max(right):
        return None
    left_mean = math.fsum(left) / len(left)
    right_mean = math.fsum(right) / len(right)
    covariance = math.fsum((a - left_mean) * (b - right_mean)
                           for a, b in zip(left, right))
    aa = math.fsum((v - left_mean) ** 2 for v in left)
    bb = math.fsum((v - right_mean) ** 2 for v in right)
    return covariance / math.sqrt(aa * bb)


def quantized_integer(value, quantum):
    """Nearest decimal-grid integer; exact binary input and ties-to-even."""
    if not math.isfinite(value) or quantum <= 0:
        raise ValueError("quantization requires finite prediction and positive grid")
    ratio = Fraction.from_float(float(value)) / quantum
    lower = ratio.numerator // ratio.denominator
    remainder = ratio.numerator - lower * ratio.denominator
    twice = 2 * remainder
    if twice < ratio.denominator:
        return lower
    if twice > ratio.denominator:
        return lower + 1
    return lower if lower % 2 == 0 else lower + 1


def source_scores(directory):
    root = pathlib.Path(directory).resolve()
    manifest_bytes = (root / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest["frozen_protocol_sha256"] != PARENT_SHA:
        raise ValueError("cohort was acquired under another parent protocol")
    checked_files(root, manifest["files"])
    scores = {}
    formatting = Counter()
    for name in sorted(manifest["files"]):
        if not name.endswith(".scores.csv"):
            continue
        with (root / name).open(newline="") as handle:
            for row in csv.DictReader(handle):
                accession, raw = row["accession"], row["score"]
                try:
                    value = Fraction(raw)
                except (ValueError, ZeroDivisionError):
                    continue
                if accession in scores:
                    raise ValueError("duplicate source accession")
                scores[accession] = (value, raw)
                if "e" in raw.lower():
                    formatting["scientific_notation"] += 1
                else:
                    formatting["fractional_digits_" + str(
                        len(raw.partition(".")[2]))] += 1
    return scores, {"manifest_sha256": sha256(manifest_bytes),
                    "source_score_formatting": dict(formatting)}


def python_position_mean(values):
    if all(value == values[0] for value in values):
        return values[0]
    return float(np.mean(np.asarray(values, dtype=np.float64)))


def build_cases(rust, independent_predictions, independent_report, sources):
    python_rows = {}
    for row in independent_predictions:
        key = (row["scheme"], row["model"], row["assay_id"], row["accession"])
        if key in python_rows:
            raise ValueError("duplicate independent prediction identity")
        python_rows[key] = row
    seen = set()
    cases = []
    raw_archive_checks = []
    for family in rust["analyses"]:
        name = family["family"]
        for model in family["models"]:
            model_name = model["model"]
            grouped = defaultdict(list)
            for row in model["predictions"]:
                grouped[row["assay_id"]].append(row)
            positions = {(row["assay_id"], row["position_1based"]): row
                         for row in model["per_position_errors"]}
            metrics = {row["assay_id"]: row for row in model["per_assay_metrics"]}
            for assay, rows in sorted(grouped.items()):
                rows.sort(key=lambda row: (row["position_1based"],
                                          row["alternate"], row["id"]))
                observed, rust_pred, python_pred = [], [], []
                position_rows = defaultdict(list)
                for row in rows:
                    key = (name, model_name, assay, row["id"])
                    peer = python_rows.get(key)
                    if peer is None or key in seen:
                        raise ValueError("missing/duplicate matched prediction")
                    seen.add(key)
                    source, raw = sources[row["id"]]
                    if (row["observed"] != peer["observed"]
                            or row["observed"] != float(source)):
                        raise ValueError("saved observation differs from original source")
                    if "observed_source_decimal" in row:
                        if Fraction(row["observed_source_decimal"]) != source:
                            raise ValueError("Rust source-decimal commitment differs")
                    if any(not math.isfinite(v) for v in
                           [row["predicted"], peer["predicted"]]):
                        raise ValueError("nonfinite prediction")
                    observed.append(source)
                    rust_pred.append(row["predicted"])
                    python_pred.append(peer["predicted"])
                    position_rows[row["position_1based"]].append(
                        (source, row["predicted"], peer["predicted"], raw))
                native_python = independent_report["evaluation"][name]["assays"][assay][model_name]
                cases.append({
                    "family": name, "model": model_name, "assay_id": assay,
                    "resolution": "row", "observed": observed,
                    "rust": rust_pred, "independent": python_pred,
                    "native_rust_rho": metrics[assay]["row_metrics"]["spearman"],
                    "native_independent_rho": native_python["spearman_rows"],
                })
                observed_positions, rust_positions, python_positions = [], [], []
                for position, members in sorted(position_rows.items()):
                    exact_mean = sum((v[0] for v in members), Fraction()) / len(members)
                    stored = positions[(assay, position)]
                    if "observed_mean_exact_numerator" in stored:
                        stored_mean = Fraction(
                            int(stored["observed_mean_exact_numerator"]),
                            int(stored["observed_mean_exact_denominator"]))
                        if stored_mean != exact_mean:
                            raise ValueError("Rust exact position mean differs from source")
                    observed_positions.append(exact_mean)
                    rust_positions.append(stored["predicted_mean"])
                    python_positions.append(python_position_mean([v[2] for v in members]))
                cases.append({
                    "family": name, "model": model_name, "assay_id": assay,
                    "resolution": "position", "observed": observed_positions,
                    "rust": rust_positions, "independent": python_positions,
                    "native_rust_rho": metrics[assay]["position_mean_metrics"]["spearman"],
                    "native_independent_rho": native_python["spearman_positions"],
                })
    if seen != set(python_rows):
        raise ValueError("implementations have different prediction identity sets")
    if len(cases) != 240:
        raise ValueError("expected all240 prescribed row/position metric pairs")
    for case in cases:
        for implementation, native_key in [
            ("rust", "native_rust_rho"), ("independent", "native_independent_rho")
        ]:
            own = rank_correlation(case["observed"], case[implementation])
            native = case[native_key]
            if own is None or native is None:
                mismatch = own != native
                error = None
            else:
                error = abs(own - native)
                mismatch = error > TOLERANCE
            raw_archive_checks.append({
                **identity(case), "implementation": implementation,
                "recomputed": own, "native_archive": native,
                "absolute_difference": error, "mismatch": mismatch,
            })
    if any(row["mismatch"] for row in raw_archive_checks):
        raise ValueError("raw rank reader does not reproduce the native archived metrics")
    return cases, raw_archive_checks, len(seen)


def identity(case):
    return {key: case[key] for key in
            ("family", "model", "assay_id", "resolution")}


def sensitivity(cases):
    results = []
    summary = []
    for level in LEVELS:
        quantum = None if level == "raw_exact" else Fraction(level)
        level_rows = []
        for case in cases:
            left = case["rust"] if quantum is None else [
                quantized_integer(v, quantum) for v in case["rust"]]
            right = case["independent"] if quantum is None else [
                quantized_integer(v, quantum) for v in case["independent"]]
            rho_left = rank_correlation(case["observed"], left)
            rho_right = rank_correlation(case["observed"], right)
            status_match = (rho_left is None) == (rho_right is None)
            difference = None if rho_left is None or rho_right is None else abs(rho_left-rho_right)
            match = status_match and (difference is None or difference <= TOLERANCE)
            item = {
                **identity(case), "level": level, "targets": len(left),
                "rust_spearman": rho_left, "independent_spearman": rho_right,
                "rust_defined": rho_left is not None,
                "independent_defined": rho_right is not None,
                "undefined_status_match": status_match,
                "absolute_correlation_difference": difference,
                "within_reporting_tolerance": match,
                "rust_unique_predicted_codes": len(set(left)),
                "independent_unique_predicted_codes": len(set(right)),
                "observed_unique_exact_values": len(set(case["observed"])),
                "maximum_raw_prediction_difference": max(
                    abs(a-b) for a, b in zip(case["rust"], case["independent"])),
                "rust_raw_prediction_range": [min(case["rust"]), max(case["rust"])],
                "independent_raw_prediction_range": [min(case["independent"]),
                                                     max(case["independent"])],
            }
            results.append(item)
            level_rows.append(item)
        finite = [r for r in level_rows if r["absolute_correlation_difference"] is not None]
        worst = sorted(level_rows, key=lambda row: (
            not row["undefined_status_match"],
            row["absolute_correlation_difference"] or 0.0), reverse=True)[:10]
        summary.append({
            "level": level, "metric_pairs": len(level_rows),
            "matched": sum(r["within_reporting_tolerance"] for r in level_rows),
            "mismatched": sum(not r["within_reporting_tolerance"] for r in level_rows),
            "undefined_status_mismatches": sum(not r["undefined_status_match"]
                                               for r in level_rows),
            "both_undefined": sum(not r["rust_defined"] and not r["independent_defined"]
                                  for r in level_rows),
            "maximum_defined_correlation_difference": max(
                (r["absolute_correlation_difference"] for r in finite), default=None),
            "worst_cases": worst,
        })
    return results, summary


def write_new(path, value):
    with pathlib.Path(path).open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def run(args):
    out = pathlib.Path(args.out)
    if out.exists():
        raise ValueError("output exists; preserve previous numerical attempts")
    protocol_path = pathlib.Path(args.protocol)
    protocol_raw = protocol_path.read_bytes()
    protocol = json.loads(protocol_raw)
    if (sha256(protocol_raw) != PROTOCOL_SHA or protocol["status"] != "frozen"
            or protocol["levels"] != LEVELS):
        raise ValueError("sensitivity protocol differs from fixed pre-run specification")
    parent_raw = pathlib.Path(args.parent_protocol).read_bytes()
    if sha256(parent_raw) != PARENT_SHA:
        raise ValueError("parent frozen protocol changed")
    independent = pathlib.Path(args.independent).resolve()
    artifacts_raw = (independent / "artifacts.json").read_bytes()
    checked_files(independent, json.loads(artifacts_raw)["files"])
    independent_report = read_json(independent / "replication.json")
    rust = read_json(args.rust)
    if (rust["protocol_sha256"] != PARENT_SHA
            or independent_report["protocol_sha256"] != PARENT_SHA):
        raise ValueError("saved implementations used different parent protocols")
    sources, source_provenance = source_scores(args.cohort)
    cases, raw_checks, matched = build_cases(
        rust, read_json(independent / "predictions.json"),
        independent_report, sources)
    print(f"Validated {matched} matched predictions and {len(cases)} metric pairs; computing all fixed levels.", flush=True)
    details, summaries = sensitivity(cases)
    provenance = {
        "sensitivity_protocol_sha256": sha256(protocol_raw),
        "parent_protocol_sha256": PARENT_SHA,
        "rust_study_path": str(pathlib.Path(args.rust).resolve()),
        "rust_study_sha256": file_hash(args.rust),
        "independent_directory": str(independent),
        "independent_artifacts_sha256": sha256(artifacts_raw),
        "source_cohort": source_provenance,
        "implementation_sha256": file_hash(__file__),
        "python": sys.version, "platform": platform.platform(),
        "numpy_version": np.__version__,
    }
    out.mkdir(parents=True)
    write_new(out / "summary.json", {
        "schema_version": 1,
        "purpose": "Post-hoc numerical sensitivity; no v1 endpoint replacement",
        "provenance": provenance,
        "matched_prediction_identities": matched,
        "raw_native_metric_checks": len(raw_checks),
        "all_levels": summaries,
        "v1_raw_exact_rank_passes": summaries[0]["mismatched"] == 0,
        "selected_quantum": None,
        "interpretation": "Each coarse grid is a separately labeled supplemental rank estimand. No quantum is selected for v1 success and no biological inference/p-value is made.",
    })
    write_new(out / "all-metric-pairs.json", details)
    write_new(out / "raw-archive-checks.json", raw_checks)
    files = [p.name for p in out.iterdir()]
    write_new(out / "artifacts.json", {"schema_version": 1,
                                      "files": {name: file_hash(out/name) for name in sorted(files)}})
    print(json.dumps([{k: v for k, v in row.items() if k != "worst_cases"}
                      for row in summaries], indent=2))


def self_test():
    q = Fraction("0.25")
    assert quantized_integer(0.125, q) == 0
    assert quantized_integer(0.375, q) == 2
    assert quantized_integer(-0.125, q) == 0
    assert quantized_integer(-0.375, q) == -2
    assert average_ranks([Fraction("0.1"), Fraction("0.10"), Fraction("0.2")]) == [1.5, 1.5, 3.0]
    assert rank_correlation([1, 2, 3], [0, 0, 0]) is None
    assert rank_correlation([1, 2, 3], [1, 2, 3]) == 1.0
    assert Fraction("0.01")+Fraction("0.02")+Fraction("0.03") == 3*Fraction("0.02")
    print("Exact ties-even quantization, source-decimal ties and undefined-rank tests passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", default="experiments/numerical-rank-sensitivity-v1.json")
    parser.add_argument("--parent-protocol", default="experiments/regulatory-prediction-v1.json")
    parser.add_argument("--cohort", default="datasets/regulatory-cohort")
    parser.add_argument("--rust")
    parser.add_argument("--independent")
    parser.add_argument("--out")
    parser.add_argument("--self-test", action="store_true")
    options = parser.parse_args()
    if options.self_test:
        self_test()
    else:
        if not options.rust or not options.independent or not options.out:
            parser.error("--rust, --independent and --out are required for the analysis")
        run(options)
