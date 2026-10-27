#!/usr/bin/env python3
"""Compare independent context diagnostics and primary numerical status."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path


def run(study, independent, comparison):
    study, independent = Path(study), Path(independent)
    digest = lambda data: hashlib.sha256(data).hexdigest()
    for root in (study, independent):
        manifest = json.loads((root / "artifacts.json").read_bytes())
        for name, expected in manifest["files"].items():
            path = (root / name).resolve()
            if root.resolve() not in path.parents or digest(path.read_bytes()) != expected:
                raise ValueError("artifact path/checksum mismatch: " + name)
    context = json.loads((study / "context.json").read_text())
    replica = json.loads((independent / "replication.json").read_text())
    compared = json.loads((Path(comparison) / "comparison.json").read_text())
    raw_study = (study / "study.json.gz").read_bytes()
    if (compared["rust_study_sha256"] != digest(raw_study)
            or compared["rust_decoded_study_sha256"] != digest(gzip.decompress(raw_study))
            or compared["independent_artifacts_sha256"] != digest((independent / "artifacts.json").read_bytes())):
        raise ValueError("comparison is not bound to the supplied parent results")
    if compared["primary_predictions_losses_and_controls_matched"] is not True:
        raise ValueError("primary numerical comparison failed")
    protocol_bytes = (study / "protocol.json").read_bytes()
    if compared["protocol_sha256"] != digest(protocol_bytes):
        raise ValueError("comparison protocol differs")
    protocol = json.loads(protocol_bytes)
    expected_pairs = {(x["left"], x["right"])
                      for x in protocol["validation_context_analysis"]["pairs"]}
    actual_pairs = {(x["left"], x["right"]) for x in context["pairs"]}
    by_pair = {(x["left"], x["right"]): x for x in replica["context_pairs"]}
    if (actual_pairs != expected_pairs or set(by_pair) != expected_pairs
            or len(context["pairs"]) != len(expected_pairs)
            or len(replica["context_pairs"]) != len(expected_pairs)):
        raise ValueError("context pair inventory differs from frozen protocol")
    for item in context["pairs"]:
        other = by_pair[(item["left"], item["right"])]
        if item["matched_variants"] != other["matched_variants"]:
            raise ValueError("matched context counts differ")
        for field in ["mse", "mae", "signed_mean_left_minus_right", "pearson", "spearman",
                      "finite_context_blind_mse_lower_bound"]:
            if item[field] is None or other[field] is None:
                if item[field] is not other[field]:
                    raise ValueError("undefined context metric differs: " + field)
            elif (not math.isfinite(item[field]) or not math.isfinite(other[field])
                  or abs(item[field] - other[field]) > 1e-10):
                raise ValueError("context metric differs: " + field)
    print("All five context-pair counts, metrics and finite source-label bounds reproduce.")
    if not compared["secondary_exact_rank_metrics_matched"]:
        print("Frozen exact-rank secondary endpoint remains unreproduced; the full comparison is failed.")
    return compared["matched"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True)
    parser.add_argument("--independent", required=True)
    parser.add_argument("--comparison", required=True)
    args = parser.parse_args()
    run(args.study, args.independent, args.comparison)
