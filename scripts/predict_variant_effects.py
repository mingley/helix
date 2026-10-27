#!/usr/bin/env python3
"""Outcome-free, zero-refit transfer of five frozen activity folds to 200nt SNVs.

The output is a reporter-scale model prediction, never a clinical interpretation.
All artifact files are checked against a caller-pinned manifest before fit parsing.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path, PurePosixPath
import re
import sys

PROTOCOL_SHA256 = "64e14037fb1f42ddc1ecaca1d23bb8096fabac1f1e892d105eaa20f525fbbb20"
SOURCE_SHA256 = "9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e"
SOURCE_MANIFEST_SHA256 = "6a3371bed4061f9ca9fa1a95bfc0947ffee91daa22fe509063708e549d453730"
CONTEXTS = ["HepG2", "K562", "WTC11"]
HEX = re.compile(r"[0-9a-f]{64}\Z")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
DIRECT_CHECK_ABSOLUTE_TOLERANCE = 1e-10


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def strict_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load_json(path, compressed=False):
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return json.load(handle, object_pairs_hook=strict_object,
                         parse_constant=lambda value: (_ for _ in ()).throw(
                             ValueError(f"nonfinite JSON literal: {value}")))


def canonical_json(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


TRIMERS = sorted({min("".join(k), reverse_complement("".join(k)))
                  for k in itertools.product("ACGT", repeat=3)})
FEATURE_NAMES = ["rc3mer_" + k for k in TRIMERS] + ["gc_fraction"]


def validate_dna(sequence):
    require(isinstance(sequence, str) and len(sequence) == 200
            and set(sequence) <= set("ACGT"), "DNA must be 200 uppercase ACGT bases")


def feature_counts(sequence):
    """Exact overlapping RC-canonical 3mer counts, followed by GC count."""
    validate_dna(sequence)
    counts = dict.fromkeys(TRIMERS, 0)
    for offset in range(198):
        word = sequence[offset:offset + 3]
        counts[min(word, reverse_complement(word))] += 1
    return [counts[word] for word in TRIMERS] + [sum(base in "GC" for base in sequence)]


def features(sequence):
    """All overlapping RC-canonical 3mer fractions, then GC; no fitting."""
    counts = feature_counts(sequence)
    return [count / 198 for count in counts[:-1]] + [counts[-1] / 200]


def numeric(value):
    return type(value) in (int, float) and math.isfinite(value)


def validate_ridge(fit):
    require(isinstance(fit, dict) and fit.get("feature_count") == 33
            and type(fit.get("training_rows")) is int and fit["training_rows"] > 0
            and fit.get("alpha") == 10.0 and numeric(fit.get("intercept")),
            "invalid saved ridge metadata")
    for key in ["training_feature_means", "training_feature_scales", "standardized_coefficients"]:
        require(isinstance(fit.get(key), list) and len(fit[key]) == 33
                and all(numeric(value) for value in fit[key]), f"invalid {key}")
    require(all(value > 0 for value in fit["training_feature_scales"]),
            "training scales must be positive")


def safe_member(root, name):
    require(isinstance(name, str), "non-string artifact member")
    part = PurePosixPath(name)
    require(not part.is_absolute() and part.parts and ".." not in part.parts
            and "\\" not in name and str(part) == name, "unsafe artifact member path")
    result = root / name
    require(not root.is_symlink() and all(not p.is_symlink()
            for p in [result, *result.parents] if p != root.parent), "symlink in artifact path")
    require(result.is_file() and result.resolve().is_relative_to(root.resolve()),
            "missing or escaped artifact member")
    return result


def load_frozen_fits(root, expected_manifest_sha256):
    """Bind manifest -> all files -> source/protocol -> five exact fitted objects."""
    root = Path(root)
    require(bool(HEX.fullmatch(expected_manifest_sha256)), "invalid expected manifest SHA256")
    manifest_path = safe_member(root, "artifacts.json")
    require(sha256(manifest_path) == expected_manifest_sha256, "artifact manifest SHA256 mismatch")
    manifest = load_json(manifest_path)
    require(manifest.get("schema_version") == 5 and isinstance(manifest.get("files"), dict),
            "expected activity-v3 artifact schema5")
    files = manifest["files"]
    needed = {"study.json.gz", "protocol.json", "source-manifest.json", "inputs.json"}
    require(needed <= files.keys(), "missing required artifact commitments")
    for name, digest in files.items():
        require(isinstance(digest, str) and bool(HEX.fullmatch(digest)), "invalid artifact SHA256")
        require(sha256(safe_member(root, name)) == digest, f"artifact SHA256 mismatch: {name}")
    require(files["protocol.json"] == PROTOCOL_SHA256, "unexpected frozen training protocol")
    require(files["source-manifest.json"] == SOURCE_MANIFEST_SHA256,
            "unexpected frozen training source manifest")
    protocol = load_json(root / "protocol.json")
    require(protocol.get("protocol_id") == "context-activity-v3", "wrong training protocol ID")
    source = load_json(root / "source-manifest.json")
    inputs = load_json(root / "inputs.json")
    require(source.get("endpoint") == inputs.get("endpoint")
            == "published_joint_log2_rna_dna_activity", "unexpected training endpoint")
    require(source.get("protocol_sha256") == inputs.get("protocol_sha256") == PROTOCOL_SHA256,
            "source/protocol binding mismatch")
    require(source.get("source", {}).get("sha256") == inputs.get("raw_source_sha256")
            == SOURCE_SHA256 and source["source"].get("bytes") == inputs.get("raw_source_bytes")
            == 53557805, "raw training source identity mismatch")
    require(inputs.get("source_manifest_sha256") == files["source-manifest.json"],
            "source manifest cross-binding mismatch")
    require(isinstance(source.get("gates"), dict) and set(source["gates"]) == {
        "archive_integrity_verified", "empirical_activity_schema_bound", "coverage_complete",
        "design_identity_verified", "replicate_aggregation_status_recorded"}
            and all(value is True for value in source["gates"].values()), "training source gate failed")
    require(isinstance(source.get("files"), dict) and all(files.get(name) == digest
            for name, digest in source["files"].items()), "dataset/artifact file binding mismatch")
    study = load_json(root / "study.json.gz", compressed=True)
    require(study.get("schema_version") == 1 and study.get("contexts") == CONTEXTS
            and study.get("alpha") == 10.0 and len(study.get("folds", [])) == 5,
            "unexpected saved study contract")
    require(study.get("panel_metadata_sha256") == inputs.get("panel_metadata_sha256")
            and study.get("outcomes_sha256") == inputs.get("outcomes_sha256"),
            "study/input binding mismatch")
    fits = []
    for index, fold in enumerate(study["folds"]):
        require(type(fold.get("fold_index_0based")) is int
                and fold["fold_index_0based"] == index, "fold order differs")
        fit = fold["fit"]
        require(fit.get("feature_names") == FEATURE_NAMES and fit.get("contexts") == CONTEXTS
                and fit.get("alpha") == 10.0
                and fit.get("panel_metadata_sha256") == study["panel_metadata_sha256"],
                "fit feature/context/source mismatch")
        require(set(fit.get("training_context_means", {})) == set(CONTEXTS)
                and all(numeric(value) for value in fit["training_context_means"].values()),
                "invalid context intercepts")
        validate_ridge(fit["shared_dna_fit"])
        require(set(fit.get("independent_context_fits", {})) == set(CONTEXTS),
                "missing context-specific fit")
        for ridge in fit["independent_context_fits"].values():
            validate_ridge(ridge)
            require(ridge["training_feature_means"] == fit["shared_dna_fit"]["training_feature_means"]
                    and ridge["training_feature_scales"] == fit["shared_dna_fit"]["training_feature_scales"],
                    "context scaler differs from shared scaler")
        fits.append(fit)
    return fits, {"artifacts_sha256": expected_manifest_sha256,
                  "study_sha256": files["study.json.gz"], "protocol_sha256": PROTOCOL_SHA256,
                  "source_manifest_sha256": files["source-manifest.json"],
                  "raw_training_source_sha256": SOURCE_SHA256}


def validate_pairs(payload):
    require(isinstance(payload, dict) and set(payload) == {"schema_version", "pairs"}
            and type(payload["schema_version"]) is int and payload["schema_version"] == 1,
            "input must contain only schema_version=1 and pairs; outcomes are not accepted")
    pairs = payload["pairs"]
    require(isinstance(pairs, list) and 0 < len(pairs) <= 1000000, "need 1..1000000 pairs")
    seen = set()
    for pair in pairs:
        require(isinstance(pair, dict) and set(pair) == {
            "id", "context", "position_0based", "reference", "alternate"},
            "pair must contain only id/context/position_0based/reference/alternate")
        for label in ["id", "context"]:
            require(isinstance(pair[label], str) and pair[label] and pair[label].strip() == pair[label],
                    f"invalid {label}")
        require(pair["id"] not in seen, "duplicate pair ID")
        seen.add(pair["id"])
        position = pair["position_0based"]
        require(type(position) is int and 0 <= position < 200, "invalid 0-based SNV position")
        validate_dna(pair["reference"])
        validate_dna(pair["alternate"])
        different = [j for j, (a, b) in enumerate(zip(pair["reference"], pair["alternate"])) if a != b]
        require(different == [position], "pair must differ by exactly one declared SNV")
    return pairs


def predict_direct(ridge, row, intercept=None):
    return (ridge["intercept"] if intercept is None else intercept) + math.fsum(
        (value - mean) / scale * beta for value, mean, scale, beta in zip(
            row, ridge["training_feature_means"], ridge["training_feature_scales"],
            ridge["standardized_coefficients"]))


def score_pair(pair, fits, model):
    require(model in {"S", "A", "I"}, "model must be S, A or I")
    require(len(fits) == 5, "exactly five frozen folds required")
    if model == "I":
        require(pair["context"] in CONTEXTS, "I requires exact training cell context")
    if model == "A":
        require(pair["context"] in CONTEXTS, "A requires exact training cell context")
    ref_counts, alt_counts = feature_counts(pair["reference"]), feature_counts(pair["alternate"])
    denominators = [198] * 32 + [200]
    ref = [count / denominator for count, denominator in zip(ref_counts, denominators)]
    alt = [count / denominator for count, denominator in zip(alt_counts, denominators)]
    # Subtract exact counts before division: the primary delta then has no
    # rounding dependence on the otherwise unchanged remote sequence context.
    delta = [(b - a) / denominator for a, b, denominator in zip(ref_counts, alt_counts, denominators)]
    folds, gc_terms, checks = [], [], []
    for fit in fits:
        ridge = (fit["independent_context_fits"][pair["context"]] if model == "I"
                 else fit["shared_dna_fit"])
        slopes = [beta / scale for beta, scale in zip(
            ridge["standardized_coefficients"], ridge["training_feature_scales"])]
        prediction = math.fsum(slope * difference for slope, difference in zip(slopes, delta))
        intercept = fit["training_context_means"][pair["context"]] if model == "A" else None
        direct = predict_direct(ridge, alt, intercept) - predict_direct(ridge, ref, intercept)
        require(abs(prediction - direct) <= DIRECT_CHECK_ABSOLUTE_TOLERANCE,
                "direct prediction subtraction cross-check failed")
        folds.append(prediction)
        gc_terms.append(slopes[-1] * delta[-1])
        checks.append(abs(prediction - direct))
    return {"id": pair["id"], "context": pair["context"],
            "position_0based": pair["position_0based"],
            "reference_sha256": hashlib.sha256(pair["reference"].encode()).hexdigest(),
            "alternate_sha256": hashlib.sha256(pair["alternate"].encode()).hexdigest(),
            "reference_feature_numerators": ref_counts,
            "alternate_feature_numerators": alt_counts,
            "feature_numerator_delta": [b - a for a, b in zip(ref_counts, alt_counts)],
            "predicted_alternate_minus_reference": math.fsum(folds) / 5,
            "fold_predictions": folds,
            "learned_gc_term_only_ablation_not_refitted_gc_model": math.fsum(gc_terms) / 5,
            "zero_delta_baseline": 0.0,
            "max_abs_direct_subtraction_check_error": max(checks)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-artifacts", type=Path, required=True)
    parser.add_argument("--expected-artifacts-sha256", required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--model", choices=["S", "A", "I"], required=True)
    parser.add_argument("--out", type=Path, required=True, help="new output directory")
    args = parser.parse_args(argv)
    require(not args.out.exists(), "output path already exists")
    fits, commitments = load_frozen_fits(args.model_artifacts, args.expected_artifacts_sha256)
    pairs = validate_pairs(load_json(args.pairs))
    rows = [score_pair(pair, fits, args.model) for pair in pairs]
    output = {"schema_version": 1, "model": args.model,
              "ensemble": "arithmetic mean of all five saved grouped-training fits; no refit",
              "endpoint_direction": "alternate minus reference",
              "units": "transferred published library-normalized log2(RNA/DNA) model scale",
              "commitments": commitments | {"pairs_sha256": sha256(args.pairs)},
              "feature_names": FEATURE_NAMES, "pairs": rows,
              "feature_numerator_definition": "overlapping RC3mer counts, then GC base count",
              "feature_denominators": [198] * 32 + [200],
              "direct_check_absolute_tolerance": DIRECT_CHECK_ABSOLUTE_TOLERANCE,
              "limits": ["No external outcome is accepted or evaluated.",
                         "Prediction is not a measured or clinical effect.",
                         "200nt same-length SNVs only; no indels, recentering or adapters.",
                         "RC3mer/GC models omit orientation and interactions beyond local 5mer for interior SNVs.",
                         "GC term is a coefficient ablation, not an independently fitted GC-only baseline.",
                         "Fold spread is not a calibrated uncertainty interval.",
                         "External assay scale/context transfer and locus independence require a separate frozen study."]}
    args.out.mkdir(parents=True)
    (args.out / "predictions.json").write_bytes(canonical_json(output))
    implementation = Path(__file__).resolve()
    receipt = {"schema_version": 1, "implementation_sha256": sha256(implementation),
               "pairs_sha256": sha256(args.pairs), "training": commitments,
               "files": {"predictions.json": sha256(args.out / "predictions.json")}}
    (args.out / "manifest.json").write_bytes(canonical_json(receipt))
    print(json.dumps({"pairs": len(rows), "manifest_sha256": sha256(args.out / "manifest.json")}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, TypeError, OSError, json.JSONDecodeError) as error:
        print(f"variant prediction failed: {error}", file=sys.stderr)
        raise SystemExit(2)
