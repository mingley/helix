#!/usr/bin/env python3
"""Outcome-free exhaustive envelope of fixed 200nt RC3mer+GC SNV predictions.

Independent numeric-coded feature changes and raw-slope arithmetic; the separate
transfer predictor is imported solely to compare the completed predictions.
No fitting or external measurements are accepted. Run from the repository root.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path, PurePosixPath
import shutil

ARTIFACT_SHA = "2ea476bc5413dc13486dec9f3d7e9162a50a0e47f1d83ce29efabc5f5cc5af9d"
TRAINING_PROTOCOL_SHA = "64e14037fb1f42ddc1ecaca1d23bb8096fabac1f1e892d105eaa20f525fbbb20"
SOURCE_MANIFEST_SHA = "6a3371bed4061f9ca9fa1a95bfc0947ffee91daa22fe509063708e549d453730"
BASES = "ACGT"
CELLS = ("HepG2", "K562", "WTC11")
MODELS = ("S", "I:HepG2", "I:K562", "I:WTC11")
TOLERANCE = 1e-10


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def unique_object(items):
    out = {}
    for key, value in items:
        require(key not in out, "duplicate JSON key")
        out[key] = value
    return out


def read_json(path, compressed=False):
    with (gzip.open if compressed else open)(path, "rt", encoding="utf-8") as handle:
        return json.load(handle, object_pairs_hook=unique_object,
                         parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def word_code(word):
    value = 0
    for base in word:
        value = 4 * value + BASES.index(base)
    return value


def code_word(value):
    return BASES[value // 16] + BASES[value // 4 % 4] + BASES[value % 4]


def reverse_code(value):
    return (3 - value % 4) * 16 + (3 - value // 4 % 4) * 4 + (3 - value // 16)


CLASSES = tuple(value for value in range(64) if value <= reverse_code(value))
CLASS_INDEX = {value: index for index, value in enumerate(CLASSES)}
FEATURE_NAMES = ["rc3mer_" + code_word(value) for value in CLASSES] + ["gc_fraction"]


def trimer_index(word):
    value = word_code(word)
    return CLASS_INDEX[min(value, reverse_code(value))]


def exact_counts(sequence):
    """Independent integer counts; no predictor feature function is used."""
    require(set(sequence) <= set(BASES), "non-DNA character")
    counts = [0] * 33
    for start in range(len(sequence) - 2):
        counts[trimer_index(sequence[start:start + 3])] += 1
    counts[-1] = sum(base in ("C", "G") for base in sequence)
    return counts


def local_delta(original, changed):
    """Only three physical windows can contain the central substituted base."""
    counts = [0] * 33
    for start in range(3):
        counts[trimer_index(original[start:start + 3])] -= 1
        counts[trimer_index(changed[start:start + 3])] += 1
    counts[-1] = int(changed[2] in ("C", "G")) - int(original[2] in ("C", "G"))
    return counts


def validate_source(root, expected):
    require(expected == ARTIFACT_SHA, "this version requires the pinned clean activity artifact")
    require(root.is_dir() and not root.is_symlink(), "invalid artifact root")
    manifest = root / "artifacts.json"
    require(digest(manifest) == expected, "source artifact-manifest SHA mismatch")
    data = read_json(manifest)
    require(data["schema_version"] == 5 and isinstance(data["files"], dict), "wrong artifact schema")
    for name, sha in data["files"].items():
        part = PurePosixPath(name)
        require(not part.is_absolute() and ".." not in part.parts and str(part) == name
                and "\\" not in name, "unsafe artifact path")
        file = root / name
        require(file.is_file() and not file.is_symlink()
                and all(not p.is_symlink() for p in file.parents if p != root.parent)
                and file.resolve().is_relative_to(root.resolve()), "missing or escaped artifact")
        require(digest(file) == sha, "source artifact byte mismatch: " + name)
    require(data["files"]["protocol.json"] == TRAINING_PROTOCOL_SHA, "training protocol mismatch")
    require(data["files"]["source-manifest.json"] == SOURCE_MANIFEST_SHA, "source dataset mismatch")
    return data


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def extract_fits(root):
    # This container also retains previously exposed training/heldout values.
    # Only the five saved fit objects are selected; no outcomes enter arithmetic.
    study = read_json(root / "study.json.gz", compressed=True)
    require(study["contexts"] == list(CELLS) and study["alpha"] == 10.0
            and len(study["folds"]) == 5, "unexpected training-study structure")
    fits = []
    for number, fold in enumerate(study["folds"]):
        require(fold["fold_index_0based"] == number, "incorrect fold order")
        fit = fold["fit"]
        require(fit["feature_names"] == FEATURE_NAMES and fit["contexts"] == list(CELLS),
                "feature or cell order mismatch")
        chosen = {"S": fit["shared_dna_fit"]} | {
            "I:" + cell: fit["independent_context_fits"][cell] for cell in CELLS}
        for ridge in chosen.values():
            require(ridge["feature_count"] == 33 and ridge["alpha"] == 10.0
                    and finite(ridge["intercept"]), "invalid ridge metadata")
            for name in ("training_feature_means", "training_feature_scales", "standardized_coefficients"):
                require(len(ridge[name]) == 33 and all(finite(x) for x in ridge[name]), "invalid saved vector")
            require(all(x > 0 for x in ridge["training_feature_scales"]), "nonpositive scale")
        fits.append(chosen)
    return fits


def raw_slopes(fits):
    return [{model: [coefficient / scale for coefficient, scale in zip(
        fold[model]["standardized_coefficients"], fold[model]["training_feature_scales"])]
        for model in MODELS} for fold in fits]


def predict_delta(delta, slopes, model):
    # Counts are subtracted exactly before division. Remote bases never enter.
    fractions = [value / (200 if index == 32 else 198) for index, value in enumerate(delta)]
    folds = [math.fsum(a * b for a, b in zip(fold[model], fractions)) for fold in slopes]
    gc = math.fsum(fold[model][-1] * fractions[-1] for fold in slopes) / 5
    return folds, math.fsum(folds) / 5, gc


def write_csv(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-artifacts", type=Path,
                        default=Path("research/results/context-activity-v3/core"))
    parser.add_argument("--expected-artifacts-sha256", default=ARTIFACT_SHA)
    parser.add_argument("--comparator", type=Path, default=Path("scripts/predict_variant_effects.py"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    require(not args.out.exists(), "output directory already exists; retain every attempt")
    args.out.mkdir(parents=True)
    shutil.copyfile(Path(__file__), args.out / "implementation.py")
    shutil.copyfile(args.comparator, args.out / "comparison-implementation.py")
    protocol = {
        "protocol_id": "variant-locality-model-envelope-v1", "schema_version": 1,
        "kind": "exhaustive synthetic architecture audit; no external outcomes or refitting",
        "artifact_manifest_sha256": ARTIFACT_SHA, "training_protocol_sha256": TRAINING_PROTOCOL_SHA,
        "training_source_manifest_sha256": SOURCE_MANIFEST_SHA,
        "implementation_sha256": digest(args.out / "implementation.py"),
        "comparison_implementation_sha256": digest(args.out / "comparison-implementation.py"),
        "enumeration": "lexicographic ACGT^5; replace central base with each other ACGT base",
        "edits": 3072, "models": list(MODELS), "folds_per_model": 5,
        "embedding": "97 A bases + local 5mer + 98 A bases; edit zero-based position99 in fixed200nt",
        "prediction": "exact local integer count changes divided by198 (trimers) or200 (GC), dot beta/scale; mean all5folds",
        "comparison": "separate predictor is used only after independent results; all5fold and ensemble deltas plus learned GC terms",
        "absolute_comparison_tolerance": TOLERANCE, "expected_exact_feature_zero_edits": 64,
        "extrema": "retain all ties for ensemble minimum, maximum and maximumabsolute permodel; allrowdata retained",
        "units": "transferred published library-normalized log2(RNA/DNA) model scale",
        "limits": ["Existing training outcomes were previously exposed; no new outcome payload accepted.",
                   "Interior same-length SNVs only; this is not a whole-sequence activity bound.",
                   "Prediction is not an observed, causal or clinical variant effect.",
                   "Five fits share one study; fold spread is not calibrated uncertainty.",
                   "Floating arithmetic is cross-checked, not formally verified.",
                   "This finite model envelope does not establish global scientific novelty."]}
    (args.out / "protocol.json").write_bytes(encoded(protocol))
    protocol_sha = digest(args.out / "protocol.json")
    try:
        source = validate_source(args.model_artifacts, args.expected_artifacts_sha256)
        fits = extract_fits(args.model_artifacts)
        slopes = raw_slopes(fits)
        fit_payload = {"schema_version": 1, "models": list(MODELS), "feature_names": FEATURE_NAMES,
                       "source_artifacts_sha256": ARTIFACT_SHA, "source_study_sha256": source["files"]["study.json.gz"],
                       "saved_ridges": fits, "raw_slopes": slopes}
        (args.out / "saved-fits.json").write_bytes(encoded(fit_payload))
        specification = importlib.util.spec_from_file_location("envelope_comparator", args.out / "comparison-implementation.py")
        comparator = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(comparator)
        # Its loader is never used for the primary result. Invoke its source
        # checks separately, then keep it solely as a comparison implementation.
        comparator_fits, comparator_bindings = comparator.load_frozen_fits(args.model_artifacts, ARTIFACT_SHA)
        predictions, comparisons, zeros = [], [], []
        exact_full_count_checks = 0
        maximum_difference = 0.0
        numeric_prediction_comparisons = 0
        all_results = {model: [] for model in MODELS}
        for bases in itertools.product(BASES, repeat=5):
            original = "".join(bases)
            for alternate in BASES:
                if alternate == original[2]:
                    continue
                changed = original[:2] + alternate + original[3:]
                delta = local_delta(original, changed)
                reference = "A" * 97 + original + "A" * 98
                mutant = "A" * 97 + changed + "A" * 98
                full_delta = [b - a for a, b in zip(exact_counts(reference), exact_counts(mutant))]
                require(delta == full_delta, "local/full independent count difference")
                exact_full_count_checks += 33
                id_value = original + ">" + changed
                zero = not any(delta)
                zero_maximum = 0.0
                for model in MODELS:
                    fold_predictions, mean_prediction, gc_prediction = predict_delta(delta, slopes, model)
                    cell = model.split(":")[1] if model != "S" else "shared"
                    row = {"edit_id": id_value, "original_5mer": original, "alternate_5mer": changed,
                           "reference_base": original[2], "alternate_base": alternate, "model": model,
                           "context": cell, "ensemble_prediction": mean_prediction,
                           "learned_gc_term": gc_prediction} | {
                               "fold_" + str(index): value for index, value in enumerate(fold_predictions)}
                    predictions.append(row)
                    all_results[model].append(row)
                    pair = {"id": id_value, "context": cell if model != "S" else "HepG2",
                            "position_0based": 99, "reference": reference, "alternate": mutant}
                    checked = comparator.score_pair(pair, comparator_fits, "S" if model == "S" else "I")
                    expected = checked["fold_predictions"] + [checked["predicted_alternate_minus_reference"],
                        checked["learned_gc_term_only_ablation_not_refitted_gc_model"]]
                    actual = fold_predictions + [mean_prediction, gc_prediction]
                    errors = [abs(a - b) for a, b in zip(actual, expected)]
                    require(all(error <= TOLERANCE for error in errors), "numeric predictor comparison failed")
                    maximum_difference = max(maximum_difference, *errors)
                    numeric_prediction_comparisons += len(errors)
                    comparisons.append({"edit_id": id_value, "model": model,
                        "maximum_absolute_prediction_difference": max(errors),
                        "comparator_direct_subtraction_error": checked["max_abs_direct_subtraction_check_error"],
                        "within_tolerance": True})
                    if zero:
                        require(all(value == 0.0 for value in actual + expected), "feature-zero edit predicted nonzero")
                        zero_maximum = max(zero_maximum, *(abs(x) for x in actual + expected))
                if zero:
                    reversed_word = "".join(BASES[3 - BASES.index(base)] for base in original[::-1])
                    require(changed == reversed_word, "unexpected feature-zero equivalence")
                    zeros.append({"edit_id": id_value, "original_5mer": original, "alternate_5mer": changed,
                                  "all_33_feature_deltas_exactly_zero": True,
                                  "alternate_is_reverse_complement": True,
                                  "all_20_fold_and_4_ensemble_predictions_exactly_zero": True,
                                  "maximum_absolute_prediction": zero_maximum})
        require(len(predictions) == 3072 * 4 and len(zeros) == 64, "incomplete enumeration")
        write_csv(args.out / "all-snv-predictions.csv",
            ["edit_id", "original_5mer", "alternate_5mer", "reference_base", "alternate_base", "model", "context"]
            + ["fold_" + str(i) for i in range(5)] + ["ensemble_prediction", "learned_gc_term"], predictions)
        write_csv(args.out / "numeric-comparisons.csv", list(comparisons[0]), comparisons)
        write_csv(args.out / "zero-feature-checks.csv", list(zeros[0]), zeros)
        extrema = {}
        for model, rows in all_results.items():
            minimum = min(row["ensemble_prediction"] for row in rows)
            maximum = max(row["ensemble_prediction"] for row in rows)
            absolute = max(abs(row["ensemble_prediction"]) for row in rows)
            extrema[model] = {"minimum": minimum, "maximum": maximum, "maximum_absolute": absolute,
                "minimum_examples": [row for row in rows if row["ensemble_prediction"] == minimum],
                "maximum_examples": [row for row in rows if row["ensemble_prediction"] == maximum],
                "maximum_absolute_examples": [row for row in rows if abs(row["ensemble_prediction"]) == absolute],
                "per_fold_maximum_absolute": [max(abs(row["fold_" + str(i)]) for row in rows) for i in range(5)]}
        summary = {"schema_version": 1, "protocol_sha256": protocol_sha,
                   "source_artifacts_sha256": ARTIFACT_SHA, "source_study_sha256": source["files"]["study.json.gz"],
                   "saved_fits_sha256": digest(args.out / "saved-fits.json"),
                   "implementation_sha256": digest(args.out / "implementation.py"),
                   "comparison_implementation_sha256": digest(args.out / "comparison-implementation.py"),
                   "oriented_5mers": 1024, "central_SNVs": 3072, "models": list(MODELS),
                   "model_edit_rows": len(predictions), "fold_predictions": len(predictions) * 5,
                   "ensemble_predictions": len(predictions), "exact_full_count_checks": exact_full_count_checks,
                   "exact_feature_zero_edits": len(zeros), "feature_zero_checks_passed": True,
                   "numeric_prediction_comparisons": numeric_prediction_comparisons,
                   "maximum_absolute_predictor_difference": maximum_difference,
                   "comparison_absolute_tolerance": TOLERANCE, "all_numeric_comparisons_passed": True,
                   "new_experimental_outcomes_acquired": False, "models_refitted": False,
                   "extrema": extrema, "limits": protocol["limits"]}
        (args.out / "summary.json").write_bytes(encoded(summary))
        print(json.dumps({"rows": len(predictions), "zero_edits": len(zeros),
                          "max_comparison_difference": maximum_difference,
                          "max_abs_by_model": {key: value["maximum_absolute"] for key, value in extrema.items()}}))
    except BaseException as error:
        (args.out / "failure.json").write_bytes(encoded({"error_type": type(error).__name__, "error": str(error)}))
        raise
    finally:
        files = {p.name: {"bytes": p.stat().st_size, "sha256": digest(p)}
                 for p in sorted(args.out.iterdir()) if p.is_file() and p.name != "manifest.json"}
        (args.out / "manifest.json").write_bytes(encoded({"schema_version": 1,
            "protocol_sha256": protocol_sha, "source_artifacts_sha256": ARTIFACT_SHA, "files": files}))


if __name__ == "__main__":
    main()
