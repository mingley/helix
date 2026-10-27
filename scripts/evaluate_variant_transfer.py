#!/usr/bin/env python3
"""Evaluate frozen signed-variant predictions from a sealed normalized handoff.

This component neither reads source spreadsheets nor fits/calibrates a model.
It verifies pinned producer artifacts, then evaluates every admitted observation.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import struct
import sys

CONTEXTS = ["K562", "HepG2"]
MODELS = ["zero", "I", "S"]
REPLICATES = 2000
SEED = 0x2026102801
PROBABILITIES = [0.025, 0.975]
MAX_DRAW_INDICES = 50_000_000
MASK = (1 << 64) - 1
HEX = re.compile(r"[0-9a-f]{64}\Z")
ANALYSIS_CONTRACT = {
    "contexts": CONTEXTS, "primary_weighting": "equal_component",
    "baseline": "zero", "candidate": "I",
    "bootstrap": {"replicates": REPLICATES, "rng": "SplitMix64", "seed": SEED,
                  "quantiles": PROBABILITIES, "quantile_type": 7},
    "minimum_variants": 100, "minimum_components": 10,
}


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


def load_json(path):
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle, object_pairs_hook=strict_object,
                         parse_constant=lambda value: (_ for _ in ()).throw(
                             ValueError(f"nonfinite JSON literal: {value}")))


def encoded(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def is_number(value):
    return type(value) in (int, float) and math.isfinite(value)


def bits(value):
    return struct.pack("<d", value).hex()


def average(values):
    require(bool(values), "cannot average an empty set")
    result = math.fsum(value / len(values) for value in values)
    require(math.isfinite(result), "nonfinite arithmetic result")
    return result


def valid_label(value):
    return isinstance(value, str) and value and value == value.strip()


class SplitMix64:
    """Same uint64 stream and modulo-index rule as the frozen existing toolkit."""
    def __init__(self, seed):
        require(type(seed) is int and 0 <= seed <= MASK, "seed must be uint64")
        self.state = seed

    def next_u64(self):
        self.state = (self.state + 0x9E3779B97F4A7C15) & MASK
        value = self.state
        value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK
        value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK
        return value ^ (value >> 31)


def quantile(values, probability):
    require(bool(values) and 0 <= probability <= 1 and all(is_number(v) for v in values),
            "invalid Type7 quantile inputs")
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    low, high = math.floor(position), math.ceil(position)
    fraction = position - low
    return ordered[low] * (1 - fraction) + ordered[high] * fraction


def validate_variants(payload):
    require(isinstance(payload, dict) and set(payload) == {
        "schema_version", "protocol_sha256", "binding_sha256", "prediction_I_sha256",
        "prediction_S_sha256", "variants"} and type(payload["schema_version"]) is int
            and payload["schema_version"] == 1, "unexpected normalized input schema")
    for field in ["protocol_sha256", "binding_sha256", "prediction_I_sha256", "prediction_S_sha256"]:
        require(isinstance(payload[field], str) and bool(HEX.fullmatch(payload[field])),
                f"invalid input commitment {field}")
    variants = payload["variants"]
    require(isinstance(variants, list) and 0 < len(variants) <= 1_000_000,
            "expected 1..1000000 admitted variants")
    seen, locus_components = set(), {}
    for row in variants:
        require(isinstance(row, dict) and set(row) == {
            "id", "component", "locus", "gc_preserving", "observed", "predictions"},
            "unexpected variant fields")
        require(all(valid_label(row[key]) for key in ["id", "component", "locus"]),
                "invalid variant/component/locus label")
        require(row["id"] not in seen, "duplicate variant ID")
        seen.add(row["id"])
        require(row["locus"] not in locus_components
                or locus_components[row["locus"]] == row["component"],
                "one exact source locus label is split across components")
        locus_components[row["locus"]] = row["component"]
        require(type(row["gc_preserving"]) is bool, "GC-preserving metadata must be boolean")
        require(isinstance(row["observed"], dict) and set(row["observed"]) == set(CONTEXTS),
                "missing or extra observed cell context")
        require(isinstance(row["predictions"], dict) and set(row["predictions"]) == {"I", "S"},
                "missing or extra prediction model")
        for values in [row["observed"], *row["predictions"].values()]:
            require(isinstance(values, dict) and set(values) == set(CONTEXTS)
                    and all(is_number(value) for value in values.values()),
                    "outcomes/predictions require ordinary finite numeric values in both cells")
    return sorted(variants, key=lambda row: row["id"])


def sign(value):
    return "positive" if value > 0 else "negative" if value < 0 else "zero"


def association(observed, predicted):
    require(len(observed) == len(predicted) and bool(observed), "association needs paired values")
    signs = ["negative", "zero", "positive"]
    confusion = {a: {b: 0 for b in signs} for a in signs}
    for actual, forecast in zip(observed, predicted):
        confusion[sign(actual)][sign(forecast)] += 1
    matching = sum(confusion[s][s] for s in signs)
    xmean, ymean = average(observed), average(predicted)
    dx = [value - xmean for value in observed]
    dy = [value - ymean for value in predicted]
    # Normalize first to avoid overflow in sums of products. Constant vectors
    # are explicitly nonestimable; no near-zero threshold or filtering is used.
    xscale, yscale = max(map(abs, dx)), max(map(abs, dy))
    # Summing rounded v/n terms can move a constant vector's mean by one ULP.
    # Detect exact equality on the original values before relying on centering.
    constant_input = (all(value == observed[0] for value in observed)
                      or all(value == predicted[0] for value in predicted))
    if constant_input or xscale == 0 or yscale == 0:
        correlation = {"estimable": False, "value": None, "reason": "constant observed or predicted vector"}
    else:
        xn, yn = [v / xscale for v in dx], [v / yscale for v in dy]
        value = math.fsum(a * b for a, b in zip(xn, yn)) / math.sqrt(
            math.fsum(a * a for a in xn) * math.fsum(b * b for b in yn))
        require(math.isfinite(value), "nonfinite Pearson result")
        correlation = {"estimable": True, "value": value, "reason": None}
    return {"observations": len(observed), "pearson": correlation,
            "exact_sign_agreement_including_zero": matching / len(observed),
            "sign_agreement_count": matching, "sign_confusion_observed_then_predicted": confusion,
            "exact_zero_observations": sum(v == 0 for v in observed),
            "exact_zero_predictions": sum(v == 0 for v in predicted)}


def losses(row):
    result = {}
    for cell in CONTEXTS:
        target = row["observed"][cell]
        result[cell] = {}
        for model in MODELS:
            prediction = 0.0 if model == "zero" else row["predictions"][model][cell]
            residual = prediction - target
            squared = residual * residual
            require(math.isfinite(squared), "loss overflow")
            result[cell][model] = {"mse": squared, "mae": abs(residual)}
    result["paired"] = {model: {metric: average([result[cell][model][metric] for cell in CONTEXTS])
                               for metric in ["mse", "mae"]} for model in MODELS}
    return result


def aggregate(rows):
    """Average variant losses; rows here are each one whole paired variant."""
    return {scope: {model: {metric: average([row["losses"][scope][model][metric] for row in rows])
                            for metric in ["mse", "mae"]} for model in MODELS}
            for scope in ["paired", *CONTEXTS]}


def weighted_component_average(components, indices, weighting):
    weights = [1 if weighting == "equal_component" else components[i]["variant_count"] for i in indices]
    denominator = sum(weights)
    return {scope: {model: {metric: math.fsum(
        components[i]["losses"][scope][model][metric] * weight / denominator
        for i, weight in zip(indices, weights)) for metric in ["mse", "mae"]} for model in MODELS}
        for scope in ["paired", *CONTEXTS]}


def gains(summary):
    return {scope: {model: {metric: summary[scope]["zero"][metric] - summary[scope][model][metric]
                            for metric in ["mse", "mae"]} for model in ["I", "S"]}
            for scope in ["paired", *CONTEXTS]}


def bootstrap_components(components, replicates=REPLICATES, seed=SEED):
    require(bool(components) and 1 <= replicates <= 100000, "invalid bootstrap size")
    require(len(components) * replicates <= MAX_DRAW_INDICES,
            "retained draws exceed 50-million-index budget; no components dropped")
    require([row["component"] for row in components] == sorted({row["component"] for row in components}),
            "bootstrap components must be unique and lexicographically sorted")
    rng = SplitMix64(seed)
    draws, statistics = [], []
    for _ in range(replicates):
        selected = [rng.next_u64() % len(components) for _ in components]
        draws.append(selected)
        record = {}
        for weighting in ["equal_component", "equal_variant"]:
            summary = weighted_component_average(components, selected, weighting)
            record[weighting] = {"losses": summary, "gains_vs_zero": gains(summary)}
        statistics.append(record)
    intervals = {}
    for weighting in ["equal_component", "equal_variant"]:
        intervals[weighting] = {}
        for kind, models in [("losses", MODELS), ("gains_vs_zero", ["I", "S"])]:
            intervals[weighting][kind] = {scope: {model: {metric: [
                quantile([row[weighting][kind][scope][model][metric] for row in statistics], p)
                for p in PROBABILITIES] for metric in ["mse", "mae"]} for model in models}
                for scope in ["paired", *CONTEXTS]}
    return {"component_order": [row["component"] for row in components],
            "replicates": replicates, "seed": seed, "rng": "SplitMix64",
            "index_rule": "next_u64 modulo component_count", "interval_probabilities": PROBABILITIES,
            "quantile_method": "Type7", "resampling_indices": draws,
            "sampled_statistics": statistics, "intervals": intervals}


def evaluate(variants, replicates=REPLICATES, seed=SEED):
    ledger = [{"id": row["id"], "component": row["component"], "locus": row["locus"],
               "gc_preserving": row["gc_preserving"], "losses": losses(row),
               "observed_float64_little_endian": {cell: bits(row["observed"][cell]) for cell in CONTEXTS},
               "prediction_float64_little_endian": {model: {cell: bits(row["predictions"][model][cell])
                                                            for cell in CONTEXTS} for model in ["I", "S"]}}
              for row in variants]
    groups = {}
    for row in ledger:
        groups.setdefault(row["component"], []).append(row)
    components = [{"component": group, "variant_count": len(rows),
                   "variant_ids": [row["id"] for row in rows],
                   "loci": sorted({row["locus"] for row in rows}), "losses": aggregate(rows)}
                  for group, rows in sorted(groups.items())]
    indices = list(range(len(components)))
    observed = {}
    for weighting in ["equal_component", "equal_variant"]:
        summary = weighted_component_average(components, indices, weighting)
        observed[weighting] = {"losses": summary, "gains_vs_zero": gains(summary)}
    associations = {}
    for scope in ["paired", *CONTEXTS]:
        cells = CONTEXTS if scope == "paired" else [scope]
        actual = [row["observed"][cell] for row in variants for cell in cells]
        associations[scope] = {model: association(actual, [
            0.0 if model == "zero" else row["predictions"][model][cell]
            for row in variants for cell in cells]) for model in MODELS}
    gc_diagnostics = {}
    for value, label in [(True, "GC_preserving"), (False, "GC_changing")]:
        selected = [row for row in ledger if row["gc_preserving"] is value]
        represented = sorted({row["component"] for row in selected})
        subset_groups = {}
        for row in selected:
            subset_groups.setdefault(row["component"], []).append(row)
        subset_losses = [aggregate(subset_groups[group]) for group in represented]
        gc_diagnostics[label] = {"variant_ids": [row["id"] for row in selected],
            "represented_component_ids": represented,
            "empty_original_component_ids": sorted(set(groups) - set(represented)),
            "variants": len(selected), "estimable": bool(selected),
            "equal_variant_losses": aggregate(selected) if selected else None,
            "equal_represented_component_losses": {scope: {model: {metric: average([
                summary[scope][model][metric] for summary in subset_losses])
                for metric in ["mse", "mae"]} for model in MODELS}
                for scope in ["paired", *CONTEXTS]} if selected else None}
    locus_groups = {}
    for row in ledger:
        locus_groups.setdefault(row["locus"], []).append(row)
    loci = [{"locus": locus, "variant_ids": [row["id"] for row in rows],
             "component_ids": sorted({row["component"] for row in rows}),
             "variant_count": len(rows), "equal_variant_losses": aggregate(rows)}
            for locus, rows in sorted(locus_groups.items())]
    resampling = bootstrap_components(components, replicates, seed)
    gain = observed["equal_component"]["gains_vs_zero"]["paired"]["I"]["mse"]
    interval = resampling["intervals"]["equal_component"]["gains_vs_zero"]["paired"]["I"]["mse"]
    support = len(variants) >= 100 and len(components) >= 10
    return {"schema_version": 1, "variant_order": [row["id"] for row in variants],
            "contexts": CONTEXTS, "model_order": MODELS, "variant_count": len(variants),
            "component_count": len(components), "observed": observed,
            "primary": {"metric": "equal_component_macroMSE_zero_minus_I",
                        "observed_improvement": gain, "conditional_95_percent_range": interval,
                        "minimum_support_passed": support, "positive_observed_gain": gain > 0,
                        "positive_lower_endpoint": interval[0] > 0,
                        "useful_transfer_rule_passed": support and gain > 0 and interval[0] > 0},
            "associations_equal_observation_weight": associations,
            "gc_diagnostics_descriptive_only": gc_diagnostics, "loci": loci,
            "per_variant": ledger, "per_component": components, "resampling": resampling,
            "limits": ["Fixed-fit conditional selected-panel component percentile range; no refitting.",
                       "Source and cohort gates are a sealed producer handoff, not independently reconstructed here.",
                       "Exact zero signs are retained; correlation is nonestimable for constant vectors.",
                       "GC subgroup macro summaries average represented original components; empty components are listed.",
                       "Secondary metrics do not replace the frozen primary comparison.",
                       "Numerical validity and biological predictive performance are separate flags."]}


def check_protocol(protocol):
    require(protocol.get("protocol_id") == "variant-transfer-v1"
            and protocol.get("status") == "frozen",
            "protocol must be variant-transfer-v1 with frozen status")
    require(protocol.get("contexts", {}).get("primary") == CONTEXTS
            and protocol.get("analysis_contract") == ANALYSIS_CONTRACT,
            "protocol analysis contract differs from implemented fixed analysis")


def verify_prediction_directory(root, model, pins, variants, protocol):
    root = Path(root)
    require(root.is_dir() and not root.is_symlink(), "prediction directory missing or symlinked")
    manifest_path, predictions_path = root / "manifest.json", root / "predictions.json"
    require(not manifest_path.is_symlink() and not predictions_path.is_symlink(), "symlinked predictor file")
    require(sha256(manifest_path) == pins[f"prediction_{model}_manifest_sha256"],
            "prediction manifest pin mismatch")
    manifest = load_json(manifest_path)
    require(manifest.get("schema_version") == 1 and isinstance(manifest.get("files"), dict)
            and set(manifest["files"]) == {"predictions.json"}, "unexpected predictor manifest schema")
    require(sha256(predictions_path) == manifest["files"]["predictions.json"]
            == pins[f"prediction_{model}_sha256"], "prediction file pin mismatch")
    prediction = load_json(predictions_path)
    require(prediction.get("schema_version") == 1 and prediction.get("model") == model
            and prediction.get("endpoint_direction") == "alternate minus reference",
            "wrong model or effect direction")
    training = manifest.get("training")
    require(isinstance(training, dict)
            and training.get("artifacts_sha256") == protocol["pretrained_model"]["source_artifacts_sha256"]
            and training.get("protocol_sha256") == protocol["pretrained_model"]["source_protocol_sha256"],
            "predictor training provenance differs from frozen protocol")
    require(prediction.get("commitments") == training | {
        "pairs_sha256": manifest.get("pairs_sha256")}, "predictor receipt cross-binding differs")
    lookup = {}
    for row in prediction.get("pairs", []):
        require(isinstance(row, dict) and valid_label(row.get("id"))
                and row["id"] not in lookup and row.get("context") in CONTEXTS
                and is_number(row.get("predicted_alternate_minus_reference")),
                "invalid or duplicate predictor row")
        lookup[row["id"]] = row
    expected = {row["id"] + "::" + context for row in variants for context in CONTEXTS}
    require(set(lookup) == expected, "predictor pair cohort differs from complete normalized cohort")
    for row in variants:
        for cell in CONTEXTS:
            prediction_row = lookup[row["id"] + "::" + cell]
            require(prediction_row["context"] == cell and bits(row["predictions"][model][cell])
                    == bits(prediction_row["predicted_alternate_minus_reference"]),
                    "normalized/predictor value or cell differs")
    return {"training": training, "pairs_sha256": manifest.get("pairs_sha256")}


def verify_handoff(args):
    pins = load_json(args.pins)
    fields = {"input_sha256", "protocol_sha256", "binding_sha256", "prediction_I_sha256",
              "prediction_S_sha256", "prediction_I_manifest_sha256", "prediction_S_manifest_sha256"}
    require(isinstance(pins, dict) and set(pins) == fields | {"schema_version"}
            and type(pins["schema_version"]) is int and pins["schema_version"] == 1
            and all(isinstance(pins[key], str) and HEX.fullmatch(pins[key]) for key in fields),
            "unexpected pins schema")
    for key, path in [("input_sha256", args.input), ("protocol_sha256", args.protocol),
                      ("binding_sha256", args.binding)]:
        require(not path.is_symlink() and sha256(path) == pins[key], f"{key} mismatch")
    protocol, binding, payload = load_json(args.protocol), load_json(args.binding), load_json(args.input)
    check_protocol(protocol)
    require(binding.get("protocol_id") == "variant-transfer-v1"
            and binding.get("protocol_sha256") == pins["protocol_sha256"]
            and binding.get("source_cohort_gate", {}).get("passed") is True,
            "sealed source/cohort handoff gate not passed")
    variants = validate_variants(payload)
    require(all(payload[key] == pins[key] for key in ["protocol_sha256", "binding_sha256",
            "prediction_I_sha256", "prediction_S_sha256"]), "normalized commitment differs from pins")
    producer_provenance = [verify_prediction_directory(root, model, pins, variants, protocol)
                          for model, root in [("I", args.prediction_I), ("S", args.prediction_S)]]
    require(producer_provenance[0] == producer_provenance[1],
            "I and S producer DNA/model provenance differs")
    return variants, pins


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ["input", "protocol", "binding", "pins", "prediction-I", "prediction-S", "out"]:
        parser.add_argument("--" + option, type=Path, required=True)
    args = parser.parse_args(argv)
    require(not args.out.exists(), "output path already exists")
    variants, pins = verify_handoff(args)
    result = evaluate(variants)
    args.out.mkdir(parents=True)
    with (args.out / "evaluation.json.gz").open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(encoded(result))
    (args.out / "summary.json").write_bytes(encoded({key: value for key, value in result.items()
        if key not in {"per_variant", "resampling"}}))
    for name, source in [("normalized-input.json", args.input), ("protocol.json", args.protocol),
                         ("binding.json", args.binding), ("pins.json", args.pins)]:
        shutil.copyfile(source, args.out / name)
    shutil.copyfile(Path(__file__), args.out / "implementation.py")
    manifest = {"schema_version": 1, "input_pins": pins,
                "files": {path.name: sha256(path) for path in sorted(args.out.iterdir())}}
    (args.out / "manifest.json").write_bytes(encoded(manifest))
    print(json.dumps({"variants": result["variant_count"], "components": result["component_count"],
                      "primary": result["primary"], "manifest_sha256": sha256(args.out / "manifest.json")},
                     sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, TypeError, KeyError, OSError, OverflowError, json.JSONDecodeError) as error:
        print(f"variant evaluation failed: {error}", file=sys.stderr)
        raise SystemExit(2)
