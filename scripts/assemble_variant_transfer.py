#!/usr/bin/env python3
"""Join authenticated observed effects and fixed predictions without arithmetic."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import struct

CELLS = ("K562", "HepG2")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def encoded(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def bits(value):
    require(type(value) in (int, float) and math.isfinite(value), "ordinary finite numeric value required")
    return struct.pack("<d", value).hex()


def verify(root):
    root = Path(root)
    require(root.is_dir() and not root.is_symlink(), "unsafe producer directory")
    manifest = load(root / "manifest.json")
    require(set(manifest["files"]) | {"manifest.json"} == {p.name for p in root.iterdir()},
            "producer inventory differs")
    for name, receipt in manifest["files"].items():
        path = root / name
        require(Path(name).name == name and path.is_file() and not path.is_symlink(), "unsafe producer member")
        expected = receipt if isinstance(receipt, str) else receipt["sha256"]
        require(sha(path) == expected, "producer member hash differs: " + name)
        if isinstance(receipt, dict):
            require(path.stat().st_size == receipt["bytes"], "producer member size differs")
    return manifest


def join(observations, predictions):
    rows, seen, complete_ids = [], set(), set()
    for row in observations["variants"]:
        require(row["id"] not in seen and type(row["complete_case"]) is bool, "duplicate or invalid observation")
        seen.add(row["id"])
        if not row["complete_case"]:
            continue
        require(set(row["observed"]) == set(CELLS), "paired cell observations required")
        for cell in CELLS:
            require(bits(row["observed"][cell]) == row["source_float64_le_hex"][cell], "source float bits changed")
        values = {}
        for model, lookup in predictions.items():
            values[model] = {}
            for cell in CELLS:
                identity = row["id"] + "::" + cell
                require(identity in lookup and lookup[identity]["context"] == cell, "prediction ID/cell differs")
                forecast = lookup[identity]["predicted_alternate_minus_reference"]
                bits(forecast)
                values[model][cell] = forecast
                complete_ids.add(identity)
        rows.append({key: row[key] for key in ("id", "component", "locus", "gc_preserving")} |
                    {"observed": row["observed"], "predictions": values})
    require(rows and set(predictions) == {"I", "S"}, "both fixed prediction models required")
    require(all(set(lookup) == complete_ids for lookup in predictions.values()), "extra prediction pairs")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("observations", "prediction-I", "prediction-S", "protocol", "binding", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    require(not args.out.exists(), "preserve previous handoff attempts")
    observe_manifest = verify(args.observations)
    binding, protocol = load(args.binding), load(args.protocol)
    require(binding["source_cohort_gate"]["passed"] is True, "metadata source gate failed")
    require(binding["protocol_id"] == protocol["protocol_id"]
            and binding["protocol_sha256"] == sha(args.protocol), "binding authority differs")
    require(observe_manifest["protocol_sha256"] == sha(args.protocol)
            and observe_manifest["binding_sha256"] == sha(args.binding), "observation authority differs")
    summary = load(args.observations / "summary.json")
    require(summary["source_observation_gate_passed"] is True, "source value gate failed")
    observations = load(args.observations / "observations.json.gz")
    predictions, pins = {}, {"schema_version": 1, "protocol_sha256": sha(args.protocol),
                             "binding_sha256": sha(args.binding)}
    for model, directory in (("I", args.prediction_I), ("S", args.prediction_S)):
        manifest = verify(directory)
        require(manifest["pairs_sha256"] == sha(args.observations / "completecasepairs.json"),
                "predictor DNA cohort differs")
        require(manifest["training"]["artifacts_sha256"] == binding["model_artifacts_sha256"],
                "frozen model differs")
        payload = load(directory / "predictions.json")
        require(payload["model"] == model and payload["endpoint_direction"] == "alternate minus reference",
                "prediction model/direction differs")
        lookup = {row["id"]: row for row in payload["pairs"]}
        require(len(lookup) == len(payload["pairs"]), "duplicate prediction ID")
        predictions[model] = lookup
        pins[f"prediction_{model}_sha256"] = sha(directory / "predictions.json")
        pins[f"prediction_{model}_manifest_sha256"] = sha(directory / "manifest.json")
    variants = join(observations, predictions)
    require(len(variants) == summary["complete_case_variants"], "complete-case count differs")
    payload = {key: pins[key] for key in ("protocol_sha256", "binding_sha256",
               "prediction_I_sha256", "prediction_S_sha256")} | {"schema_version": 1, "variants": variants}
    args.out.mkdir(parents=True)
    (args.out / "normalized-input.json").write_bytes(encoded(payload))
    pins["input_sha256"] = sha(args.out / "normalized-input.json")
    (args.out / "pins.json").write_bytes(encoded(pins))
    shutil.copyfile(__file__, args.out / "implementation.py")
    (args.out / "manifest.json").write_bytes(encoded({"schema_version": 1,
        "source_observation_manifest_sha256": sha(args.observations / "manifest.json"),
        "files": {p.name: sha(p) for p in sorted(args.out.iterdir())}}))
    print(json.dumps({"variants": len(variants), "input_sha256": pins["input_sha256"]}))


if __name__ == "__main__":
    main()
