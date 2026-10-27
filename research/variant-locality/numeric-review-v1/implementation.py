#!/usr/bin/env python3
"""Independent frozen-model numerical audit, with no predictor imports or refit."""
import argparse
import collections
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys

DNA = "ACGT"
TRANS = str.maketrans("ACGT", "TGCA")
PIN = "2ea476bc5413dc13486dec9f3d7e9162a50a0e47f1d83ce29efabc5f5cc5af9d"
CELLS = ["HepG2", "K562", "WTC11"]


def rc(s):
    return s.translate(TRANS)[::-1]


KMERS = sorted({min("".join(t), rc("".join(t))) for t in itertools.product(DNA, repeat=3)})


def cnt(s):
    counts = collections.Counter(min(s[i:i + 3], rc(s[i:i + 3])) for i in range(len(s) - 2))
    return tuple(counts[k] for k in KMERS) + (s.count("G") + s.count("C"),)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(p, value):
    p.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def all_cases():
    cases = []
    for chars in itertools.product(DNA, repeat=5):
        word = "".join(chars)
        for alt in DNA:
            if alt != word[2]:
                reference = "A" * 97 + word + "C" * 98
                mutant = reference[:99] + alt + reference[100:]
                cases.append(("interior", word + ">" + alt, 99, reference, mutant))
    for position, size in ((0, 3), (1, 4)):
        for chars in itertools.product(DNA, repeat=size):
            word = "".join(chars)
            for alt in DNA:
                if alt == word[position]:
                    continue
                reference = word + "A" * (200 - size)
                mutant = reference[:position] + alt + reference[position + 1:]
                name = str(position) + ":" + word + ">" + alt
                cases.append(("boundary", name, position, reference, mutant))
                cases.append(("boundary", "rc:" + name, 199 - position, rc(reference), rc(mutant)))
    return cases


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if a.out.exists():
        raise SystemExit("Preserve existing attempt")
    core = a.root / "research/results/context-activity-v3/core"
    assert sha(core / "artifacts.json") == PIN
    manifest = json.loads((core / "artifacts.json").read_bytes())
    for name, expected in manifest["files"].items():
        assert sha(core / name) == expected
    with gzip.open(core / "study.json.gz", "rt") as handle:
        study = json.load(handle)
    assert len(study["folds"]) == 5
    assert study["contexts"] == CELLS
    assert all(f["fit"]["feature_names"] == ["rc3mer_" + k for k in KMERS] + ["gc_fraction"]
               for f in study["folds"])
    fits = {"S": [f["fit"]["shared_dna_fit"] for f in study["folds"]]}
    fits.update({cell: [f["fit"]["independent_context_fits"][cell] for f in study["folds"]]
                 for cell in CELLS})
    # The exact interpretation of each persisted IEEE binary64 coefficient and
    # scale is its exact rational value, not a rounded decimal approximation.
    weights = {name: [[Fraction.from_float(float(b)) / Fraction.from_float(float(s))
                     / (198 if j < 32 else 200) for j, (b, s) in enumerate(zip(
                         fit["standardized_coefficients"], fit["training_feature_scales"]))]
                      for fit in group] for name, group in fits.items()}
    cases = all_cases()
    assert len(cases) == 4992
    deltas = {}
    for kind, name, pos, reference, mutant in cases:
        assert len(reference) == len(mutant) == 200
        assert [i for i in range(200) if reference[i] != mutant[i]] == [pos]
        deltas[name] = tuple(b - c for b, c in zip(cnt(mutant), cnt(reference)))
    a.out.mkdir(parents=True)
    overall = []
    bounds = {}
    for model in ("I", "S", "A"):
        cells = CELLS if model == "I" else ["K562"]
        pairs = [{"id": cell + "/" + name, "context": cell, "position_0based": pos,
                  "reference": reference, "alternate": mutant}
                 for cell in cells for kind, name, pos, reference, mutant in cases]
        infile = a.out / (model + "-pairs.json")
        dump(infile, {"schema_version": 1, "pairs": pairs})
        subprocess.run([sys.executable, str(a.root / "scripts/predict_variant_effects.py"),
            "--model-artifacts", str(core), "--expected-artifacts-sha256", PIN,
            "--pairs", str(infile), "--model", model, "--out", str(a.out / (model + "-primary"))],
            check=True)
        actual = json.loads((a.out / (model + "-primary/predictions.json")).read_bytes())
        assert len(actual["pairs"]) == len(pairs)
        max_ensemble_error = 0.0
        max_fold_error = 0.0
        max_direct_error = 0.0
        exact_by_cell = {cell: {"interior": [], "boundary": []} for cell in cells}
        actual_by_cell = {cell: {"interior": [], "boundary": []} for cell in cells}
        for output, pair in zip(actual["pairs"], pairs):
            assert output["id"] == pair["id"]
            cell, name = pair["id"].split("/", 1)
            vector = deltas[name]
            selected = weights[cell if model == "I" else "S"]
            exact_folds = [sum((weight * number for weight, number in zip(fold, vector) if number),
                               Fraction(0)) for fold in selected]
            exact_mean = sum(exact_folds, Fraction(0)) / 5
            error = abs(float(exact_mean) - output["predicted_alternate_minus_reference"])
            max_ensemble_error = max(error, max_ensemble_error)
            for exact, value in zip(exact_folds, output["fold_predictions"]):
                max_fold_error = max(max_fold_error, abs(float(exact) - value))
            max_direct_error = max(max_direct_error, output["max_abs_direct_subtraction_check_error"])
            assert error < 1e-12
            if not any(vector):
                assert exact_mean == output["predicted_alternate_minus_reference"] == 0
            category = "interior" if pair["position_0based"] == 99 else "boundary"
            exact_by_cell[cell][category].append(exact_mean)
            actual_by_cell[cell][category].append(output["predicted_alternate_minus_reference"])
        for cell in cells:
            bounds[model + "/" + cell] = {}
            for category, values in exact_by_cell[cell].items():
                lo, hi = min(values), max(values)
                assert lo == -hi
                bounds[model + "/" + cell][category] = {
                    "cases": len(values), "exact_distinct_predictions": len(set(values)),
                    "binary64_distinct_predictions": len(set(actual_by_cell[cell][category])),
                    "exact_zero_predictions": values.count(Fraction(0)),
                    "minimum_float": float(lo), "maximum_float": float(hi),
                    "minimum_rational": str(lo), "maximum_rational": str(hi)}
        overall.append({"model": model, "pairs": len(pairs),
            "max_abs_ensemble_error_vs_exact_rational_binary64_coefficients": max_ensemble_error,
            "max_abs_fold_error_vs_exact_rational_binary64_coefficients": max_fold_error,
            "max_abs_primary_direct_subtraction_error": max_direct_error})
    s = json.loads((a.out / "S-primary/predictions.json").read_bytes())["pairs"]
    additive = json.loads((a.out / "A-primary/predictions.json").read_bytes())["pairs"]
    assert [v["predicted_alternate_minus_reference"] for v in s] == [
        v["predicted_alternate_minus_reference"] for v in additive]
    summary = {"primary_implementation_imported": False,
        "experimental_external_outcomes_used": False, "models_refitted": False,
        "training_artifacts_sha256": PIN, "primary_implementation_sha256": sha(
            a.root / "scripts/predict_variant_effects.py"),
        "cases_per_context": {"interior": 3072, "boundary": 1920},
        "comparisons": overall, "bounds": bounds,
        "limits": ["Exact rational calculations interpret finite saved binary64 coefficients/scales mathematically; predictor operations still round.",
            "Bounds cover all 200nt single substitutions under frozen coefficients; they are predicted model-scale bounds, not observed or clinical effects.",
            "Boundary duplicates include both orientations/ends intentionally.",
            "No proof of IEEE754 exact mixed-second-difference zero is implied."]}
    dump(a.out / "summary.json", summary)
    (a.out / "implementation.py").write_bytes(Path(__file__).read_bytes())
    dump(a.out / "manifest.json", {"files": {str(p.relative_to(a.out)): {
        "sha256": sha(p), "bytes": p.stat().st_size} for p in sorted(a.out.rglob("*")) if p.is_file()}})
    print(json.dumps({"comparisons": overall,
        "bounds_float": {k: {c: {field: v for field, v in record.items() if "rational" not in field}
                              for c, record in group.items()} for k, group in bounds.items()}}, indent=2))


if __name__ == "__main__":
    main()
