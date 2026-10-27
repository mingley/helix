#!/usr/bin/env python3
"""Retain and verify exact variant-transfer evidence with the standard library."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil

V1 = "811930c2d3a5ee763014ea611acda4d2c84404b7d5ae399c7843c62d0a828149"
V1_BINDING = "9c26b3fc628f16bd6bb677dee8eed01e509b6102d44f95c0037ae2f55676ba79"
V2 = "26861a28a1c049bcb6da66d17ff25bb13143f45cbf166a224af2a5aa27f2e04c"
V2_BINDING = "a54badc4645d7ac88c31b7722fbbd89aaefa2b995bc0a02befd15123d83d4fba"
V2_NAMES = ("observations", "I", "S", "handoff", "evaluation", "independent", "comparison", "figures")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def files(root, excluded=()):
    require(root.is_dir() and not root.is_symlink(), "unsafe evidence root")
    result = {}
    for path in sorted(root.rglob("*")):
        require(not path.is_symlink(), "symlink in evidence")
        if path.is_file():
            name = path.relative_to(root).as_posix()
            if name not in excluded:
                require(path.stat().st_size < 100_000_000, "file exceeds Git publication limit")
                result[name] = {"bytes": path.stat().st_size, "sha256": sha(path)}
    return result


def checked(root, manifest="manifest.json"):
    document = load(root / manifest)
    actual = files(root, (manifest,))
    require(set(document["files"]) == set(actual), "evidence inventory differs: " + str(root))
    for name, receipt in document["files"].items():
        part = PurePosixPath(name)
        require(not part.is_absolute() and ".." not in part.parts and "\\" not in name,
                "unsafe evidence member")
        require(actual[name] == receipt if isinstance(receipt, dict)
                else actual[name]["sha256"] == receipt, "evidence SHA/size differs: " + name)
    return document


def failed_status(primary, independent, comparison):
    for root in (primary, independent, comparison):
        checked(root)
    require((primary / "failure.json").is_file() and (independent / "failure.json").is_file(),
            "both source failures must be retained")
    require(all(not (root / "evaluation.json.gz").exists() for root in (primary, independent)),
            "failed source publication must not contain a performance evaluation")
    require(sha(primary / "protocol.json") == sha(independent / "protocol.json") == V1
            and sha(primary / "binding.json") == sha(independent / "binding.json") == V1_BINDING,
            "v1 frozen authority differs")
    record = load(comparison / "comparison.json")
    require(record["protocol_sha256"] == V1 and record["binding_sha256"] == V1_BINDING
            and record["source_gate_failure_reproduced"] is True
            and record["source_observation_gate_passed"] is False
            and record["failure_count"] == 0 and record["predictive_performance_evaluated"] is False
            and record["bootstrap_computed"] is False and record["performance_result"] is None,
            "v1 failure comparison/status differs")
    require(record["primary_manifest"]["sha256"] == sha(primary / "manifest.json")
            and record["independent_manifest"]["sha256"] == sha(independent / "manifest.json"),
            "failure comparison identifies other runs")
    require(record["selected_cell_comparisons"] == 241
            and record["ordinary_numeric_float_bits_compared_exactly"] == 240
            and record["first_failure"]["literal_text"] == "NA"
            and record["first_failure"]["coordinate"] == "AO42972", "source failure evidence differs")
    return {"kind": "variant-transfer-v1-source-failure", "protocol_sha256": V1,
            "binding_sha256": V1_BINDING, "second_freeze_commit": "a21fd6f44f493b20cddc51e6cde7c7ae14e490e4",
            "source_observation_gate_passed": False, "failure_reproduced": True,
            "transfer_performance_computed": False, "selected_cell_comparisons": 241,
            "exact_numeric_float_cells": 240, "fully_decoded_paired_variants": 120,
            "first_failure": record["first_failure"], "independent_fixed_predictions_retained": True,
            "gate_comparison_sha256": sha(comparison / "comparison.json"),
            "limits": ["Reproduced source-gate failure, not a model-accuracy result.",
                       "Independent fixed predictions exist; no target-performance/bootstrap was computed.",
                       "Openpyxl can instantiate nonselected numeric cells internally."]}


def completed_status(roots):
    for name in V2_NAMES[:-1]:
        checked(roots[name])
    evaluation, independent = roots["evaluation"], roots["independent"]
    summary = load(evaluation / "summary.json")
    observations = load(roots["observations"] / "summary.json")
    protocol, binding = load(evaluation / "protocol.json"), load(evaluation / "binding.json")
    require(sha(evaluation / "protocol.json") == sha(independent / "protocol.json")
            == sha(roots["observations"] / "protocol.json") == V2
            and sha(evaluation / "binding.json") == sha(independent / "binding.json")
            == sha(roots["observations"] / "binding.json") == V2_BINDING, "v2 authority differs")
    require(protocol["analysis_classification"] == summary["analysis_classification"]
            == "exploratory_after_v1_source_failure"
            and protocol["operational_missingness"]["source_semantics_established"] is False
            and protocol["parent_protocol_sha256"] == V1, "exploratory classification was changed")
    require(observations["source_observation_gate_passed"] is True
            and load(independent / "status.json")["status"] == "complete", "numeric source/run gates incomplete")
    comparison = load(roots["comparison"] / "comparison.json")
    require(comparison["matched"] is True and comparison["failure_count"] == 0
            and comparison["protocol_sha256"] == V2 and comparison["failures"] == []
            and all(row["matched"] is True for row in comparison["checks"]), "numeric comparison did not pass")
    require(comparison["primary_manifest"]["sha256"] == sha(evaluation / "manifest.json")
            and comparison["independent_manifest"]["sha256"] == sha(independent / "manifest.json"),
            "comparison identifies other evaluation runs")
    raw = comparison["raw_source_and_fixed_prediction_comparison"]
    require(raw["matched"] is True and raw["failure_count"] == 0
            and comparison["normalized_handoff_comparison"]["matched"] is True, "source/prediction comparison failed")
    for key, name in (("primary_observations", "observations"), ("prediction_I", "I"), ("prediction_S", "S")):
        require(raw["parent_manifests"][key]["sha256"] == sha(roots[name] / "manifest.json"),
                "raw comparison identifies another producer")
    require(sha(roots["handoff"] / "normalized-input.json") == sha(evaluation / "normalized-input.json")
            and sha(roots["handoff"] / "pins.json") == sha(evaluation / "pins.json"), "normalized handoff differs")
    figure = load(roots["figures"] / "manifest.json")
    require(files(roots["figures"], ("manifest.json",)) == figure["outputs"], "figure inventory differs")
    require(figure["producer_manifest_sha256"] == sha(evaluation / "manifest.json")
            and figure["producer_files"] == load(evaluation / "manifest.json")["files"],
            "figures identify another evaluation")
    grouping = binding["source_cohort_gate"]["independent_grouping"]
    require(grouping["status"] == "matched" and grouping["failed_checks"] == [], "global grouping comparison failed")
    losses = summary["observed"]["equal_component"]["losses"]["paired"]
    require(summary["primary"]["observed_improvement"] == losses["zero"]["mse"] - losses["I"]["mse"],
            "primary gain differs from archived losses")
    return {"kind": "variant-transfer-v2-exploratory", "protocol_sha256": V2,
        "binding_sha256": V2_BINDING, "freeze_commit": "156c5676c4f79328153edc517a737a2d83c97a48",
        "analysis_classification": summary["analysis_classification"], "numeric_reproduction": True,
        "comparison_categories": len(comparison["checks"]), "global_grouping_checks": grouping["checks_passed"],
        "max_absolute_evaluation_error": comparison["max_absolute_numeric_error"],
        "exact_integer_feature_numerators": raw["exact_integer_feature_numerators_compared"],
        "fixed_predictions_compared": raw["fixed_prediction_numbers_compared"],
        "metadata_variants": observations["metadata_eligible_variants"],
        "complete_case_variants": summary["variant_count"], "components": summary["component_count"],
        "missing_variants": observations["missing_variants"], "literal_NA_cells": observations["missing_cells"],
        "macro_losses": losses, "primary": summary["primary"], "v1_source_gate_passed": False,
        "operational_numeric_source_gate_passed": True, "source_missingness_semantics_established": False,
        "full_protocol_independent_reproduction_claimed": comparison["full_protocol_independent_reproduction_claimed"],
        "limits": ["Exploratory complete-case analysis after the retained v1 source failure.",
                   "Literal NA is operationally unreported; its biological/producer reason remains unknown.",
                   "No refit or calibration on new outcomes; conditional intervals omit assay/refit uncertainty.",
                   "Global grouping was independently reconstructed; assembly/source-class preparation remains a bound handoff.",
                   "No independent endogenous mechanism, clinical interpretation or global novelty is established."]}


def verify(root):
    checked(root, "publication.json")
    published = load(root / "publication.json")
    if published["status"]["kind"] == "variant-transfer-v1-source-failure":
        status = failed_status(root / "primary-source", root / "independent-source", root / "gate-comparison")
    else:
        require(published["status"]["kind"] == "variant-transfer-v2-exploratory", "unknown publication kind")
        status = completed_status({name: root / name for name in V2_NAMES})
    require(load(root / "publication.json")["status"] == status, "publication status was changed")
    print(json.dumps({"publication_verified": True, **status}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--primary", type=Path)
    parser.add_argument("--independent", type=Path)
    parser.add_argument("--comparison", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--v2", action="store_true")
    parser.add_argument("--prefix", type=Path)
    args = parser.parse_args()
    if args.verify:
        verify(args.verify)
        return
    if args.v2:
        require(args.prefix and args.out and not args.out.exists(), "fresh prefix/output required")
        roots = {name: Path(str(args.prefix) + "-" + name) for name in V2_NAMES}
        status = completed_status(roots)
        args.out.mkdir(parents=True)
        for name, source in roots.items():
            shutil.copytree(source, args.out / name)
        shutil.copyfile(__file__, args.out / "packaging-implementation.py")
        (args.out / "README.md").write_text('''# Exploratory external variant-effect transfer

The fixed five-fit model does not beat zero-effect prediction on 2,153 complete
variants in 172 untouched components. Macro MSE is 0.2033556 versus 0.2033201;
zero-minus-model gain is -0.0000354532, conditional 95% range
[-0.00144180, 0.00140062]. Useful-transfer rule: **false**.

This is explicitly exploratory after the preserved v1 source-gate failure.
Exactly 20 literal NA cells in 19 variants are operationally unreported, with no
claim about why they are missing. No new-outcome fitting or calibration occurred.
The metadata cohort, source/model/context/grouping/metrics and resampling remain
fixed. V1 is retained separately; its original source gate remains failed.

All 17 independent evaluation comparison categories pass, including 852,588 exact
feature numerators, all 51,672 fold/ensemble predictions and 344,000 bootstrap
indices. Of 4,344 selected source fields, 4,324 are bit-exact numeric values and 20 are
the literal NA/missing records. Full global grouping separately passes 4,689 checks.
Raw assembly/source-class preparation is an authenticated handoff; the numerical
comparison does not claim full independent reconstruction of every source step.

Every original source ledger, missingness record, prediction, input, loss,
bootstrap draw, implementation, independent comparison and standalone figure is
retained in its named directory. `publication.json` hashes every byte except
itself and distinguishes numerical agreement from predictive failure.

See the [study report](../../variant-transfer-v2-study.md) and
[original failed-v1 report](../../variant-transfer-study.md).

Verify with standard-library Python, without source download or refitting:

```bash
python3 -S research/results/variant-transfer-v2/packaging-implementation.py \\
  --verify research/results/variant-transfer-v2
```
''')
        (args.out / "publication.json").write_bytes(encoded({"schema_version": 1, "status": status,
            "original_runs": {name: str(path) for name, path in roots.items()}, "files": files(args.out)}))
        verify(args.out)
        return
    require(all((args.primary, args.independent, args.comparison, args.out)), "all source paths are required")
    require(not args.out.exists(), "preserve existing publications")
    status = failed_status(args.primary, args.independent, args.comparison)
    args.out.mkdir(parents=True)
    for name, source in (("primary-source", args.primary), ("independent-source", args.independent),
                         ("gate-comparison", args.comparison)):
        shutil.copytree(source, args.out / name)
    shutil.copyfile(__file__, args.out / "packaging-implementation.py")
    (args.out / "README.md").write_text('''# Retained v1 source-gate failure

Both original source readers stop on the same literal `NA` in cell AO42972,
K562 variant `chr12:72652619:A:G`. All 241 selected-cell records agree, including
240 ordinary numeric float bit patterns from 120 fully decoded paired variants.
The frozen policy accepted no undocumented missing token. **The source gate
failed; no model-accuracy comparison or bootstrap was computed.**

`primary-source/`, `independent-source/` and `gate-comparison/` preserve every
original artifact and manifest. The independent reader exported fixed predictions
before the source failure; those are retained and are not performance evidence.
Openpyxl may instantiate other worksheet values internally. The separately
versioned exploratory v2 does not change this v1 failure.

See the [study report](../../variant-transfer-study.md) for source/assembly,
control/exposure exclusions, the two committed freezes and all limitations.
The final metadata cohort had 2,172 variants in 172 original components.

Verify without new source access or dependencies:

```bash
python3 -S research/results/variant-transfer-v1/packaging-implementation.py \\
  --verify research/results/variant-transfer-v1
```

`publication.json` inventories all published bytes except itself. It binds both
original run manifests, their failed source comparison and the frozen authority.
''')
    (args.out / "publication.json").write_bytes(encoded({"schema_version": 1, "status": status,
        "original_runs": {"primary": str(args.primary), "independent": str(args.independent),
                          "comparison": str(args.comparison)}, "files": files(args.out)}))
    verify(args.out)


if __name__ == "__main__":
    main()
