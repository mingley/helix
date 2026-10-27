#!/usr/bin/env python3
"""Retain and verify exact variant-transfer evidence with the standard library."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil

V1 = "811930c2d3a5ee763014ea611acda4d2c84404b7d5ae399c7843c62d0a828149"
V1_BINDING = "9c26b3fc628f16bd6bb677dee8eed01e509b6102d44f95c0037ae2f55676ba79"


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


def verify(root):
    checked(root, "publication.json")
    status = failed_status(root / "primary-source", root / "independent-source", root / "gate-comparison")
    require(load(root / "publication.json")["status"] == status, "publication status was changed")
    print(json.dumps({"publication_verified": True, **status}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--primary", type=Path)
    parser.add_argument("--independent", type=Path)
    parser.add_argument("--comparison", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.verify:
        verify(args.verify)
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
