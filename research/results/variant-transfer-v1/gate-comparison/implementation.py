#!/usr/bin/env python3
"""Compare independently retained source-gate failures, not model performance."""
import argparse
import ast
import gzip
import hashlib
import json
from pathlib import Path
import shutil

BINDING = "9c26b3fc628f16bd6bb677dee8eed01e509b6102d44f95c0037ae2f55676ba79"
PROTOCOL = "811930c2d3a5ee763014ea611acda4d2c84404b7d5ae399c7843c62d0a828149"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return json.load(handle)


def verify(root):
    manifest = load(root / "manifest.json")
    actual = {p.relative_to(root).as_posix(): p for p in root.rglob("*")
              if p.is_file() and p.relative_to(root).as_posix() != "manifest.json"}
    assert set(actual) == set(manifest["files"])
    for name, item in manifest["files"].items():
        path = actual[name]
        assert not path.is_symlink() and path.resolve().is_relative_to(root.resolve())
        assert sha(path) == (item["sha256"] if isinstance(item, dict) else item)
        if isinstance(item, dict):
            assert path.stat().st_size == item["bytes"]
    return {"sha256": sha(root / "manifest.json"), "bytes": (root / "manifest.json").stat().st_size}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--independent", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "Preserve previous attempts"
    pins = {"primary_manifest": verify(args.primary), "independent_manifest": verify(args.independent)}
    own_binding, primary_binding = load(args.independent / "binding.json"), load(args.primary / "binding.json")
    assert sha(args.independent / "binding.json") == sha(args.primary / "binding.json") == BINDING
    assert own_binding == primary_binding and own_binding["missing_policy"]["allowed_tokens"] == []
    assert sha(args.independent / "protocol.json") == sha(args.primary / "protocol.json") == PROTOCOL
    with gzip.open(args.primary / "source-type-float-ledger.jsonl.gz", "rt") as handle:
        primary = [json.loads(line) for line in handle]
    independent = load(args.independent / "source-cell-type-ledger.json.gz")
    checks, failures = [], []
    for number, (a, b) in enumerate(zip(primary, independent)):
        matched = (a["id"] == b["id"] and a["context"] == b["context"]
                   and a["row_1based"] == b["source_row_1based"] and a["xml_type"] == b["cell_data_type"])
        if a["status"] == "ordinary_finite_numeric":
            matched = matched and b["status"] == "ordinary_numeric" and a["float64_le_hex"] == b["float64_le_hex"]
        elif a["status"] == "gate_failure":
            matched = matched and b["status"] == "gate_failure" and a["text"] == ast.literal_eval(
                b["openpyxl_value_repr_not_raw_XML_lexeme"])
        else:
            matched = False
        record = {"selected_cell_ordinal_0based": number, "id": a["id"], "context": a["context"],
            "source_row_1based": a["row_1based"], "status": a["status"], "matched": matched}
        checks.append(record)
        if not matched:
            failures.append(record)
    assert len(primary) == len(independent) == len(checks)
    assert len(primary) > 0 and primary[-1]["status"] == independent[-1]["status"] == "gate_failure"
    assert all(a["status"] == "ordinary_finite_numeric" for a in primary[:-1])
    assert not (args.primary / "evaluation.json.gz").exists()
    assert not (args.independent / "evaluation.json.gz").exists()
    failure = primary[-1]
    paired = len(primary[:-1]) // 2
    assert len(primary[:-1]) == 2 * paired
    for i in range(paired):
        assert primary[2 * i]["id"] == primary[2 * i + 1]["id"]
        assert [primary[2 * i]["context"], primary[2 * i + 1]["context"]] == ["K562", "HepG2"]
    summary = {"schema_version": 1, "kind": "independently_reproduced_variant_transfer_v1_source_gate_failure",
        "protocol_sha256": PROTOCOL, "binding_sha256": BINDING, **pins,
        "source_gate_failure_reproduced": not failures, "source_observation_gate_passed": False,
        "selected_cell_comparisons": len(checks), "ordinary_numeric_float_bits_compared_exactly": len(primary) - 1,
        "fully_decoded_paired_variants_before_failure": paired,
        "failure_count": len(failures), "failures": failures,
        "first_failure": {"id": failure["id"], "context": failure["context"],
            "coordinate": failure["coordinate"], "row_1based": failure["row_1based"],
            "xml_type": failure["xml_type"], "literal_text": failure["text"]},
        "allowed_missing_tokens": [], "predictive_performance_evaluated": False,
        "bootstrap_computed": False, "performance_result": None,
        "full_protocol_independent_reproduction_claimed": False,
        "limits": ["This reproduces a failed source gate; it does not estimate transfer performance.",
            "The first unsupported string is retained without reinterpreting it as source-documented missingness.",
            "Openpyxl may instantiate nonselected worksheet cell values internally; only selected-cell use is audited here.",
            "No numeric agreement can make the frozen source gate pass."]}
    args.out.mkdir(parents=True)
    (args.out / "comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
    (args.out / "selected-cell-checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    shutil.copyfile(Path(__file__), args.out / "implementation.py")
    manifest = {p.name: {"sha256": sha(p), "bytes": p.stat().st_size} for p in sorted(args.out.iterdir())}
    (args.out / "manifest.json").write_text(json.dumps({"schema_version": 1, "files": manifest}, indent=2) + "\n")
    verify(args.out)
    print(json.dumps(summary, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
