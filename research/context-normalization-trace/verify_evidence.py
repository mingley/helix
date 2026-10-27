#!/usr/bin/env python3
"""Offline byte/blob and code-AST checks; never executes producer code or reads assay data."""
import ast
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def inventory():
    result = {}
    for path in sorted(ROOT.rglob("*")):
        if path.is_symlink():
            raise ValueError("symlink in evidence: " + str(path))
        if not path.is_file() or path.name == "manifest.json" or "__pycache__" in path.parts:
            continue
        raw = path.read_bytes()
        result[path.relative_to(ROOT).as_posix()] = {"bytes": len(raw), "sha256": sha(raw)}
    return result


def function_ast(path):
    module = ast.parse(path.read_text())
    function = next(node for node in module.body if isinstance(node, ast.FunctionDef)
                    and node.name == "validate_joint_ess")
    return sha(ast.dump(function, include_attributes=False).encode())


def main():
    manifest = json.loads((ROOT / "manifest.json").read_bytes())
    if inventory() != manifest["files"]:
        raise ValueError("evidence inventory/size/SHA-256 differs")
    trace = json.loads((ROOT / "trace.json").read_bytes())
    for entry in trace["git_blob_bindings"]:
        path = ROOT / entry["file"]
        if ROOT not in path.resolve().parents:
            raise ValueError("unsafe evidence path")
        raw = path.read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        if blob != entry["git_blob_sha1"]:
            raise ValueError("Git blob differs: " + entry["file"])
    expected = trace["downstream_hypothesis"]["function_ast_sha256"]
    for relative in ("dalal/LentiMPRA/process.py", "dalal-initial/LentiMPRA/process.py"):
        if function_ast(ROOT / relative) != expected:
            raise ValueError("initial/head code hypothesis AST differs")
    if ((ROOT / "mpraflow-current/src/make_master_tables.R").read_bytes()
            != (ROOT / "mpraflow-2020/src/make_master_tables.R").read_bytes()):
        raise ValueError("checked MPRAflow master-table snapshots differ")
    if any(Path(relative).suffix in (".xlsx", ".tsv", ".parquet", ".zip", ".gz")
           for relative in manifest["files"]):
        raise ValueError("unexpected quantitative/archive payload in code-only evidence")
    print(json.dumps({"files_verified": len(manifest["files"]),
                      "git_blobs_verified": len(trace["git_blob_bindings"]),
                      "initial_head_function_ast_identical": True,
                      "assay_payloads_in_this_snapshot": False,
                      "source_numeric_identity_verified_by_this_script": False,
                      "original_producer_history_certified": False}))


if __name__ == "__main__":
    main()
