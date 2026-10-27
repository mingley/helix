#!/usr/bin/env python3
"""Publish integrity-bound context evidence without promoting unresolved semantics."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
import shlex
import shutil
import tarfile


ENDPOINT = "published_element_specificity_score"
CORE_FILES = {
    "features.json.gz", "findings.md", "grouping.json", "inputs.json",
    "ledger.json", "protocol.json", "reference.json", "resampling.json.gz",
    "source-manifest.json", "study.json.gz", "summary.json",
}
SEMANTIC_SCOPE = (
    "Historical semantic_endpoint_acceptance:true receipts identify the published "
    "columns only; their formula acceptance is superseded by the unresolved "
    "source-formula audit. They are preserved as historical evidence."
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def relative_name(name):
    require(isinstance(name, str) and name and "\\" not in name,
            "invalid relative artifact name")
    path = PurePosixPath(name)
    require(not path.is_absolute() and all(p not in (".", "..") for p in name.split("/"))
            and str(path) == name, "unsafe artifact path: " + name)
    return name


def digest(handle):
    result = hashlib.sha256()
    size = 0
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        size += len(block)
        result.update(block)
    return {"bytes": size, "sha256": result.hexdigest()}


class Files:
    """Verified regular files, from a directory or a non-extracted Python tar."""
    def __init__(self, root, archive=False):
        self.root, self.archive = Path(root), archive
        self.files = {}
        if archive:
            with tarfile.open(self.root, "r:gz") as bundle:
                for member in bundle:
                    require(member.isfile() and member.name.startswith("independent-python/"),
                            "Python archive must contain regular files under independent-python/")
                    name = relative_name(member.name[len("independent-python/"):])
                    require(name not in self.files, "duplicate Python archive member: " + name)
                    with bundle.extractfile(member) as handle:
                        self.files[name] = digest(handle)
        else:
            require(self.root.is_dir() and not self.root.is_symlink(),
                    "expected regular artifact directory: " + str(root))
            for path in sorted(self.root.rglob("*")):
                require(not path.is_symlink(), "symlink in evidence: " + str(path))
                if path.is_dir():
                    continue
                require(path.is_file(), "nonregular evidence: " + str(path))
                with path.open("rb") as handle:
                    self.files[relative_name(path.relative_to(self.root).as_posix())] = digest(handle)

    def read(self, name):
        relative_name(name)
        require(name in self.files, "missing artifact: " + name)
        if self.archive:
            with tarfile.open(self.root, "r:gz") as bundle:
                return bundle.extractfile("independent-python/" + name).read()
        return (self.root / name).read_bytes()

    def json(self, name):
        return json.loads(self.read(name))

    def sha(self, name):
        require(name in self.files, "missing artifact: " + name)
        return self.files[name]["sha256"]


def manifest(files, name, sized=False):
    document = files.json(name)
    expected = document if sized else document["files"]
    require(isinstance(expected, dict) and expected, "empty artifact manifest: " + name)
    require(set(files.files) == set(expected) | {name},
            "manifest inventory differs from directory/archive: " + str(files.root))
    for relative, commitment in expected.items():
        relative_name(relative)
        actual = files.files.get(relative)
        require(actual == commitment if sized else actual is not None and actual["sha256"] == commitment,
                "artifact size/checksum mismatch: " + str(files.root) + "/" + relative)
    return document


def raw_binding(document, key, files, name):
    commitment = document[key]
    require(commitment.get("raw_sha256") == files.sha(name)
            and commitment.get("raw_bytes") == files.files[name]["bytes"],
            "comparison raw artifact binding differs: " + key)


def validate(core, independent, comparison_files, figures_files, arithmetic_files, endpoint_files=None):
    artifacts = manifest(core, "artifacts.json")
    require(artifacts.get("schema_version") == 4 and CORE_FILES <= set(artifacts["files"]),
            "context core must contain the complete schema-4 artifacts")
    inputs, summary = core.json("inputs.json"), core.json("summary.json")
    source = core.json("source-manifest.json")
    protocol_sha, panel_sha = inputs["protocol_sha256"], inputs["panel_sha256"]
    require(inputs.get("endpoint") == summary.get("endpoint") == source.get("endpoint") == ENDPOINT,
            "core endpoint differs from published ESS")
    require(core.sha("protocol.json") == protocol_sha
            and core.sha("source-manifest.json") == inputs["source_manifest_sha256"]
            and core.sha("grouping.json") == inputs["grouping_sha256"]
            and core.sha("reference.json") == source["files"].get("reference-validation.json")
            and source.get("protocol_sha256") == protocol_sha
            and source["files"].get("panel.json") == panel_sha,
            "core protocol/source/panel/grouping bindings differ")
    require(core.json("reference.json").get("reference_verified") is True
            and source.get("record_count") == summary["sequence_count"]
            and source.get("contexts") == summary["contexts"],
            "core reference/count/context commitments differ")
    semantics = summary["source_semantics"]
    require(semantics.get("full_scientific_semantic_validation") is False
            and str(semantics.get("status", "")).startswith("unresolved")
            and semantics.get("raw_source_values_preserved") is True,
            "this publication requires preserved targets and explicit unresolved scientific semantics")

    independent_manifest = manifest(independent, "manifest.json")
    own_inputs = independent.json("inputs.json")
    own_status = independent.json("status.json")
    require(independent_manifest.get("protocol_sha256") == protocol_sha
            and independent_manifest.get("endpoint") == ENDPOINT
            and own_inputs.get("protocol_sha256") == protocol_sha
            and own_inputs.get("panel_sha256") == panel_sha
            and own_inputs.get("manifest_sha256") == core.sha("source-manifest.json")
            and own_inputs.get("endpoint") == ENDPOINT
            and own_status.get("status") == "complete"
            and own_status.get("protocol_sha256") == protocol_sha
            and own_status.get("endpoint") == ENDPOINT
            and own_status.get("source_elements") == summary["sequence_count"]
            and own_status.get("source_measurements") == summary["measurement_count"],
            "independent run is incomplete or binds different inputs")
    require(own_inputs.get("implementation_sha256") == independent.sha("implementation.py"),
            "independent implementation commitment differs")

    comparison_manifest = manifest(comparison_files, "manifest.json")
    comparison = comparison_files.json("comparison.json")
    require(comparison.get("matched") is True and comparison.get("failure_count") == 0
            and comparison.get("failures") == [] and comparison.get("checks")
            and all(check.get("matched") is True for check in comparison["checks"]),
            "complete numerical reproduction has not matched")
    require(comparison_manifest.get("protocol_sha256") == comparison.get("protocol_sha256") == protocol_sha
            and comparison.get("endpoint") == ENDPOINT
            and comparison["rust_artifacts_manifest"].get("sha256") == core.sha("artifacts.json")
            and comparison["independent_manifest"].get("sha256") == independent.sha("manifest.json")
            and comparison.get("implementation_sha256") == comparison_files.sha("implementation.py"),
            "comparison manifest/implementation/parent bindings differ")
    for key, files, name in (
            ("rust_study", core, "study.json.gz"),
            ("rust_features", core, "features.json.gz"),
            ("rust_resampling", core, "resampling.json.gz"),
            ("independent_study", independent, "evaluation.json.gz"),
            ("independent_resampling", independent, "resampling.json.gz"),
            ("independent_raw_source_records", independent, "source-records.json.gz")):
        raw_binding(comparison, key, files, name)

    figures = manifest(figures_files, "figures.json")
    require(figures.get("input_study_manifest_sha256") == core.sha("artifacts.json"),
            "figures bind another Rust study manifest")
    manifest(arithmetic_files, "artifacts.json", sized=True)
    arithmetic = arithmetic_files.json("summary.json")
    if figures.get("arithmetic_artifacts_sha256") is not None:
        require(figures["arithmetic_artifacts_sha256"] == arithmetic_files.sha("artifacts.json"),
                "figures bind another arithmetic audit")
    hashes = arithmetic["hashes"]
    require(hashes.get("parent_protocol_sha256") == protocol_sha
            and hashes.get("panel_sha256") == panel_sha
            and hashes.get("audit_protocol_sha256") == arithmetic_files.sha("protocol.json")
            and hashes.get("implementation_sha256") == arithmetic_files.sha("implementation.py")
            and hashes.get("source_workbook_sha256") == source["files"].get("raw/Table-S9-folds-performance.xlsx")
            and arithmetic.get("historical_normalization_verified") is False,
            "arithmetic audit has different inputs or claims verified historical normalization")

    if endpoint_files is not None:
        endpoint_manifest = manifest(endpoint_files, "manifest.json")
        endpoint = endpoint_files.json("endpoint-audit.json")
        commitments = endpoint["input_commitments"]
        require(endpoint_manifest.get("protocol_sha256") == endpoint.get("protocol_sha256") == protocol_sha
                and endpoint.get("endpoint") == ENDPOINT
                and commitments.get("panel_sha256") == panel_sha
                and commitments.get("source_manifest_sha256") == core.sha("source-manifest.json")
                and commitments["independent_manifest"].get("sha256") == independent.sha("manifest.json"),
                "endpoint audit binds different parent inputs")
        raw_binding(commitments, "source_records", independent, "source-records.json.gz")

    return {
        "schema_version": 1,
        "endpoint": ENDPOINT,
        "software": inputs["software"],
        "protocol_sha256": protocol_sha,
        "panel_sha256": panel_sha,
        "source_manifest_sha256": core.sha("source-manifest.json"),
        "core_artifacts_manifest_sha256": core.sha("artifacts.json"),
        "independent_manifest_sha256": independent.sha("manifest.json"),
        "complete_numerical_reproduction": comparison["matched"],
        "full_scientific_semantic_validation": semantics["full_scientific_semantic_validation"],
        "scientific_semantic_status": semantics["status"],
        "historical_normalization_verified": arithmetic["historical_normalization_verified"],
        "overall_frozen_protocol_success": comparison["matched"] and semantics["full_scientific_semantic_validation"],
        "historical_semantic_acceptance_scope": SEMANTIC_SCOPE,
        "counts": {key: summary[key] for key in
                   ("sequence_count", "measurement_count", "group_count", "contexts")},
    }


def deterministic_tar(source, destination):
    with destination.open("xb") as output:
        with gzip.GzipFile(fileobj=output, mode="wb", filename="", mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode="w|", format=tarfile.PAX_FORMAT) as bundle:
                for name, commitment in sorted(source.files.items()):
                    info = tarfile.TarInfo("independent-python/" + name)
                    info.size = commitment["bytes"]
                    info.mode = 0o644
                    info.mtime = info.uid = info.gid = 0
                    with (source.root / name).open("rb") as handle:
                        bundle.addfile(info, handle)


def verify(root):
    published = Files(root)
    publication = published.json("publication.json")
    expected = {name: entry for name, entry in published.files.items() if name != "publication.json"}
    require(publication.get("files") == expected, "publication file inventory/size/checksum mismatch")
    status = validate(Files(root / "core"), Files(root / "independent-python.tar.gz", archive=True),
                      Files(root / "comparison"), Files(root / "figures"), Files(root / "arithmetic"),
                      Files(root / "endpoint-audit") if (root / "endpoint-audit").exists() else None)
    require({k: v for k, v in publication.items() if k != "files"} == status,
            "publication statuses differ from bound evidence")
    print(json.dumps({"publication_verified": True, "published_files": len(expected),
                      "complete_numerical_reproduction": status["complete_numerical_reproduction"],
                      "full_scientific_semantic_validation": status["full_scientific_semantic_validation"],
                      "scientific_semantic_status": status["scientific_semantic_status"]}))


def readme(out, status):
    location = shlex.quote(out.as_posix())
    software, counts = status["software"], status["counts"]
    return f"""# Context-specificity evidence snapshot

Complete numerical reproduction: **{str(status['complete_numerical_reproduction']).lower()}**.
Full scientific semantic validation: **{str(status['full_scientific_semantic_validation']).lower()}**;
status: `{status['scientific_semantic_status']}`. Overall frozen-protocol success
is false. Rust/Python agreement does not resolve the source ESS formula mismatch.
Historical normalization verified: **{str(status['historical_normalization_verified']).lower()}**.

This is computational reanalysis of published Element Specificity Score columns,
not raw activity, independent laboratory replication, or a new biological mechanism.
The source's whole-panel normalization includes held-out elements; it is not
prospective or training-only target preprocessing. Published values remain unchanged.
{SEMANTIC_SCOPE}

The bound core reports {counts['sequence_count']:,} sequences,
{counts['measurement_count']:,} measurements, {counts['group_count']:,} groups,
and contexts {', '.join(counts['contexts'])}. Source code commit:
`{software['commit']}`; working-tree dirty flag: `{software['working_tree_dirty']}`.

- `core/`: schema-4 Rust fits, predictions, features, resampling, ledger, frozen
  protocol, source commitments, reference validation and explicit semantic limits.
- `independent-python.tar.gz`: deterministic archive of the complete independent
  implementation snapshot, source checks, independently fitted outputs and manifests.
- `comparison/`: complete numerical comparison and its implementation snapshot.
- `figures/`: standalone figures, plotted data and exact parent-result commitment.
- `arithmetic/`: post-hoc source arithmetic diagnostic. A fitted weighted plane
  does not identify historical normalization and never changes benchmark targets.
- `endpoint-audit/` and `history/`, when present: independent endpoint audit and
  retained earlier attempts. Superseded results are evidence, not alternate successes.
- `publication.json`: byte sizes and SHA-256 for every other published file,
  plus machine-readable numerical and scientific statuses.

From the repository root, audit the core and verify the complete package:

```bash
source scripts/research-env.sh
cargo run --release --frozen -- audit {location}/core
python3 scripts/package_context_study.py --verify {location}
```

Replay only the saved numerical comparison, using a new replay directory:

```bash
mkdir -p runs/context-published-replay
tar -xzf {location}/independent-python.tar.gz -C runs/context-published-replay
python3 scripts/replicate_context_study.py \\
  --compare-rust {location}/core/study.json.gz \\
  --independent-run runs/context-published-replay/independent-python \\
  --out runs/context-published-replay/comparison
```

Numerical replay returns 0 when agreement holds. To refit both implementations,
regenerate diagnostics and retain the distinct semantic failure, run
`bash scripts/check_context_study.sh runs/<new-prefix>` (replace the placeholder).
The full workflow is expected to return **2** while source semantic validation
remains unresolved, even when its complete numerical comparison matches.

The raw source dataset and its licenses/receipts remain in repository
`datasets/context-panel/`; this bundle commits to them without duplicating them.
The roughly 800 MB hg38 reference cache is separate, under the configured tools
directory; `bash scripts/setup_reference.sh` acquires and verifies it for full replay.
The snapshot is not a standalone replacement for that repository and cache.

Attribute source workbooks to Agarwal et al., *Massively parallel characterization
of transcriptional regulatory elements*, DOI
[10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9), and the
preprint-associated deposit [Zenodo 10558183](https://zenodo.org/records/10558183).
The deposit declares [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/);
retain attribution and the source-version distinction recorded in the manifests.
Source material is not relicensed under the repository's code license. Hashes
establish integrity relative to commitments, not empirical truth or source authorship.
"""


def package(prefix, out, arithmetic, endpoint_audit=None, history=None):
    require(not out.exists() and not out.is_symlink(), "output exists; preserve previous snapshots")
    paths = {key: Path(str(prefix) + "-" + key) for key in ("rust", "python", "comparison", "figures")}
    sources = {key: Files(path) for key, path in paths.items()}
    arithmetic_files = Files(arithmetic)
    endpoint_files = Files(endpoint_audit) if endpoint_audit is not None else None
    if history is not None:
        Files(history)  # Reject symlinks/nonregular entries before publication.
    for source in list(paths.values()) + [arithmetic, endpoint_audit, history]:
        if source is not None:
            require(out.resolve() != source.resolve() and source.resolve() not in out.resolve().parents,
                    "publication output cannot be inside an input directory")
    status = validate(sources["rust"], sources["python"], sources["comparison"],
                      sources["figures"], arithmetic_files, endpoint_files)
    out.mkdir(parents=True)
    for key, target in (("rust", "core"), ("comparison", "comparison"), ("figures", "figures")):
        shutil.copytree(paths[key], out / target)
    shutil.copytree(arithmetic, out / "arithmetic")
    deterministic_tar(sources["python"], out / "independent-python.tar.gz")
    if endpoint_audit is not None:
        shutil.copytree(endpoint_audit, out / "endpoint-audit")
    if history is not None:
        shutil.copytree(history, out / "history")
    (out / "README.md").write_text(readme(out, status), encoding="utf-8")
    status["files"] = Files(out).files
    (out / "publication.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    verify(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("prefix", "out", "arithmetic", "endpoint-audit", "history", "verify"):
        parser.add_argument("--" + name, type=Path)
    args = parser.parse_args()
    if args.verify is not None:
        if any(value is not None for key, value in vars(args).items() if key != "verify"):
            parser.error("--verify cannot be combined with packaging arguments")
        verify(args.verify)
    elif all(value is not None for value in (args.prefix, args.out, args.arithmetic)):
        package(args.prefix, args.out, args.arithmetic, args.endpoint_audit, args.history)
    else:
        parser.error("provide --verify, or --prefix, --out and --arithmetic")


if __name__ == "__main__":
    main()
