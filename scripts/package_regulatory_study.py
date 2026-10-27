#!/usr/bin/env python3
"""Package complete study evidence, or verify a published snapshot (stdlib only)."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inventory(root):
    return {
        str(path.relative_to(root)): {"bytes": path.stat().st_size, "sha256": sha(path)}
        for path in sorted(root.rglob("*")) if path.is_file()
        and path != root / "publication.json"
    }


def verify_manifest(root, name="artifacts.json"):
    document = json.loads((root / name).read_bytes())
    for relative, expected in document["files"].items():
        path = (root / relative).resolve()
        if root.resolve() not in path.parents or not path.is_file() or sha(path) != expected:
            raise ValueError("artifact path/checksum mismatch: " + relative)
    return document


def deterministic_tar(source, destination):
    with destination.open("xb") as output:
        with gzip.GzipFile(fileobj=output, mode="wb", filename="", mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode="w|", format=tarfile.PAX_FORMAT) as archive:
                for path in sorted(source.rglob("*")):
                    if not path.is_file() or path.is_symlink():
                        raise ValueError("expected flat regular-file independent bundle")
                    info = tarfile.TarInfo("independent-python/" + str(path.relative_to(source)))
                    info.size = path.stat().st_size
                    info.mode = 0o644
                    info.mtime = info.uid = info.gid = 0
                    with path.open("rb") as handle:
                        archive.addfile(info, handle)


def verify(root):
    document = json.loads((root / "publication.json").read_bytes())
    actual = inventory(root)
    if document["files"] != actual:
        raise ValueError("published file inventory/checksum mismatch")
    verify_manifest(root / "core")
    verify_manifest(root / "rank-sensitivity")
    print(f"Verified {len(actual)} published files; integrity is not biological validation.")


def package(prefix, out, history):
    paths = {suffix: Path(str(prefix) + "-" + suffix)
             for suffix in ("rust", "python", "comparison", "figures", "rank-sensitivity")}
    if out.exists():
        raise ValueError("publication output exists; preserve the previous snapshot")
    for suffix in ("rust", "python", "rank-sensitivity"):
        verify_manifest(paths[suffix])
    figures = verify_manifest(paths["figures"], "figures.json")
    if figures["input_study_manifest_sha256"] != sha(paths["rust"] / "artifacts.json"):
        raise ValueError("figures do not identify these parent results")
    comparison = json.loads((paths["comparison"] / "comparison.json").read_bytes())
    if not comparison["primary_predictions_losses_and_controls_matched"]:
        raise ValueError("primary numerical comparison failed; review before publication")
    if comparison["failure_count"] != len(comparison["failures"]):
        raise ValueError("comparison failure inventory differs from its count")
    rank_checks = sum(comparison["checked"][key] for key in ("spearman_rows", "spearman_positions"))
    rank_failures = sum(item["category"] in ("spearman_rows", "spearman_positions")
                        for item in comparison["failures"])
    full_status = ("passed" if comparison["matched"] else "remains failed")
    replay_status = ("returns **0** for full numerical agreement" if comparison["matched"]
                     else "returns **2** for the retained comparison failures")
    rust_sha = sha(paths["rust"] / "study.json.gz")
    independent_sha = sha(paths["python"] / "artifacts.json")
    if (comparison["rust_study_sha256"] != rust_sha
            or comparison["independent_artifacts_sha256"] != independent_sha):
        raise ValueError("comparison does not identify these parent results")
    supplement = json.loads((paths["rank-sensitivity"] / "summary.json").read_bytes())
    if (supplement["provenance"]["rust_study_sha256"] != rust_sha
            or supplement["provenance"]["independent_artifacts_sha256"] != independent_sha):
        raise ValueError("supplement does not identify these parent results")
    out.mkdir(parents=True)
    for source, target in (("rust", "core"), ("comparison", "comparison"),
                           ("figures", "figures"), ("rank-sensitivity", "rank-sensitivity")):
        shutil.copytree(paths[source], out / target)
    deterministic_tar(paths["python"], out / "independent-python.tar.gz")
    if history is not None:
        shutil.copytree(history, out / "history")
    inputs = json.loads((paths["rust"] / "inputs.json").read_bytes())
    (out / "README.md").write_text(f"""# Regulatory study v1: published evidence

This is computational reanalysis of one published MPRA study, not independent
laboratory replication or a new biological mechanism. Read the
[research report](../../regulatory-study.md) for findings, methods and limits.

The core was generated from code commit `{inputs['software']['commit']}`;
working-tree dirty flag: `{inputs['software']['working_tree_dirty']}`.

**Primary predictions/losses/controls agree; full frozen numerical reproduction
{full_status}, with {rank_failures}/{rank_checks} secondary exact-rank failures.** Every failure is in
`comparison/comparison.json`. The separate post-hoc rank supplement reports all
five declared grids; it does not replace the failed endpoint.

- `core/`: complete compressed Rust predictions/fits, scoped proof, source and
  protocol commitments, losses, fold memberships, controls and revision ledger.
- `independent-python.tar.gz`: full independently fitted outputs, features,
  predictions, resampling, software record and exact implementation snapshot.
- `figures/`: five standalone SVG/PNG figures plus exact CSV data and provenance.
- `rank-sensitivity/`: all levels and all metric pairs, including raw failures.
- `history/` (when present): compressed earlier attempts and their original
  manifests/implementations where captured. Superseded artifacts are evidence
  of the investigation, not alternative successful endpoints.
- `publication.json`: SHA-256/size inventory of every published file, including
  histories. It records the partial-validation status explicitly.

From the repository root, verify integrity and replay the saved comparison:

```bash
source scripts/research-env.sh
python3 scripts/package_regulatory_study.py --verify research/results/regulatory-v1
cargo run --release --frozen -- audit research/results/regulatory-v1/core
mkdir -p runs/published-replay
tar -xzf research/results/regulatory-v1/independent-python.tar.gz -C runs/published-replay
python3 scripts/replicate_regulatory_study.py \\
  --compare-rust research/results/regulatory-v1/core/study.json.gz \\
  --independent-run runs/published-replay/independent-python \\
  --out runs/published-replay/comparison
```

The last command {replay_status}. Use a new
output directory for another replay. `bash scripts/check_study.sh runs/<new-prefix>`
refits both implementations and regenerates diagnostics and figures.

Source terms remain attached to the source bundles: MaveDB score sets declare
CC0; JASPAR priors are CC BY4.0. Retained author code/documentation are GPL-3.0
reference evidence, and the original OSF table has separately documented producer
terms. Those sources are not relabeled under Helix's MIT code license. Hashes
establish integrity relative to commitments, not authentication or empirical truth.
""")
    publication = {
        "schema_version": 1,
        "software": inputs["software"],
        "protocol_sha256": inputs["protocol_sha256"],
        "primary_predictions_losses_and_controls_matched": True,
        "full_frozen_numerical_reproduction_matched": comparison["matched"],
        "secondary_exact_rank_metrics_matched": comparison["secondary_exact_rank_metrics_matched"],
        "secondary_failure_count": rank_failures,
        "secondary_exact_rank_checks": rank_checks,
        "posthoc_supplement_replaces_parent_endpoint": False,
        "files": inventory(out),
    }
    (out / "publication.json").write_text(json.dumps(publication, indent=2) + "\n")
    verify(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--history", type=Path)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    if args.verify is not None:
        verify(args.verify)
    elif args.prefix is not None and args.out is not None:
        package(args.prefix, args.out, args.history)
    else:
        parser.error("provide --verify, or --prefix and --out")


if __name__ == "__main__":
    main()
