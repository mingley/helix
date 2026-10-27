from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from importlib import metadata
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "manifest.json"
MANIFEST_VERSION = 1
ENVIRONMENT_NAME = "lentimpra-enhancer-processing"
ENVIRONMENT_PLATFORM = "linux-64"
ENVIRONMENT_LOCK = "environment-lock-linux-64.txt"
SOURCE_OUTPUT_SEMANTICS = "source_output_semantics"
SEMANTIC_FINGERPRINT_ALGORITHM = "arrow-ipc-canonical-column-sha256-v1"
PROCESSING_REPOSITORY = (
    "https://github.com/TaykhoomDalal/LentiMPRA-Enhancer-Processing"
)
PAPER_TITLE = (
    "Massively parallel characterization of transcriptional regulatory elements"
)
PAPER_DOI = "10.1038/s41586-024-08430-9"
PAPER_DOI_URL = f"https://doi.org/{PAPER_DOI}"
ARCHIVE_NAME = "41586_2024_8430_MOESM5_ESM.zip"
ARCHIVE_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41586-024-08430-9/"
    "MediaObjects/41586_2024_8430_MOESM5_ESM.zip"
)
ARCHIVE_SHA256 = "aefb1974e9ebe5c196a89a00f4e868c0abaf69e3d3952e79c74bfb42f217f6a4"
ARCHIVE_BYTES = 342_234_089
ARCHIVE_PREFIX = "2023-03-03628-s5/"
WORKBOOKS = {
    "S3": {
        "filename": "2023-03-03628C-Table S3 - large-scale lib design (1).xlsx",
        "bytes": 56_198_290,
        "sha256": "cc5164f62629a218d2d99485f220e9f4f29ccd4c8f7879003760b8163975fde3",
    },
    "S4": {
        "filename": "2023-03-03628C-Table S4 - largescaledata (1).xlsx",
        "bytes": 133_841_008,
        "sha256": "5c24214988a7551186e62426d4d86040dc406988a4580b99d677f4b2b7120c7c",
    },
    "S8": {
        "filename": "2023-03-03628C-Table S8 - folds and performance.xlsx",
        "bytes": 6_126_664,
        "sha256": "5788f1d7343b80746525035bcc6901241dc0f9da7d6e4354f1f36a6e5c5d0575",
    },
    "S10": {
        "filename": "2023-03-03628C-Table S10 - joint lib design.xlsx",
        "bytes": 8_301_566,
        "sha256": "f9ec4c6c2d230b0ef507a39e5d5c73b3b36f90ccc5758c33b442c7561e20e7e5",
    },
    "S11": {
        "filename": "2023-03-03628C-Table S11 - jointdata.xlsx",
        "bytes": 53_557_805,
        "sha256": "9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e",
    },
    "S12": {
        "filename": "2023-03-03628C-Table S12 - folds and performance.xlsx",
        "bytes": 3_564_742,
        "sha256": "0f9ff22d4b46fab36d419399e6ab743dc08a0813f812dc7bd96f90eb7e0f2b8a",
    },
}
OUTPUT_PATHS = (
    "LentiMPRA/lentimpra-hepg2.parquet",
    "LentiMPRA/lentimpra-k562.parquet",
    "LentiMPRA/lentimpra-wtc11.parquet",
    "LentiMPRA/lentimpra-joint.parquet",
    "LentiMPRA/published-performance-large-scale.parquet",
    "LentiMPRA/published-performance-joint.parquet",
)
HF_FILES = {
    "README.md": "LentiMPRA/HF_DATASET_CARD.md",
    **{Path(path).name: path for path in OUTPUT_PATHS},
}
PROCESSING_SOURCE_PATHS = (
    "LentiMPRA/process.py",
    "scripts/build_hf_release.py",
    "scripts/common.py",
    "scripts/create_environment.sh",
    "scripts/download_data.py",
    "scripts/refresh_manifest.py",
    "scripts/release_staged_outputs.py",
    "scripts/run_all.sh",
    "scripts/verify_outputs.py",
    "environment.yml",
    ENVIRONMENT_LOCK,
)
DIRECT_PACKAGES = {
    "numpy": "2.4.6",
    "pandas": "3.0.3",
    "pyarrow": "24.0.0",
    "python": "3.13.13",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, allow_nan=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def require_exact_keys(value: dict[str, Any], expected: set[str], where: str) -> None:
    actual = set(value)
    require(
        actual == expected,
        f"{where}: expected keys {sorted(expected)}, got {sorted(actual)}",
    )


def source_specification() -> dict[str, Any]:
    return {
        "paper": {
            "title": PAPER_TITLE,
            "doi": PAPER_DOI,
            "doi_url": PAPER_DOI_URL,
        },
        "archive": {
            "data_path": ARCHIVE_NAME,
            "url": ARCHIVE_URL,
            "bytes": ARCHIVE_BYTES,
            "sha256": ARCHIVE_SHA256,
        },
        "workbooks": {
            role: {
                **record,
                "data_path": f"raw/{record['filename']}",
                "archive_member": ARCHIVE_PREFIX + record["filename"],
            }
            for role, record in WORKBOOKS.items()
        },
    }


def explicit_artifacts(path: Path) -> dict[str, str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    require("@EXPLICIT" in lines, f"{path}: not a valid explicit Conda lock")
    artifacts = {}
    for line in lines:
        if not line or line.startswith(("#", "@")):
            continue
        require("#" in line, f"{path}: explicit artifact lacks SHA-256: {line}")
        url, checksum = line.rsplit("#", 1)
        require(
            re.fullmatch(r"[0-9a-f]{64}", checksum) is not None,
            f"{path}: invalid artifact SHA-256: {line}",
        )
        require(url not in artifacts, f"{path}: duplicate explicit artifact: {url}")
        artifacts[url] = checksum
    require(artifacts, f"{path}: explicit lock has no artifacts")
    return artifacts


def explicit_urls(path: Path) -> set[str]:
    return set(explicit_artifacts(path))


def verify_locked_environment(root: Path = ROOT) -> None:
    lock = root / ENVIRONMENT_LOCK
    require(lock.is_file(), f"missing explicit environment lock: {lock}")
    text = lock.read_text(encoding="utf-8")
    require(
        f"# platform: {ENVIRONMENT_PLATFORM}" in text,
        f"{lock}: expected platform {ENVIRONMENT_PLATFORM}",
    )
    require(
        sys.platform == "linux" and platform.machine().lower() in {"x86_64", "amd64"},
        f"expected {ENVIRONMENT_PLATFORM}, got {sys.platform}-{platform.machine()}",
    )
    require(
        Path(sys.prefix).name == ENVIRONMENT_NAME,
        f"expected environment {ENVIRONMENT_NAME}, got {sys.prefix}",
    )
    observed_versions = {
        "python": platform.python_version(),
        "numpy": metadata.version("numpy"),
        "pandas": metadata.version("pandas"),
        "pyarrow": metadata.version("pyarrow"),
    }
    require(
        observed_versions == DIRECT_PACKAGES,
        f"direct package versions changed: {observed_versions}",
    )
    completed = subprocess.run(
        ["conda", "list", "--prefix", sys.prefix, "--explicit", "--sha256"],
        check=True,
        capture_output=True,
        text=True,
    )
    installed = explicit_artifacts_from_lines(
        completed.stdout.splitlines(),
        "installed Conda explicit metadata",
    )
    locked = explicit_artifacts(lock)
    require(
        installed == locked,
        "installed Conda artifact URL/SHA-256 set does not exactly match "
        f"{ENVIRONMENT_LOCK}: missing={sorted(set(locked) - set(installed))}, "
        f"extra={sorted(set(installed) - set(locked))}, "
        f"checksum_mismatch={sorted(url for url in set(locked) & set(installed) if locked[url] != installed[url])}",
    )
    metadata_artifacts = {}
    metadata_root = Path(sys.prefix) / "conda-meta"
    for record_path in sorted(metadata_root.glob("*.json")):
        record = json.loads(record_path.read_text(encoding="utf-8"))
        url = record.get("url")
        checksum = record.get("sha256")
        require(
            isinstance(url, str)
            and re.fullmatch(r"[0-9a-f]{64}", checksum or "") is not None,
            f"{record_path}: missing installed artifact URL/SHA-256 metadata",
        )
        require(
            url not in metadata_artifacts,
            f"{metadata_root}: duplicate installed artifact URL: {url}",
        )
        metadata_artifacts[url] = checksum
    require(
        metadata_artifacts == locked,
        "installed conda-meta artifact identities do not exactly match "
        f"{ENVIRONMENT_LOCK}",
    )


def explicit_artifacts_from_lines(
    lines: list[str],
    label: str,
) -> dict[str, str]:
    artifacts = {}
    for line in lines:
        if not line or line.startswith(("#", "@")):
            continue
        require("#" in line, f"{label}: artifact lacks SHA-256: {line}")
        url, checksum = line.rsplit("#", 1)
        require(
            re.fullmatch(r"[0-9a-f]{64}", checksum) is not None,
            f"{label}: invalid artifact SHA-256: {line}",
        )
        require(url not in artifacts, f"{label}: duplicate artifact: {url}")
        artifacts[url] = checksum
    require(artifacts, f"{label}: no artifacts")
    return artifacts


def processing_record(root: Path = ROOT) -> dict[str, Any]:
    source_files = {
        path: sha256(root / path)
        for path in PROCESSING_SOURCE_PATHS
    }
    serialized = json.dumps(
        source_files,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return {
        "repository_url": PROCESSING_REPOSITORY,
        "revision": hashlib.sha256(serialized).hexdigest(),
        "source_files": source_files,
    }


def describe_parquet(path: Path) -> dict[str, Any]:
    parquet = pq.ParquetFile(path)
    return {
        "bytes": path.stat().st_size,
        "rows": parquet.metadata.num_rows,
        "sha256": sha256(path),
        "schema": [
            {
                "name": field.name,
                "type": str(field.type),
                "nullable": field.nullable,
            }
            for field in parquet.schema_arrow
        ],
    }


def task_contract() -> dict[str, Any]:
    return {
        "task": "lentimpra-enhancer-activity",
        "directory": "LentiMPRA",
        "input_columns": ["sequence"],
        "sequence": {
            "column": "sequence",
            "length": 200,
            "definition": (
                "Exact assayed insert after removing the two constant 15 nt "
                "synthesis adapters; use unchanged as model input."
            ),
            "orientation": (
                "Reverse constructs already contain the assayed reverse-complemented "
                "insert and must not be reverse-complemented again."
            ),
        },
        "coordinates": {
            "semantics": "zero-based half-open BED-like",
            "assemblies": ["hg19", "hg38"],
            "column_groups": {
                "hg19": ["hg19_chrom", "hg19_start", "hg19_end", "hg19_strand"],
                "hg38": ["hg38_chrom", "hg38_start", "hg38_end", "hg38_strand"],
            },
            "chromosome_names": "UCSC chr-prefixed",
            "exceptional_source_interval_widths": {
                "semantics": "preserved without rewriting",
                "hg19_total": 0,
                "hg38_total": 345,
                "categories": {
                    "negative control (Inoue et al 2019)": 99,
                    "positive control (Inoue et al 2019)": 49,
                    "potential enhancer": 38,
                    "potential enhancer, HepG2": 6,
                    "promoter": 153,
                },
            },
            "strand": "source genomic strand",
            "construct_orientation_column": "orientation",
            "assembly_limit": (
                "Only source labels hg19 and hg38 are claimed; no patch or "
                "sequence accession is inferred."
            ),
        },
        "identifiers": {
            "construct": "name",
            "orientation_collapsed_element": "element_id",
            "category": "category",
        },
        "targets": {
            "large_scale": {
                "primary": "observed_log2_rna_dna",
                "fold": "fold",
                "activity_call": "active",
                "source_mean": "mean_log2_rna_dna",
            },
            "joint": {
                "primary": [
                    "hepg2_observed_element_specificity_score",
                    "k562_observed_element_specificity_score",
                    "wtc11_observed_element_specificity_score",
                ],
                "source_element_specificity": [
                    "hepg2_element_specificity_score",
                    "k562_element_specificity_score",
                    "wtc11_element_specificity_score",
                ],
                "fold": "fold",
                "activity_calls": [
                    "hepg2_active",
                    "k562_active",
                    "wtc11_active",
                ],
            },
            "published_performance": {
                "metric": "pearson_r",
                "grouping": ["dataset", "cell_line", "test_fold", "model"],
                "rows": 330,
            },
        },
        "cohorts": {
            "construct_rows_across_four_libraries": 543_629,
            "joint_potential_enhancers_complete_activity_labels": 54_000,
            "joint_potential_enhancers_complete_activity_labels_and_fold": 53_151,
        },
        "verification": {
            "all_output_columns": (
                "Independently reconstructed from the six source workbooks; "
                "S11 element-specificity scores are copied by exact source "
                "identifier because their upstream normalization is unpublished."
            )
        },
        "outputs": list(OUTPUT_PATHS),
    }


def _arrow_column_sha256(table: pa.Table, column: str) -> str:
    values = table.column(column).combine_chunks()
    if pa.types.is_string(values.type) or pa.types.is_large_string(values.type):
        replacement = ""
    elif pa.types.is_floating(values.type):
        replacement = 0.0
    elif pa.types.is_integer(values.type):
        replacement = 0
    elif pa.types.is_boolean(values.type):
        replacement = False
    else:
        raise RuntimeError(f"{column}: unsupported semantic fingerprint type {values.type}")
    column_table = pa.table(
        [
            pc.is_null(values),
            pc.fill_null(values, pa.scalar(replacement, type=values.type)),
        ],
        names=[f"{column}__is_null", column],
    )
    sink = pa.BufferOutputStream()
    with pa.ipc.new_stream(sink, column_table.schema) as writer:
        writer.write_table(column_table)
    return hashlib.sha256(sink.getvalue()).hexdigest()


def build_source_output_semantics(
    frames: dict[str, Any],
) -> dict[str, Any]:
    require(set(frames) == set(OUTPUT_PATHS), "semantic output membership changed")
    outputs = {}
    for relative in OUTPUT_PATHS:
        table = pa.Table.from_pandas(frames[relative], preserve_index=False)
        row_count = table.num_rows
        vector_indices = sorted({0, row_count // 2, row_count - 1})
        outputs[relative] = {
            "rows": row_count,
            "ordered_columns": table.column_names,
            "column_sha256": {
                column: _arrow_column_sha256(table, column)
                for column in table.column_names
            },
            "test_vectors": [
                {
                    "row_index": index,
                    "fields": table.slice(index, 1).to_pylist()[0],
                }
                for index in vector_indices
            ],
        }
    return {
        "fingerprint_algorithm": SEMANTIC_FINGERPRINT_ALGORITHM,
        "outputs": outputs,
    }


def validate_source_output_semantics(
    contract: Any,
    output_records: dict[str, Any],
) -> None:
    require(isinstance(contract, dict), "source/output semantic contract must be an object")
    require_exact_keys(
        contract,
        {"fingerprint_algorithm", "outputs"},
        "source/output semantic contract",
    )
    require(
        contract["fingerprint_algorithm"] == SEMANTIC_FINGERPRINT_ALGORITHM,
        "source/output semantic fingerprint algorithm changed",
    )
    outputs = contract["outputs"]
    require(
        isinstance(outputs, dict) and set(outputs) == set(OUTPUT_PATHS),
        "source/output semantic output membership changed",
    )
    for relative in OUTPUT_PATHS:
        declaration = outputs[relative]
        require_exact_keys(
            declaration,
            {"rows", "ordered_columns", "column_sha256", "test_vectors"},
            f"{relative} semantic contract",
        )
        expected_columns = [
            field["name"]
            for field in output_records[relative]["schema"]
        ]
        require(
            declaration["rows"] == output_records[relative]["rows"],
            f"{relative}: semantic row count changed",
        )
        require(
            declaration["ordered_columns"] == expected_columns,
            f"{relative}: semantic ordered columns changed",
        )
        hashes = declaration["column_sha256"]
        require(
            isinstance(hashes, dict) and set(hashes) == set(expected_columns),
            f"{relative}: semantic column fingerprint membership changed",
        )
        for column, checksum in hashes.items():
            require(
                re.fullmatch(r"[0-9a-f]{64}", checksum or "") is not None,
                f"{relative} {column}: invalid semantic SHA-256",
            )
        vectors = declaration["test_vectors"]
        require(
            isinstance(vectors, list) and vectors,
            f"{relative}: semantic test vectors missing",
        )
        for vector in vectors:
            require_exact_keys(
                vector,
                {"row_index", "fields"},
                f"{relative} semantic test vector",
            )
            require(
                isinstance(vector["row_index"], int)
                and 0 <= vector["row_index"] < declaration["rows"],
                f"{relative}: invalid semantic test-vector row",
            )
            require(
                isinstance(vector["fields"], dict)
                and set(vector["fields"]) == set(expected_columns),
                f"{relative}: semantic test-vector fields changed",
            )


def assert_source_output_semantics(
    frames: dict[str, Any],
    contract: dict[str, Any],
) -> None:
    require(set(frames) == set(OUTPUT_PATHS), "semantic frame membership changed")
    for relative in OUTPUT_PATHS:
        assert_output_semantics(
            relative,
            frames[relative],
            contract["outputs"][relative],
        )


def assert_output_semantics(
    relative: str,
    frame: Any,
    declaration: dict[str, Any],
) -> None:
    table = pa.Table.from_pandas(frame, preserve_index=False)
    require(
        table.num_rows == declaration["rows"]
        and table.column_names == declaration["ordered_columns"],
        f"{relative}: source-derived semantic shape changed",
    )
    for vector in declaration["test_vectors"]:
        observed = table.slice(vector["row_index"], 1).to_pylist()[0]
        require(
            observed == vector["fields"],
            f"{relative}: fixed source-to-output test vector changed at "
            f"row {vector['row_index']}",
        )
    for column, expected in declaration["column_sha256"].items():
        require(
            _arrow_column_sha256(table, column) == expected,
            f"{relative} {column}: source-derived semantic fingerprint changed",
        )


def build_contracts(
    source_output_semantics: dict[str, Any],
    root: Path = ROOT,
) -> dict[str, Any]:
    card = root / HF_FILES["README.md"]
    return {
        "environment": {
            "name": ENVIRONMENT_NAME,
            "platform": ENVIRONMENT_PLATFORM,
            "lock_file": ENVIRONMENT_LOCK,
            "lock_sha256": sha256(root / ENVIRONMENT_LOCK),
            "explicit_package_count": len(explicit_urls(root / ENVIRONMENT_LOCK)),
            "direct_packages": DIRECT_PACKAGES,
        },
        "processing": processing_record(root),
        "task": task_contract(),
        "hugging_face": {
            "allowlist": list(HF_FILES),
            "files": {
                "README.md": {
                    "source": HF_FILES["README.md"],
                    "bytes": card.stat().st_size,
                    "sha256": sha256(card),
                },
                **{
                    name: {"source": source}
                    for name, source in HF_FILES.items()
                    if name != "README.md"
                },
            },
        },
        SOURCE_OUTPUT_SEMANTICS: source_output_semantics,
    }


def build_manifest(
    output_dir: Path,
    source_output_semantics: dict[str, Any],
    root: Path = ROOT,
) -> dict[str, Any]:
    outputs = {
        path: describe_parquet(output_dir / Path(path).name)
        for path in OUTPUT_PATHS
    }
    return {
        "manifest_version": MANIFEST_VERSION,
        "sources": source_specification(),
        "outputs": outputs,
        "contracts": build_contracts(source_output_semantics, root),
    }


def load_manifest(
    path: Path = MANIFEST_PATH,
    *,
    require_current_contracts: bool = True,
) -> dict[str, Any]:
    require(path.is_file(), f"missing manifest: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    validate_manifest(
        payload,
        require_current_contracts=require_current_contracts,
    )
    return payload


def validate_manifest(
    payload: Any,
    root: Path = ROOT,
    *,
    require_current_contracts: bool = True,
) -> None:
    require(isinstance(payload, dict), "manifest root must be an object")
    require_exact_keys(
        payload,
        {"manifest_version", "sources", "outputs", "contracts"},
        "manifest",
    )
    require(payload["manifest_version"] == MANIFEST_VERSION, "manifest version changed")
    require(
        payload["sources"] == source_specification(),
        "manifest source declarations changed",
    )
    require(
        set(payload["outputs"]) == set(OUTPUT_PATHS),
        "manifest output allowlist changed",
    )
    for path, record in payload["outputs"].items():
        require_exact_keys(record, {"bytes", "rows", "sha256", "schema"}, path)
        require(isinstance(record["schema"], list), f"{path}: schema must be a list")
    contracts = payload["contracts"]
    require(
        isinstance(contracts, dict)
        and set(contracts)
        == {"environment", "processing", "task", "hugging_face", SOURCE_OUTPUT_SEMANTICS},
        "manifest contract membership changed",
    )
    validate_source_output_semantics(
        contracts[SOURCE_OUTPUT_SEMANTICS],
        payload["outputs"],
    )
    if require_current_contracts:
        require(
            contracts
            == build_contracts(
                contracts[SOURCE_OUTPUT_SEMANTICS],
                root,
            ),
            "manifest contracts do not match the current environment, processing, "
            "task, and Hugging Face declarations; use explicit manifest refresh",
        )


def assert_flat_regular_files(
    directory: Path,
    expected: set[str],
    label: str,
) -> None:
    require(directory.is_dir(), f"missing {label} directory: {directory}")
    actual = {path.name for path in directory.iterdir()}
    require(
        actual == expected,
        f"{directory}: {label} membership changed; "
        f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}",
    )
    for filename in expected:
        path = directory / filename
        require(
            path.is_file() and not path.is_symlink(),
            f"{path}: {label} entry must be a regular file",
        )


def relative_payload_files(directory: Path) -> set[str]:
    require(directory.is_dir(), f"not a directory: {directory}")
    files = set()
    for path in directory.rglob("*"):
        relative = path.relative_to(directory)
        if relative.parts[0] == ".git":
            require(not path.is_symlink(), f"release metadata contains symlink: {path}")
            continue
        if relative.as_posix() == ".gitattributes":
            require(
                path.is_file() and not path.is_symlink(),
                f"release metadata must be a regular file: {path}",
            )
            continue
        require(not path.is_symlink(), f"release payload contains symlink: {path}")
        require(path.is_file(), f"release payload contains directory: {path}")
        files.add(relative.as_posix())
    return files


def assert_hf_release_tree(directory: Path, manifest: dict[str, Any]) -> None:
    expected = set(HF_FILES)
    actual = relative_payload_files(directory)
    require(
        actual == expected,
        f"HF payload mismatch: missing={sorted(expected - actual)}, "
        f"extra={sorted(actual - expected)}",
    )
    contract = manifest["contracts"]["hugging_face"]
    require(contract["allowlist"] == list(HF_FILES), "HF allowlist changed")
    require(set(contract["files"]) == expected, "HF file membership changed")
    for name, declaration in contract["files"].items():
        path = directory / name
        if name == "README.md":
            expected_record = {
                "bytes": declaration["bytes"],
                "sha256": declaration["sha256"],
            }
        else:
            expected_record = manifest["outputs"][declaration["source"]]
        require(path.stat().st_size == expected_record["bytes"], f"{name}: bytes changed")
        require(sha256(path) == expected_record["sha256"], f"{name}: checksum changed")
