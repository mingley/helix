#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import sys
import urllib.request
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import assert_flat_regular_files, load_manifest, require, sha256


def verify(path: Path, record: dict) -> None:
    require(path.is_file(), f"missing source file: {path}")
    require(
        path.stat().st_size == record["bytes"],
        f"{path}: expected {record['bytes']} bytes, got {path.stat().st_size}",
    )
    observed = sha256(path)
    require(
        observed == record["sha256"],
        f"{path}: expected sha256 {record['sha256']}, got {observed}",
    )


def prepare_archive(path: Path, record: dict, reuse: Path | None) -> None:
    if path.is_file():
        verify(path, record)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    if reuse is not None:
        reuse = reuse.expanduser().resolve()
        verify(reuse, record)
        path.symlink_to(reuse)
        return
    temporary = path.with_suffix(path.suffix + ".tmp")
    request = urllib.request.Request(
        record["url"],
        headers={"User-Agent": "LentiMPRA-Enhancer-Processing/1.0"},
    )
    with urllib.request.urlopen(request, timeout=300) as source:
        with temporary.open("wb") as target:
            shutil.copyfileobj(source, target)
    verify(temporary, record)
    temporary.replace(path)


def extract(archive: Path, output: Path, workbooks: dict[str, dict]) -> None:
    output.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as source:
        for record in workbooks.values():
            target = output / record["filename"]
            if target.is_file():
                verify(target, record)
                continue
            temporary = target.with_suffix(target.suffix + ".tmp")
            with source.open(record["archive_member"]) as member:
                with temporary.open("wb") as handle:
                    shutil.copyfileobj(member, handle)
            verify(temporary, record)
            temporary.replace(target)
    assert_flat_regular_files(
        output,
        {record["filename"] for record in workbooks.values()},
        "source workbook",
    )


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=root / "data")
    parser.add_argument("--reuse-archive", type=Path)
    parser.add_argument("--manifest", type=Path, default=root / "manifest.json")
    args = parser.parse_args()
    data_root = args.data_root.expanduser().resolve()
    manifest = load_manifest(args.manifest)
    sources = manifest["sources"]

    archive_record = sources["archive"]
    archive = data_root / archive_record["data_path"]
    prepare_archive(archive, archive_record, args.reuse_archive)
    extract(archive, data_root / "raw", sources["workbooks"])
    print(f"prepared six manifest-pinned supplementary workbooks under {data_root}")


if __name__ == "__main__":
    main()
