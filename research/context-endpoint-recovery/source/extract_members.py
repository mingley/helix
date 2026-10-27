#!/usr/bin/env python3
"""Inventory a pinned source ZIP and extract only explicitly named members.

This tool never opens a workbook or interprets quantitative cells. Protocol
selection and semantic authorization remain the caller's responsibility.
"""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import tempfile
import zipfile


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def safe_member(info):
    name = info.orig_filename
    path = PurePosixPath(name)
    mode = info.external_attr >> 16
    if (
        not name
        or "\x00" in name
        or "\\" in name
        or path.is_absolute()
        or any(part in (".", "..", "") for part in name.split("/"))
        or re.match(r"^[A-Za-z]:", name)
        or info.is_dir()
        or stat.S_ISLNK(mode)
        or (stat.S_IFMT(mode) and not stat.S_ISREG(mode))
        or info.flag_bits & 1
    ):
        raise ValueError(f"Unsafe or unsupported selected member: {name!r}")
    return path


def inventory(archive):
    return [
        {
            "index": index,
            "name": info.orig_filename,
            "directory": info.is_dir(),
            "compressed_bytes": info.compress_size,
            "uncompressed_bytes": info.file_size,
            "compression_method": info.compress_type,
            "crc32_hex": f"{info.CRC:08x}",
            "local_header_offset": info.header_offset,
            "flag_bits": info.flag_bits,
            "external_attributes": info.external_attr,
        }
        for index, info in enumerate(archive.infolist())
    ]


def extract(archive, names, output, max_member, max_total):
    if len(names) != len(set(names)):
        raise ValueError("Repeated selection name")
    selected = []
    for name in names:
        matches = [info for info in archive.infolist() if info.orig_filename == name]
        if len(matches) != 1:
            raise ValueError(f"Selection must match exactly one literal member: {name!r}")
        info = matches[0]
        relative = safe_member(info)
        if info.file_size > max_member:
            raise ValueError(f"Selected member exceeds byte bound: {name!r}")
        selected.append((info, relative))
    if sum(info.file_size for info, _ in selected) > max_total:
        raise ValueError("Selected members exceed total uncompressed byte bound")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    receipts = []
    for info, relative in selected:
        target = output.joinpath(*relative.parts)
        # Check each existing ancestor before making any child directory. A
        # single mkdir(parents=True) could create directories through a symlink
        # before a later resolve/confinement check noticed it.
        parent = output
        for part in relative.parts[:-1]:
            parent = parent / part
            if parent.is_symlink():
                raise ValueError(f"Symbolic-link output ancestor: {parent}")
            if parent.exists():
                if not parent.is_dir():
                    raise ValueError(f"Non-directory output ancestor: {parent}")
            else:
                parent.mkdir()
        if target.parent.resolve() != target.parent or target.exists() or target.is_symlink():
            raise ValueError(f"Existing output or symbolic-link parent: {target}")
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=target.parent, prefix=".extract-", delete=False) as handle:
                temporary = Path(handle.name)
                size = 0
                digest = hashlib.sha256()
                with archive.open(info) as source:
                    for chunk in iter(lambda: source.read(1024 * 1024), b""):
                        size += len(chunk)
                        if size > info.file_size or size > max_member:
                            raise ValueError("Extracted bytes exceed declared/bounded size")
                        digest.update(chunk)
                        handle.write(chunk)
                if size != info.file_size:
                    raise ValueError("Extracted member size mismatch")
            # ZipExtFile validates CRC at EOF; link publishes without overwrite.
            os.link(temporary, target)
            receipts.append({
                "literal_member_name": info.orig_filename,
                "output_path": str(target),
                "uncompressed_bytes": size,
                "compressed_bytes": info.compress_size,
                "zip_crc32_hex": f"{info.CRC:08x}",
                "crc32_verified_by_zipfile": True,
                "sha256": digest.hexdigest(),
                "workbook_or_quantitative_cells_inspected": False,
            })
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    return receipts


def self_test():
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        source = directory / "source.zip"
        with zipfile.ZipFile(source, "w") as archive:
            archive.writestr("tables/Table 10.xlsx", b"synthetic payload")
            archive.writestr("../escape.xlsx", b"bad")
            link = zipfile.ZipInfo("link.xlsx")
            link.create_system = 3
            link.external_attr = (stat.S_IFLNK | 0o777) << 16
            archive.writestr(link, b"target")
            archive.writestr("unsafe/newdir/member.xlsx", b"must not create outside")
        with zipfile.ZipFile(source) as archive:
            assert len(inventory(archive)) == 4
            receipts = extract(archive, ["tables/Table 10.xlsx"], directory / "out", 100, 100)
            assert receipts[0]["sha256"] == hashlib.sha256(b"synthetic payload").hexdigest()
            for names, limit in [(["../escape.xlsx"], 100), (["link.xlsx"], 100), (["missing.xlsx"], 100), (["tables/Table 10.xlsx"], 2)]:
                try:
                    extract(archive, names, directory / "reject", limit, 100)
                except ValueError:
                    continue
                raise AssertionError(f"Expected rejection: {names}")
            outside = directory / "outside"
            outside.mkdir()
            reject = directory / "ancestor-check"
            reject.mkdir()
            (reject / "unsafe").symlink_to(outside, target_is_directory=True)
            try:
                extract(archive, ["unsafe/newdir/member.xlsx"], reject, 100, 100)
            except ValueError:
                pass
            else:
                raise AssertionError("Symbolic-link output ancestor was accepted")
            assert not (outside / "newdir").exists()
        damaged = directory / "damaged.zip"
        raw = source.read_bytes().replace(b"synthetic payload", b"synthetic payloaX", 1)
        damaged.write_bytes(raw)
        try:
            with zipfile.ZipFile(damaged) as archive:
                extract(archive, ["tables/Table 10.xlsx"], directory / "crc", 100, 100)
        except zipfile.BadZipFile:
            pass
        else:
            raise AssertionError("Corrupted CRC was accepted")
        assert not (directory / "crc/tables/Table 10.xlsx").exists()
    print("Synthetic extraction, traversal, symlink member/ancestor, missing-selection, byte-bound and CRC checks passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--expected-archive-sha256")
    parser.add_argument("--protocol", type=Path)
    parser.add_argument("--protocol-sha256")
    parser.add_argument("--inventory-json", type=Path)
    parser.add_argument("--select", action="append", default=[])
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--receipt-json", type=Path)
    parser.add_argument("--max-member-bytes", type=int, default=64 * 1024 * 1024)
    parser.add_argument("--max-total-bytes", type=int, default=192 * 1024 * 1024)
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not all((args.archive, args.expected_archive_sha256, args.protocol, args.protocol_sha256, args.inventory_json)):
        parser.error("Archive, expected archive SHA, protocol, protocol SHA and inventory path required")
    if sha256(args.protocol) != args.protocol_sha256:
        raise ValueError("Frozen protocol hash mismatch")
    observed = sha256(args.archive)
    if observed != args.expected_archive_sha256:
        raise ValueError("Full source archive hash mismatch")
    if min(args.max_member_bytes, args.max_total_bytes) <= 0:
        raise ValueError("Positive byte bounds required")
    if args.select and not all((args.output_dir, args.receipt_json)):
        parser.error("Explicit extraction output and receipt paths required")
    receipt = {"started_at_utc": now(), "protocol_path": str(args.protocol.resolve()), "protocol_sha256": args.protocol_sha256, "archive_path": str(args.archive.resolve()), "archive_bytes": args.archive.stat().st_size, "archive_sha256": observed, "max_member_bytes": args.max_member_bytes, "max_total_bytes": args.max_total_bytes, "scope": "ZIP metadata and exact-literal member extraction only; no workbook parsing"}
    with zipfile.ZipFile(args.archive) as archive:
        entries = inventory(archive)
        with args.inventory_json.open("x") as handle:
            json.dump({**receipt, "members": entries}, handle, indent=2)
            handle.write("\n")
        if args.select:
            receipt["selected_members"] = extract(archive, args.select, args.output_dir, args.max_member_bytes, args.max_total_bytes)
            receipt["completed_at_utc"] = now()
            with args.receipt_json.open("x") as handle:
                json.dump(receipt, handle, indent=2)
                handle.write("\n")
    print(json.dumps({"inventory_members": len(entries), "extracted_members": len(args.select), "archive_sha256": observed}))


if __name__ == "__main__":
    main()
