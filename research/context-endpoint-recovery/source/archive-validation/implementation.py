#!/usr/bin/env python3
"""Bounded archive-integrity stage; no extraction or scientific cell parsing.

All member bytes are streamed for CRC/SHA verification. Reading those bytes is
not an assertion of source semantics or independent experimental validation.
The active downloader is a separate, immutable implementation/receipt stage.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import struct
import zipfile
import zlib


ALLOWED_MIME = {"application/zip", "application/x-zip-compressed", "application/octet-stream"}
LOCAL = struct.Struct("<4s5H3I2H")
END = struct.Struct("<4s4H2IH")
CHUNK = 1024 * 1024


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, value):
    with Path(path).open("x") as handle:
        json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")


def valid_digest(value):
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError("explicit lowercase SHA256 required")


def safe_name(info):
    name = info.orig_filename
    parts = name.split("/")
    if info.is_dir():
        parts = parts[:-1]
    mode = info.external_attr >> 16
    kind = stat.S_IFMT(mode)
    allowed_kind = stat.S_IFDIR if info.is_dir() else stat.S_IFREG
    if (not name or not parts or "\x00" in name or "\\" in name
            or PurePosixPath(name).is_absolute()
            or any(part in ("", ".", "..") for part in parts)
            or re.match(r"^[A-Za-z]:", name)
            or kind not in (0, allowed_kind) or info.flag_bits & 0x41):
        raise ValueError(f"unsafe/encrypted/nonregular ZIP entry: {name!r}")
    if info.is_dir() and info.file_size:
        raise ValueError("directory entry has a payload")
    return "/".join(parts)


def end_record(handle, size):
    handle.seek(max(0, size - 65557))
    tail = handle.read()
    position = tail.rfind(b"PK\x05\x06")
    if position < 0 or len(tail) - position < END.size:
        raise ValueError("missing/truncated ZIP end record")
    fields = END.unpack_from(tail, position)
    _, disk, central_disk, disk_count, count, central_bytes, central_offset, comment_bytes = fields
    absolute = size - len(tail) + position
    if absolute + END.size + comment_bytes != size:
        raise ValueError("ZIP end/comment length mismatch or trailing bytes")
    if disk or central_disk or disk_count != count:
        raise ValueError("multivolume ZIP is unsupported")
    if count == 0xFFFF or central_bytes == 0xFFFFFFFF or central_offset == 0xFFFFFFFF:
        raise ValueError("ZIP64 is unsupported by this explicit integrity stage")
    if central_offset + central_bytes != absolute:
        raise ValueError("central-directory bounds, prefix or concatenated ZIP mismatch")
    return {"entries": count, "central_offset": central_offset,
            "central_bytes": central_bytes, "end_offset": absolute,
            "comment_bytes": comment_bytes}


def local_bounds(handle, info, central_offset):
    start = info.header_offset
    if start < 0 or start + LOCAL.size > central_offset:
        raise ValueError("local header is outside member area")
    handle.seek(start)
    header = handle.read(LOCAL.size)
    if len(header) != LOCAL.size:
        raise ValueError("truncated local ZIP header")
    signature, version, flags, compression, _, _, crc, compressed, size, name_bytes, extra_bytes = LOCAL.unpack(header)
    if signature != b"PK\x03\x04":
        raise ValueError("invalid local-header magic")
    if flags != info.flag_bits or compression != info.compress_type:
        raise ValueError("central/local ZIP flags or compression disagree")
    data_start = start + LOCAL.size + name_bytes + extra_bytes
    data_end = data_start + info.compress_size
    if data_start > central_offset or data_end > central_offset:
        raise ValueError("member payload exceeds central-directory boundary")
    name = handle.read(name_bytes).decode("utf-8" if flags & 0x800 else "cp437")
    if name != info.orig_filename:
        raise ValueError("central/local literal member names disagree")
    if version >= 45 or compressed == 0xFFFFFFFF or size == 0xFFFFFFFF:
        raise ValueError("ZIP64/newer extraction version unsupported by this stage")
    if flags & 8:
        # Both classic descriptor encodings are legal; the signature is optional.
        handle.seek(data_end)
        descriptor = handle.read(min(16, central_offset - data_end))
        candidates = []
        if len(descriptor) >= 12:
            candidates.append((12, struct.unpack_from("<III", descriptor)))
        if len(descriptor) >= 16 and descriptor[:4] == b"PK\x07\x08":
            candidates.append((16, struct.unpack_from("<III", descriptor, 4)))
        lengths = [length for length, values in candidates
                   if values == (info.CRC, info.compress_size, info.file_size)]
        if not lengths:
            raise ValueError("missing/inconsistent ZIP data descriptor")
        data_end += min(lengths)
    elif (crc, compressed, size) != (info.CRC, info.compress_size, info.file_size):
        raise ValueError("central/local member CRC or size disagree")
    return {"header_start": start, "data_start": data_start, "record_end": data_end}


def inspect_and_scan(path, selected=(), *, max_member_bytes=1024**3, max_total_bytes=1024**3):
    """Validate every entry's metadata/bounds, then every byte's CRC and SHA.

    This implementation intentionally supports a single classic ZIP, stored or
    deflated entries, and classic data descriptors. Unsupported variants abort;
    they are never silently skipped or interpreted as a different archive.
    """
    if min(max_member_bytes, max_total_bytes) <= 0:
        raise ValueError("positive decompression bounds required")
    if len(selected) != len(set(selected)):
        raise ValueError("duplicate selected name")
    size = path.stat().st_size
    with path.open("rb") as handle:
        if handle.read(4) != b"PK\x03\x04":
            raise ValueError("archive must begin with ZIP local-header magic, not HTML/prefix")
        end = end_record(handle, size)
        with zipfile.ZipFile(path) as archive:
            entries = archive.infolist()
            if len(entries) != end["entries"] or archive.start_dir != end["central_offset"]:
                raise ValueError("central-directory count/offset mismatch")
            if not entries:
                raise ValueError("empty archive")
            names = [entry.orig_filename for entry in entries]
            if len(names) != len(set(names)):
                raise ValueError("duplicate literal ZIP member names")
            if not set(selected).issubset(names):
                raise ValueError("selected literal member absent")
            canonical = set()
            bounds = []
            declared_total = 0
            for info in entries:
                key = safe_name(info)
                if key in canonical:
                    raise ValueError("conflicting file/directory paths")
                canonical.add(key)
                if info.compress_type not in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED):
                    raise ValueError("unsupported compression method")
                declared_total += info.file_size
                if info.file_size > max_member_bytes or declared_total > max_total_bytes:
                    raise ValueError("declared decompression budget exceeded")
                bounds.append(local_bounds(handle, info, end["central_offset"]))
            ordered = sorted(bounds, key=lambda bound: bound["header_start"])
            if any(left["record_end"] > right["header_start"]
                   for left, right in zip(ordered, ordered[1:])):
                raise ValueError("overlapping ZIP member records")
            receipts = []
            observed_total = 0
            for info, bound in zip(entries, bounds):
                count, crc, digest = 0, 0, hashlib.sha256()
                with archive.open(info) as source:
                    for chunk in iter(lambda: source.read(CHUNK), b""):
                        count += len(chunk)
                        observed_total += len(chunk)
                        if count > info.file_size or count > max_member_bytes or observed_total > max_total_bytes:
                            raise ValueError("observed decompression budget exceeded")
                        crc = zlib.crc32(chunk, crc)
                        digest.update(chunk)
                if count != info.file_size or crc != info.CRC:
                    raise ValueError("member size/CRC mismatch at EOF")
                receipts.append({"literal_member_name": info.orig_filename, "directory": info.is_dir(),
                    "compressed_bytes": info.compress_size, "uncompressed_bytes": count,
                    "crc32_hex": f"{crc:08x}", "sha256": digest.hexdigest(),
                    "crc_verified": True, "selected_for_authorized_schema_inspection": info.orig_filename in selected,
                    "scientific_values_interpreted": False, **bound})
            return {"end_record": end, "member_count": len(entries),
                "all_members_crc_verified": True, "total_uncompressed_bytes": observed_total,
                "max_member_bytes": max_member_bytes, "max_total_bytes": max_total_bytes,
                "selected_literal_names": list(selected), "members": receipts,
                "scope": "ALL members streamed for CRC/SHA, no extraction, no workbook/cell parsing"}


def verify(archive, expected_archive_sha256, protocol, protocol_sha256, acquisition,
           out, selected=(), predicted_archive_sha256=None, max_member_bytes=1024**3,
           max_total_bytes=1024**3):
    archive, protocol, acquisition, out = map(Path, (archive, protocol, acquisition, out))
    out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, out / "implementation.py")
    receipt = {"schema_version": 1, "started_at_utc": now(), "complete": False,
        "status": "failed", "archive_path": str(archive.resolve()),
        "protocol_sha256": protocol_sha256, "expected_archive_sha256": expected_archive_sha256,
        "implementation_sha256": sha256(out / "implementation.py"),
        "producer_checksum_asserted": False, "scientific_values_interpreted": False}
    try:
        for value in (expected_archive_sha256, protocol_sha256):
            valid_digest(value)
        if sha256(protocol) != protocol_sha256:
            raise ValueError("frozen protocol hash mismatch")
        contract = json.loads(protocol.read_bytes())
        shutil.copyfile(protocol,out/"protocol.json")
        source = contract["authoritative_source"]["primary_archive"]
        receipt["protocol_path"] = str(protocol.resolve())
        acquired = json.loads(acquisition.read_bytes())
        receipt["acquisition_receipt_sha256"] = sha256(acquisition)
        shutil.copyfile(acquisition,out/"acquisition.json")
        mime = acquired.get("headers", {}).get("content-type", "").split(";", 1)[0].strip().lower()
        if mime not in ALLOWED_MIME:
            raise ValueError(f"unsupported/HTML MIME: {mime!r}")
        if (acquired.get("complete") is not True or acquired.get("protocol_sha256") != protocol_sha256
                or acquired.get("url") != source["url"] or acquired.get("sha256") != expected_archive_sha256
                or acquired.get("status") != 200 or acquired.get("bytes") != source["expected_size_bytes"]):
            raise ValueError("download receipt/source contract mismatch")
        before = archive.stat()
        if before.st_size != source["expected_size_bytes"]:
            raise ValueError("archive byte count differs from frozen primary source")
        observed = sha256(archive)
        if observed != expected_archive_sha256:
            raise ValueError("archive hash differs from complete acquisition commitment")
        receipt.update(archive_bytes=before.st_size, archive_sha256=observed,
            download_content_type=acquired["headers"]["content-type"])
        if predicted_archive_sha256 is not None:
            valid_digest(predicted_archive_sha256)
            receipt["cross_source_code_prediction"] = {"sha256": predicted_archive_sha256,
                "matches_official_download": observed == predicted_archive_sha256,
                "producer_digest": False,
                "interpretation": "Frozen downstream code prediction only; a discrepancy is a provenance finding, not permission to select a different payload."}
        inventory = inspect_and_scan(archive, selected, max_member_bytes=max_member_bytes, max_total_bytes=max_total_bytes)
        after = archive.stat()
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
                after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
            raise ValueError("source archive changed while verifying")
        save(out / "inventory.json", {"archive_sha256": observed, **inventory})
        receipt.update(complete=True, status="archive_integrity_verified", all_members_crc_verified=True,
            member_count=inventory["member_count"], total_uncompressed_bytes=inventory["total_uncompressed_bytes"],
            inventory_sha256=sha256(out / "inventory.json"))
    except Exception as error:
        receipt.update(error_type=type(error).__name__, error=str(error))
        raise
    finally:
        receipt["finished_at_utc"] = now()
        save(out / "verification.json", receipt)
        names = ["implementation.py", "verification.json"]
        for name in ("protocol.json","acquisition.json","inventory.json"):
            if (out/name).exists():
                names.append(name)
        save(out / "manifest.json", {"schema_version": 1,
            "files": {name: sha256(out / name) for name in names}})
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--expected-archive-sha256", required=True)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--acquisition-receipt", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--select", action="append", default=[])
    parser.add_argument("--predicted-archive-sha256")
    parser.add_argument("--max-member-bytes", type=int, default=1024**3)
    parser.add_argument("--max-total-bytes", type=int, default=1024**3)
    args = parser.parse_args()
    print(json.dumps(verify(args.archive,args.expected_archive_sha256,args.protocol,args.protocol_sha256,
        args.acquisition_receipt,args.out,args.select,args.predicted_archive_sha256,
        args.max_member_bytes,args.max_total_bytes),indent=2))


if __name__ == "__main__":
    main()
