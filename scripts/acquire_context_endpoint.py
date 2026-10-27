#!/usr/bin/env python3
"""Acquire the version-linked final supplement; inspect ZIP metadata only.

The raw archive stays outside Git. No member is executed or extracted here.
Use a new output directory for every attempt, including failed downloads.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request
from zipfile import ZipFile


URL = ("https://media.springernature.com/original/springer-static/esm/"
       "art%3A10.1038%2Fs41586-024-08430-9/MediaObjects/"
       "41586_2024_8430_MOESM5_ESM.zip")
SOURCE_BYTES = 342234089
HEADERS = {"content-type", "content-length", "etag", "last-modified",
           "accept-ranges", "content-range", "date", "server"}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def acquire(protocol, expected, out):
    protocol, out = Path(protocol), Path(out)
    if digest(protocol) != expected:
        raise ValueError("frozen acquisition protocol SHA256 mismatch")
    contract = json.loads(protocol.read_bytes())
    # Bind the exact destination to text present in the frozen contract. This
    # does not infer semantic authority from an arbitrary supplied destination.
    if URL not in json.dumps(contract):
        raise ValueError("publisher URL absent from frozen protocol")
    out.mkdir(parents=True, exist_ok=False)
    receipt = {"protocol_sha256": expected, "url": URL, "method": "GET",
               "expected_bytes": SOURCE_BYTES,
               "started_at_utc": datetime.now(timezone.utc).isoformat(),
               "complete": False, "new_member_payloads_read": False,
               "publisher_checksum_available": False,
               "etag_is_asserted_digest": False}
    target = out / "final-supplement.zip.partial"
    try:
        request = urllib.request.Request(URL, headers={"Accept-Encoding": "identity"})
        with urllib.request.urlopen(request, timeout=90) as response:
            receipt["status"] = response.status
            receipt["final_url"] = response.geturl()
            receipt["headers"] = {k.lower(): v for k, v in response.headers.items()
                                  if k.lower() in HEADERS}
            if response.status != 200:
                raise ValueError("expected full-file HTTP200")
            declared = response.headers.get("Content-Length")
            if declared is not None and int(declared) != SOURCE_BYTES:
                raise ValueError("publisher byte count changed; retain failure and review")
            count = 0
            with target.open("xb") as handle:
                while block := response.read(1048576):
                    count += len(block)
                    if count > SOURCE_BYTES:
                        raise ValueError("source exceeded frozen byte count")
                    handle.write(block)
        if count != SOURCE_BYTES:
            raise ValueError("truncated source archive")
        receipt["bytes"] = count
        receipt["sha256"] = digest(target)
        # Listing reads the central directory, not quantitative member payloads.
        with ZipFile(target) as archive:
            members = [{"path": i.filename, "compressed_bytes": i.compress_size,
                        "bytes": i.file_size, "crc32": f"{i.CRC:08x}",
                        "compression": i.compress_type, "flags": i.flag_bits,
                        "header_offset": i.header_offset}
                       for i in archive.infolist()]
        save(out / "member-inventory.json", {"archive_sha256": receipt["sha256"],
                                             "members": members})
        target.rename(out / "final-supplement.zip")
        receipt["complete"] = True
    except Exception as error:
        receipt["error"] = str(error)
        if target.exists():
            receipt["retained_partial_bytes"] = target.stat().st_size
            receipt["retained_partial_sha256"] = digest(target)
        raise
    finally:
        receipt["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        receipt["implementation_sha256"] = digest(__file__)
        save(out / "acquisition.json", receipt)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    print(json.dumps(acquire(args.protocol, args.protocol_sha256, args.out), indent=2))
