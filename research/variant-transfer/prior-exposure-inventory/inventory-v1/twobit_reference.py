#!/usr/bin/env python3
"""Strict, memory-mapped UCSC twoBit v0 reader using only Python's stdlib.

All interval coordinates are zero-based and half-open. The reader deliberately
does not alias chromosomes, lift assemblies, clamp coordinates or replace Ns.
It supports either file byte order, N blocks, soft masks and reverse strands.
"""
from __future__ import annotations

import argparse
from array import array
from bisect import bisect_right
from dataclasses import dataclass
import json
import mmap
import os
from pathlib import Path
import struct
import sys


SIGNATURE = 0x1A412743
_BASES = b"TCAG"
_DECODE = tuple(bytes(_BASES[(byte >> shift) & 3] for shift in (6, 4, 2, 0))
                for byte in range(256))
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class TwoBitError(ValueError):
    """Invalid or unsupported twoBit structure, without inferred repairs."""


@dataclass(frozen=True)
class _Record:
    length: int
    n_starts: array
    n_sizes: array
    mask_starts: array
    mask_sizes: array
    packed_offset: int
    record_start: int
    record_end: int


class TwoBitReader:
    """Read validated twoBit v0 records without copying whole chromosomes.

    Every record and block table is structurally validated when opening. Only
    requested packed sequence bytes are decoded. The small tables use uint32
    arrays, and the large reference remains a read-only memory mapping.
    """

    def __init__(self, path: str | os.PathLike):
        self.path = Path(path)
        self._file = self.path.open("rb")
        self._mapping = None
        try:
            self.file_bytes = os.fstat(self._file.fileno()).st_size
            if self.file_bytes < 16:
                raise TwoBitError("truncated twoBit header")
            self._mapping = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)
            self._parse()
        except BaseException:
            self.close()
            raise

    def _u32(self, offset: int, limit: int) -> tuple[int, int]:
        if offset < 0 or offset + 4 > limit:
            raise TwoBitError("truncated uint32 field")
        return struct.unpack_from(self._endian + "I", self._mapping, offset)[0], offset + 4

    def _array(self, offset: int, count: int, limit: int) -> tuple[array, int]:
        end = offset + count * 4
        if count < 0 or end > limit:
            raise TwoBitError("block array exceeds its record boundary")
        values = array("I")
        if values.itemsize != 4:
            raise TwoBitError("this platform's unsigned-int array is not uint32")
        values.frombytes(self._mapping[offset:end])
        if self.byte_order != sys.byteorder:
            values.byteswap()
        return values, end

    @staticmethod
    def _check_blocks(starts: array, sizes: array, length: int, kind: str) -> None:
        previous_end = 0
        for start, size in zip(starts, sizes):
            if size == 0 or start < previous_end or start + size > length:
                raise TwoBitError(f"{kind} blocks must be positive, sorted, disjoint and in bounds")
            previous_end = start + size

    def _parse_record(self, start: int, limit: int) -> _Record:
        length, offset = self._u32(start, limit)
        n_count, offset = self._u32(offset, limit)
        if n_count > length:
            raise TwoBitError("N block count exceeds sequence length")
        n_starts, offset = self._array(offset, n_count, limit)
        n_sizes, offset = self._array(offset, n_count, limit)
        mask_count, offset = self._u32(offset, limit)
        if mask_count > length:
            raise TwoBitError("mask block count exceeds sequence length")
        mask_starts, offset = self._array(offset, mask_count, limit)
        mask_sizes, offset = self._array(offset, mask_count, limit)
        reserved, offset = self._u32(offset, limit)
        if reserved != 0:
            raise TwoBitError("sequence reserved field must be zero")
        self._check_blocks(n_starts, n_sizes, length, "N")
        self._check_blocks(mask_starts, mask_sizes, length, "mask")
        end = offset + (length + 3) // 4
        if end > limit:
            raise TwoBitError("packed DNA exceeds sequence record boundary")
        return _Record(length, n_starts, n_sizes, mask_starts, mask_sizes,
                       offset, start, end)

    def _parse(self) -> None:
        signature = self._mapping[:4]
        if signature == struct.pack(">I", SIGNATURE):
            self._endian, self.byte_order = ">", "big"
        elif signature == struct.pack("<I", SIGNATURE):
            self._endian, self.byte_order = "<", "little"
        else:
            raise TwoBitError("twoBit signature mismatch")
        version, count, reserved = struct.unpack_from(self._endian + "III", self._mapping, 4)
        if version != 0:
            raise TwoBitError(f"unsupported twoBit version {version}; only v0 uint32 offsets are supported")
        if reserved != 0:
            raise TwoBitError("header reserved field must be zero")
        if count > (self.file_bytes - 16) // 6:
            raise TwoBitError("sequence count cannot fit in index")
        self.version = version
        offset = 16
        index = {}
        for _ in range(count):
            if offset >= self.file_bytes:
                raise TwoBitError("truncated sequence index")
            name_length = self._mapping[offset]
            offset += 1
            if name_length == 0 or offset + name_length > self.file_bytes:
                raise TwoBitError("sequence name must be nonempty and fit inside index")
            raw_name = self._mapping[offset:offset + name_length]
            if any(character < 33 or character > 126 for character in raw_name):
                raise TwoBitError("sequence names must be printable nonwhitespace ASCII")
            name = raw_name.decode("ascii")
            offset += name_length
            record_offset, offset = self._u32(offset, self.file_bytes)
            if name in index:
                raise TwoBitError(f"duplicate sequence name: {name}")
            index[name] = record_offset
        if len(set(index.values())) != len(index):
            raise TwoBitError("distinct sequences cannot share record offsets")
        physical = sorted((record_offset, name) for name, record_offset in index.items())
        records = {}
        for i, (record_offset, name) in enumerate(physical):
            if not offset <= record_offset < self.file_bytes:
                raise TwoBitError("sequence record offset points inside index or outside file")
            limit = physical[i + 1][0] if i + 1 < len(physical) else self.file_bytes
            records[name] = self._parse_record(record_offset, limit)
        self._records = {name: records[name] for name in index}

    @property
    def sequence_sizes(self) -> dict[str, int]:
        return {name: record.length for name, record in self._records.items()}

    def fetch(self, name: str, start0: int = 0, end0: int | None = None,
              *, strand: str = "+", masked: bool = False) -> str:
        """Extract an exact interval, optionally preserving UCSC lowercase masks.

        Defaults to uppercase A/C/G/T/N. Masking is applied after N blocks, so
        overlapping N+mask blocks become lowercase n when masked=True. Strand
        reversal happens last. Unknown contig names fail without aliasing.
        """
        if self._mapping is None:
            raise ValueError("twoBit reader is closed")
        if name not in self._records:
            raise KeyError(f"sequence is absent from twoBit index: {name}")
        record = self._records[name]
        if end0 is None:
            end0 = record.length
        if any(isinstance(value, bool) or not isinstance(value, int) for value in (start0, end0)):
            raise ValueError("interval coordinates must be integers")
        if not 0 <= start0 <= end0 <= record.length:
            raise ValueError("interval must satisfy 0 <= start <= end <= sequence length")
        if strand not in {"+", "-"} or not isinstance(masked, bool):
            raise ValueError("strand must be +/- and masked must be boolean")
        if start0 == end0:
            return ""
        packed_start = record.packed_offset + start0 // 4
        packed_end = record.packed_offset + (end0 + 3) // 4
        unpacked = b"".join(_DECODE[byte] for byte in self._mapping[packed_start:packed_end])
        shift = start0 % 4
        result = bytearray(unpacked[shift:shift + end0 - start0])

        def intersecting(starts, sizes):
            index = max(0, bisect_right(starts, start0) - 1)
            while index < len(starts) and starts[index] < end0:
                left = max(start0, starts[index])
                right = min(end0, starts[index] + sizes[index])
                if left < right:
                    yield left - start0, right - start0
                index += 1

        for left, right in intersecting(record.n_starts, record.n_sizes):
            result[left:right] = b"N" * (right - left)
        if masked:
            for left, right in intersecting(record.mask_starts, record.mask_sizes):
                result[left:right] = result[left:right].lower()
        sequence = result.decode("ascii")
        return sequence if strand == "+" else sequence.translate(_COMPLEMENT)[::-1]

    def metadata(self) -> dict:
        return {"format": "UCSC twoBit", "version": self.version,
                "byte_order": self.byte_order, "file_bytes": self.file_bytes,
                "sequence_count": len(self._records),
                "sequences": {name: {"length": record.length,
                                     "n_blocks": len(record.n_starts),
                                     "mask_blocks": len(record.mask_starts)}
                              for name, record in self._records.items()}}

    def close(self) -> None:
        if self._mapping is not None:
            self._mapping.close()
            self._mapping = None
        if self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self) -> "TwoBitReader":
        return self

    def __exit__(self, *_args) -> None:
        self.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--interval", nargs=3, metavar=("CHROM", "START0", "END0"))
    parser.add_argument("--strand", choices=["+", "-"], default="+")
    parser.add_argument("--masked", action="store_true")
    args = parser.parse_args()
    with TwoBitReader(args.reference) as reference:
        if args.interval:
            chromosome, start, end = args.interval
            print(reference.fetch(chromosome, int(start), int(end),
                                  strand=args.strand, masked=args.masked))
        else:
            print(json.dumps(reference.metadata(), indent=2))


if __name__ == "__main__":
    main()
