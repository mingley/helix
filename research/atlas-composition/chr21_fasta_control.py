#!/usr/bin/env python3
"""Compare the NC_000021.9 genomic FASTA to chr21 in the frozen hg38 twoBit."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from twobit_reference import TwoBitReader

ACCESSION = "NC_000021.9"
CHROMOSOME = "chr21"
IUPAC_DNA = b"ACGTRYSWKMBDHVN"
CHUNK_BYTES = 1 << 20


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source_site(source_check):
    sites = [site for site in source_check["exonic_n_sites"]
             if site["chromosome"] == CHROMOSOME]
    require(len(sites) == 1, "expected exactly one traced chr21 exonic N")
    site = sites[0]
    start = site["start0"]
    require(type(start) is int and start >= 0 and site["end0"] == start + 1
            and site["reference_base"] == "N", "invalid traced chr21 interval")
    return start


def fasta_chunks(stream):
    header = stream.readline(CHUNK_BYTES)
    label = b">" + ACCESSION.encode("ascii")
    require(header.startswith(label) and len(header) > len(label)
            and header[len(label)] in b" \t\r\n", f"expected FASTA accession {ACCESSION}")
    while not header.endswith(b"\n"):
        header = stream.readline(CHUNK_BYTES)
        require(header, "missing newline after FASTA header")

    for block in iter(lambda: stream.read(CHUNK_BYTES), b""):
        require(b">" not in block, "unexpected additional FASTA header")
        bases = block.translate(None, b"\r\n").upper()
        invalid = bases.translate(None, IUPAC_DNA)
        if invalid:
            raise ValueError(f"invalid genomic FASTA symbol 0x{invalid[0]:02x}")
        if bases:
            yield bases


def compare(fasta_path, twobit_path, source_check_path):
    source_bytes = source_check_path.read_bytes()
    source_check = json.loads(source_bytes)
    position = source_site(source_check)
    identity = json.loads((ROOT / "research/context-reference/hg38-reference.json").read_bytes())
    with twobit_path.open("rb") as raw_reference:
        reference_hash = hashlib.file_digest(raw_reference, "sha256").hexdigest()
    require(reference_hash == identity["sha256"] == source_check["reference_sha256"],
            "twoBit checksum differs from the hg38 reference or source check")

    fasta_digest = hashlib.sha256()
    twobit_digest = hashlib.sha256()
    fasta_counts = Counter()
    differences = Counter()
    fasta_context = bytearray()
    offset = 0
    with TwoBitReader(twobit_path) as reference, fasta_path.open("rb") as raw:
        require(reference.file_bytes == identity["bytes"], "twoBit size differs from hg38")
        length = reference.sequence_sizes[CHROMOSOME]
        require(position < length, "traced site lies beyond chr21")
        left, right = max(0, position - 8), min(length, position + 9)
        twobit_context = reference.fetch(CHROMOSOME, left, right)
        require(twobit_context[position - left] == "N", "traced site differs in the twoBit")
        compressed_hash = hashlib.file_digest(raw, "sha256").hexdigest()
        compressed_bytes = raw.seek(0, 2)
        raw.seek(0)
        with gzip.GzipFile(fileobj=raw, mode="rb") as stream:
            for sequence in fasta_chunks(stream):
                end = offset + len(sequence)
                require(end <= length, f"FASTA sequence exceeds {CHROMOSOME} length")
                reference_bases = reference.fetch(CHROMOSOME, offset, end).encode("ascii")
                require(len(reference_bases) == len(sequence), "truncated twoBit window")
                fasta_digest.update(sequence)
                twobit_digest.update(reference_bases)
                fasta_counts.update(sequence)
                if sequence != reference_bases:
                    differences.update((a, b) for a, b in zip(sequence, reference_bases) if a != b)
                start_overlap, end_overlap = max(left, offset), min(right, end)
                if start_overlap < end_overlap:
                    fasta_context.extend(sequence[start_overlap - offset:end_overlap - offset])
                offset = end

    require(offset == length, f"FASTA has {offset} bases, expected {length}")
    require(len(fasta_context) == right - left, "incomplete FASTA context at the traced site")
    changed = sum(differences.values())
    return {
        "schema_version": 1,
        "chromosome": CHROMOSOME,
        "source_check_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "genomic_fasta": {
            "accession": ACCESSION,
            "gzip_bytes": compressed_bytes,
            "gzip_sha256": compressed_hash,
            "uppercase_sequence_sha256": fasta_digest.hexdigest(),
            "alphabet_counts": {chr(base): fasta_counts[base] for base in IUPAC_DNA},
            "non_acgtn_bases": sum(fasta_counts[base] for base in IUPAC_DNA if base not in b"ACGTN"),
        },
        "twobit": {
            "reference_sha256": reference_hash,
            "uppercase_sequence_sha256": twobit_digest.hexdigest(),
        },
        "compared_bases": length,
        "identical_bases": length - changed,
        "differing_bases": changed,
        "differences_fasta_to_twobit": {
            f"{chr(a)}>{chr(b)}": count for (a, b), count in sorted(differences.items())
        },
        "exonic_n_site": {
            "start0": position,
            "end0": position + 1,
            "context_start0": left,
            "context_end0": right,
            "fasta_base": chr(fasta_context[position - left]),
            "twobit_base": twobit_context[position - left],
            "fasta_context": fasta_context.decode("ascii"),
            "twobit_context": twobit_context,
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path, help="gzipped primary chromosome 21 genomic FASTA")
    parser.add_argument("twobit", type=Path, help="pinned hg38 twoBit")
    parser.add_argument("source_check", type=Path, help="independent annotation-composition receipt")
    args = parser.parse_args()
    print(json.dumps(compare(args.fasta, args.twobit, args.source_check), indent=2, sort_keys=True))
