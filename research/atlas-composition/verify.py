#!/usr/bin/env python3
"""Check primary-chromosome composition from raw UCSC genePred and twoBit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from twobit_reference import TwoBitReader

ATLAS = ROOT / "datasets" / "genome-atlas"
GENEPRED = "ncbiRefSeqCurated.txt.gz"
GENEPRED_SHA256 = "67090326b606658bdabdfe5f8dcb97b7943ac802c7945399522be7174ec021ca"
BASES = "ACGTN"
CHROMOSOMES = tuple(f"chr{i}" for i in range(1, 23)) + ("chrX", "chrY", "chrM")
CHUNK_BASES = 1 << 20
UNION_WIDTHS = ("transcript_span_union_bases", "exon_union_bases", "cds_union_bases")
FROZEN_LENGTHS = {
    "cds": 34_561_429,
    "exon_outside_cds": 58_522_466,
    "transcript_span_outside_exon": 1_328_191_871,
    "outside_transcript_span": 1_667_010_635,
}


def file_sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def require(condition, description):
    if not condition:
        raise ValueError(description)


def merged(spans):
    """Union of one coverage layer on the forward reference coordinate axis."""
    result = []
    for start, end in sorted(spans):
        if result and start <= result[-1][1]:
            if end > result[-1][1]:
                result[-1] = (result[-1][0], end)
        else:
            result.append((start, end))
    return result


def read_genepred(path, sizes):
    # Retain unique coordinates. Isoforms, genes and strands may overlap.
    geometry = {chrom: (set(), set(), set()) for chrom in CHROMOSOMES}
    retained = 0
    with gzip.open(path, "rt", encoding="utf-8") as source:
        for row_number, line in enumerate(source, 1):
            fields = line.rstrip("\r\n").split("\t")
            require(len(fields) == 16, f"invalid genePred row {row_number}")
            chrom = fields[2]
            if chrom not in geometry:
                continue
            tx_start, tx_end, cds_start, cds_end, n_exons = map(int, fields[4:9])
            starts = [int(token) for token in fields[9].rstrip(",").split(",")]
            ends = [int(token) for token in fields[10].rstrip(",").split(",")]
            require(0 <= tx_start < tx_end <= sizes[chrom]
                    and tx_start <= cds_start <= cds_end <= tx_end
                    and len(starts) == len(ends) == n_exons and n_exons > 0,
                    f"invalid transcript geometry at genePred row {row_number}")
            spans, exons, cds = geometry[chrom]
            spans.add((tx_start, tx_end))
            previous_end = tx_start
            for start, end in zip(starts, ends):
                require(previous_end <= start < end <= tx_end,
                        f"invalid exon geometry at genePred row {row_number}")
                previous_end = end
                exons.add((start, end))
                left, right = max(start, cds_start), min(end, cds_end)
                if left < right:
                    cds.add((left, right))
            retained += 1
    return geometry, retained


def add_counts(totals, sequence, left, right):
    for i, base in enumerate(BASES):
        totals[i] += sequence.count(base, left, right)


def chromosome_counts(reference, chrom, length, raw_layers, chunk_size=CHUNK_BASES):
    # Each union is counted independently, then nested cumulative counts are
    # subtracted per nucleotide. DNA is decoded once per bounded window.
    layers = [merged(spans) for spans in raw_layers]
    cursors = [0, 0, 0]
    counts = [[0] * len(BASES) for _ in layers]
    whole = [0] * len(BASES)
    for start in range(0, length, chunk_size):
        end = min(length, start + chunk_size)
        sequence = reference.fetch(chrom, start, end)
        add_counts(whole, sequence, 0, len(sequence))
        for level, intervals in enumerate(layers):
            index = cursors[level]
            while index < len(intervals):
                left, right = intervals[index]
                if left >= end:
                    break
                add_counts(counts[level], sequence, max(start, left) - start,
                           min(end, right) - start)
                if right > end:
                    break
                index += 1
            cursors[level] = index

    require(sum(whole) == length, f"{chrom}: reference alphabet or length differs")
    for observed, spans in zip(counts, layers):
        require(sum(observed) == sum(end - start for start, end in spans),
                f"{chrom}: incomplete reference count for annotation union")
    return whole, counts


def trace_exon_ns(reference, chrom, raw_exons):
    positions = []
    for left, right in merged(raw_exons):
        for start in range(left, right, CHUNK_BASES):
            sequence = reference.fetch(chrom, start, min(right, start + CHUNK_BASES))
            positions.extend(start + match.start() for match in re.finditer("N", sequence))
    return positions


def attach_raw_exons(path, sites):
    by_chrom = {}
    for site in sites:
        by_chrom.setdefault(site["chromosome"], []).append(site)
    if not by_chrom:
        return
    seen_rows = set()
    with gzip.open(path, "rt", encoding="utf-8") as source:
        for row_number, line in enumerate(source, 1):
            fields = line.rstrip("\r\n").split("\t")
            chrom = fields[2]
            if chrom not in by_chrom or line in seen_rows:
                continue
            tx_start, tx_end, cds_start, cds_end = map(int, fields[4:8])
            enclosed = [site for site in by_chrom[chrom] if tx_start <= site["start0"] < tx_end]
            if not enclosed:
                continue
            # The pinned UCSC genePred includes the bin column. name2 carries
            # the source gene symbol in this RefSeq track.
            starts = map(int, fields[9].rstrip(",").split(","))
            ends = map(int, fields[10].rstrip(",").split(","))
            for start, end in zip(starts, ends):
                for site in enclosed:
                    position = site["start0"]
                    if start <= position < end:
                        site["raw_exons"].append({
                            "source_row_1based": row_number,
                            "transcript_accession": fields[1],
                            "source_gene_symbol": fields[12],
                            "strand": fields[3],
                            "exon_start0": start,
                            "exon_end0": end,
                            "site_in_cds": cds_start <= position < cds_end,
                        })
            # Distinct records sharing an accession or exon retain their rows.
            seen_rows.add(line)
    for site in sites:
        require(site["raw_exons"], f"no raw genePred exon for {site['chromosome']}:{site['start0']}")


def exclusive_counts(whole, layers):
    spans, exons, cds = layers
    def minus(outer, inner):
        return [a - b for a, b in zip(outer, inner)]

    return {
        "cds": dict(zip("acgtn", cds)),
        "exon_outside_cds": dict(zip("acgtn", minus(exons, cds))),
        "transcript_span_outside_exon": dict(zip("acgtn", minus(spans, exons))),
        "outside_transcript_span": dict(zip("acgtn", minus(whole, spans))),
    }


def nuclear_called_gc_ordering(by_chrom):
    class_order = ("cds", "exon_outside_cds", "transcript_span_outside_exon",
                   "outside_transcript_span")
    nuclear = [chrom for chrom in CHROMOSOMES if chrom != "chrM"]
    ordered = 0
    gaps = []
    for chrom in nuclear:
        fractions = {}
        for label in class_order:
            counts = by_chrom[chrom][label]
            called = sum(counts[base] for base in "acgt")
            require(called > 0, f"{chrom}/{label}: undefined GC fraction among called bases")
            fractions[label] = Fraction(counts["g"] + counts["c"], called)
        adjacent = [(fractions[first] - fractions[second], chrom, first, second)
                    for first, second in zip(class_order, class_order[1:])]
        ordered += all(gap > 0 for gap, _, _, _ in adjacent)
        gaps.extend(adjacent)

    gap, chrom, first, second = min(gaps, key=lambda entry: entry[0])
    return {
        "class_order": list(class_order),
        "chromosomes_compared": len(nuclear),
        "chromosomes_in_strict_descending_order": ordered,
        "smallest_adjacent_difference": {
            "chromosome": chrom,
            "first_class": first,
            "second_class": second,
            # Keep the reduced fraction exact for consumers with 53-bit JSON integers.
            "numerator": str(gap.numerator),
            "denominator": str(gap.denominator),
        },
    }


def compare_classes(label, observed, expected):
    require(observed.keys() == expected.keys(), f"{label}: annotation class names differ")
    for key, counts in expected.items():
        require(observed[key].keys() == counts.keys() and
                all(type(n) is int and n >= 0 for n in observed[key].values()),
                f"{label}/{key}: expected nonnegative integer A/C/G/T/N counts")
        require(observed[key] == counts,
                f"{label}/{key}: raw-source counts {counts}, Rust supplied {observed[key]}")


def verify(reference_path, report_path):
    report = json.loads(report_path.read_bytes())
    require(report["schema_version"] == 1 and report["operation"] == "atlas_composition",
            "unknown atlas composition report")

    manifest_bytes = (ATLAS / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    require(tuple(manifest["primary_chromosomes"]) == CHROMOSOMES,
            "atlas primary chromosome order differs")
    source_hash = file_sha256(ATLAS / GENEPRED)
    require(source_hash == manifest["files"][GENEPRED] == GENEPRED_SHA256,
            "frozen raw genePred checksum differs")
    binding = {field: manifest[field] for field in ("assembly", "annotation", "coordinate_system")}
    binding.update(manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),
                   catalog_sha256=manifest["files"]["genes.json"])
    require(report["atlas"] == binding, "reported atlas identity differs from the frozen manifest")

    summary_bytes = (ATLAS / "summary.json").read_bytes()
    require(hashlib.sha256(summary_bytes).hexdigest() == manifest["files"]["summary.json"],
            "frozen coverage summary checksum differs")
    summary = json.loads(summary_bytes)
    annotation_rows = summary["chromosomes"]
    report_rows = report["chromosomes"]
    for rows in (annotation_rows, report_rows):
        require([row["name"] for row in rows] == list(CHROMOSOMES),
                "primary chromosome coverage or order differs")

    identity = json.loads((ROOT / "research/context-reference/hg38-reference.json").read_bytes())
    reference_hash = file_sha256(reference_path)
    require(reference_hash == identity["sha256"], "hg38 reference checksum differs")
    with TwoBitReader(reference_path) as reference:
        expected_reference = {
            "format": "UCSC twoBit", "format_version": 0,
            "file_bytes": identity["bytes"], "byte_order": reference.byte_order,
            "declared_assembly": "GRCh38", "sha256": reference_hash,
        }
        require(reference.file_bytes == identity["bytes"]
                and report["reference"] == expected_reference, "reported reference identity differs")
        sizes = reference.sequence_sizes
        geometry, retained = read_genepred(ATLAS / GENEPRED, sizes)
        require(retained == summary["retained_transcript_count"],
                "retained raw model count differs from the frozen atlas")

        total = {label: dict.fromkeys("acgtn", 0) for label in FROZEN_LENGTHS}
        exonic_n_sites = []
        independent_counts = {}
        for chrom, width_row, row in zip(CHROMOSOMES, annotation_rows, report_rows):
            length = sizes[chrom]
            require(type(row["length"]) is int and length == row["length"] == width_row["length"],
                    f"{chrom}: chromosome length differs")
            raw_layers = geometry.pop(chrom)
            whole, layers = chromosome_counts(reference, chrom, length, raw_layers)
            exon_ns = layers[1][BASES.index("N")]
            if exon_ns:
                positions = trace_exon_ns(reference, chrom, raw_layers[1])
                require(len(positions) == exon_ns, f"{chrom}: localized exonic N count differs")
                exonic_n_sites.extend({"chromosome": chrom, "start0": pos, "end0": pos + 1,
                                       "reference_base": "N", "raw_exons": []} for pos in positions)
            for counts, width in zip(layers, UNION_WIDTHS):
                require(sum(counts) == width_row[width], f"{chrom}: {width} differs")
            expected = exclusive_counts(whole, layers)
            independent_counts[chrom] = expected
            compare_classes(chrom, row["classes"], expected)
            for category, counts in expected.items():
                for base, count in counts.items():
                    total[category][base] += count

    compare_classes("all primary chromosomes", report["total"], total)
    attach_raw_exons(ATLAS / GENEPRED, exonic_n_sites)
    coding_sites = sum(any(exon["site_in_cds"] for exon in site["raw_exons"])
                       for site in exonic_n_sites)
    require(coding_sites == total["cds"]["n"]
            and len(exonic_n_sites) - coding_sites == total["exon_outside_cds"]["n"],
            "raw site assignments disagree with the verified class counts")
    widths = {label: sum(counts.values()) for label, counts in total.items()}
    require(widths == FROZEN_LENGTHS, "exclusive lengths differ from the frozen atlas")
    return {
        "verified_chromosomes": len(CHROMOSOMES),
        "verified_bases": sum(widths.values()),
        "retained_raw_models": retained,
        "genepred_sha256": source_hash,
        "reference_sha256": reference_hash,
        "total": total,
        "exonic_n_sites": exonic_n_sites,
        "nuclear_called_gc_ordering": nuclear_called_gc_ordering(independent_counts),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.reference, args.report), indent=2, sort_keys=True))
