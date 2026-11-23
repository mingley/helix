#!/usr/bin/env python3
"""Replay primary hg38 annotation counts with UCSC twoBit soft masking preserved."""
import argparse
import hashlib
import json
from pathlib import Path

from verify import (
    ATLAS, BASES, CHROMOSOMES, FROZEN_LENGTHS, GENEPRED, GENEPRED_SHA256, ROOT,
    TwoBitReader, chromosome_counts, compare_classes, exclusive_counts, file_sha256,
    read_genepred, require,
)

# The source reader emits A/C/G/T/N outside mask blocks and a/c/g/t/n inside.
STRATA = ("unmasked", "masked")
NUCLEOTIDES = BASES.lower()


def check_counts(label, observed, unstratified, whole, layers):
    require(observed.keys() == FROZEN_LENGTHS.keys(), f"{label}: annotation classes differ")
    for category, counts in observed.items():
        require(set(counts) == set(STRATA), f"{label}/{category}: mask strata differ")

    for index, stratum in enumerate(STRATA):
        columns = slice(index * len(BASES), (index + 1) * len(BASES))
        expected = exclusive_counts(whole[columns], [layer[columns] for layer in layers])
        reported = {category: counts[stratum] for category, counts in observed.items()}
        compare_classes(f"{label}/{stratum}", reported, expected)

    recombined = {
        category: {base: counts["masked"][base] + counts["unmasked"][base]
                   for base in NUCLEOTIDES}
        for category, counts in observed.items()
    }
    compare_classes(f"{label}/retained unstratified report", unstratified, recombined)


def verify(reference_path, report_path):
    report = json.loads(report_path.read_bytes())
    frozen_bytes = Path(__file__).with_name("summary.json").read_bytes()
    frozen = json.loads(frozen_bytes)
    for source, operation in ((report, "atlas_composition_by_soft_mask"),
                              (frozen, "atlas_composition")):
        require(type(source["schema_version"]) is int and source["schema_version"] == 1
                and source["operation"] == operation, f"expected {operation} schema version 1")
        require(tuple(row["name"] for row in source["chromosomes"]) == CHROMOSOMES,
                f"{operation}: primary chromosome coverage or order differs")

    manifest_bytes = (ATLAS / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    require(tuple(manifest["primary_chromosomes"]) == CHROMOSOMES,
            "atlas primary chromosome order differs")
    raw_hash = file_sha256(ATLAS / GENEPRED)
    require(raw_hash == GENEPRED_SHA256 == manifest["files"][GENEPRED],
            "raw genePred checksum differs")
    atlas = {key: manifest[key] for key in ("assembly", "annotation", "coordinate_system")}
    atlas.update(manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),
                 catalog_sha256=file_sha256(ATLAS / "genes.json"))
    require(atlas["catalog_sha256"] == manifest["files"]["genes.json"]
            and atlas == report["atlas"] == frozen["atlas"], "atlas identity differs")
    coverage_bytes = (ATLAS / "summary.json").read_bytes()
    require(hashlib.sha256(coverage_bytes).hexdigest() == manifest["files"]["summary.json"],
            "atlas coverage summary checksum differs")
    coverage = json.loads(coverage_bytes)

    identity = json.loads((ROOT / "research/context-reference/hg38-reference.json").read_bytes())
    reference_hash = file_sha256(reference_path)
    require(reference_hash == identity["sha256"], "hg38 reference checksum differs")
    totals = [[0] * (len(STRATA) * len(BASES)) for _ in range(4)]  # whole, spans, exons, CDS
    with TwoBitReader(reference_path) as reference:
        require(reference.file_bytes == identity["bytes"], "hg38 file size differs")
        for source in (report, frozen):
            declared = source["reference"]["declared_assembly"]
            expected = {"format": "UCSC twoBit", "format_version": reference.version,
                        "file_bytes": reference.file_bytes, "byte_order": reference.byte_order,
                        "sha256": reference_hash, "declared_assembly": declared}
            require(declared in (None, "GRCh38") and source["reference"] == expected,
                    "reference identity differs")

        sizes = reference.sequence_sizes
        geometry, retained = read_genepred(ATLAS / GENEPRED, sizes)
        require(retained == coverage["retained_transcript_count"], "raw transcript count differs")
        for chrom, row, old in zip(CHROMOSOMES, report["chromosomes"], frozen["chromosomes"],
                                   strict=True):
            length = sizes[chrom]
            require(type(row["length"]) is int and length == row["length"] == old["length"],
                    f"{chrom}: chromosome length differs")
            whole, layers = chromosome_counts(reference, chrom, length, geometry.pop(chrom),
                                              split_soft_mask=True)
            check_counts(chrom, row["classes"], old["classes"], whole, layers)
            for target, counts in zip(totals, [whole, *layers], strict=True):
                for column, count in enumerate(counts):
                    target[column] += count

    whole, *layers = totals
    check_counts("all primary chromosomes", report["total"], frozen["total"], whole, layers)
    require(sum(whole) == coverage["totals"]["length"], "primary base count differs")
    mask_totals = {}
    for index, stratum in enumerate(STRATA):
        start = index * len(BASES)
        counts = dict(zip(NUCLEOTIDES, whole[start:start + len(BASES)], strict=True))
        mask_totals[stratum] = {"called_bases": sum(counts[b] for b in "acgt"),
                              "n_bases": counts["n"]}
    cells = len(FROZEN_LENGTHS) * len(BASES)
    return {
        "schema_version": 1,
        "reference_sha256": reference_hash,
        "genepred_sha256": raw_hash,
        "atlas_manifest_sha256": atlas["manifest_sha256"],
        "catalog_sha256": atlas["catalog_sha256"],
        "unstratified_summary_sha256": hashlib.sha256(frozen_bytes).hexdigest(),
        "retained_raw_models": retained,
        "verified_chromosomes": len(CHROMOSOMES),
        "verified_bases": sum(whole),
        "verified_cells": {
            "chromosome_stratified": len(CHROMOSOMES) * len(STRATA) * cells,
            "aggregate_stratified": len(STRATA) * cells,
            "chromosome_recombined": len(CHROMOSOMES) * cells,
            "aggregate_recombined": cells,
        },
        "mask_totals": mask_totals,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.reference, args.report), indent=2, sort_keys=True))
