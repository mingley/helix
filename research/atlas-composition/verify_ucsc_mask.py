#!/usr/bin/env python3
"""Reconcile two masked hg38 reports with UCSC's published whole-image faSize."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

CONTEXT = Path(__file__).resolve().parents[1] / "context-reference"
KENT = "https://github.com/ucscGenomeBrowser/kent/blob/fa1ec8f7a7071301796d3050f7bb566679b55bbe/"
BUILD_LOG = KENT + "src/hg/makeDb/doc/hg38/hg38.txt#L335-L380"
FASIZE_SOURCE = KENT + "src/utils/faSize/faSize.c#L193-L213"
# faSize combines N/n, then counts the remaining upper- and lowercase bases.
UCSC_FASIZE = {
    "sequences": 455, "bases": 3_209_286_105, "n": 159_970_322,
    "unmasked_called": 1_460_684_798, "masked_called": 1_588_630_985,
}
STRATA = ("unmasked", "masked")
BASES = "acgtn"
CLASSES = ("cds", "exon_outside_cds", "transcript_span_outside_exon", "outside_transcript_span")
PRIMARY = tuple(f"chr{number}" for number in range(1, 23)) + ("chrX", "chrY", "chrM")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def bins(groups, label):
    require(set(groups) == set(STRATA), f"{label}: invalid mask strata")
    counts = Counter()
    for stratum, values in groups.items():
        require(set(values) == set(BASES) and all(type(n) is int and n >= 0 for n in values.values()),
                f"{label}/{stratum}: expected nonnegative integer A/C/G/T/N counts")
        counts.update({(stratum, base): count for base, count in values.items()})
    return counts


def scope(counts, sequences):
    return {
        "sequences": sequences,
        "bases": counts.total(),
        "n": counts["unmasked", "n"] + counts["masked", "n"],
        "unmasked_called": sum(counts["unmasked", base] for base in "acgt"),
        "masked_called": sum(counts["masked", base] for base in "acgt"),
    }


def verify(full_path, atlas_path):
    full = json.loads(full_path.read_bytes())
    atlas = json.loads(atlas_path.read_bytes())
    manifest = json.loads((CONTEXT / "manifest.json").read_bytes())
    catalog_data = (CONTEXT / "twobit-structure.json").read_bytes()
    catalog_sha = hashlib.sha256(catalog_data).hexdigest()
    require(catalog_sha == manifest["files"]["twobit-structure.json"], "twoBit catalog checksum differs")
    catalog = json.loads(catalog_data)
    for report, operation in ((full, "summary_by_soft_mask"), (atlas, "atlas_composition_by_soft_mask")):
        require(report["schema_version"] == 1 and report["operation"] == operation,
                f"expected schema 1 {operation}")
        ref = report["reference"]
        require(ref["sha256"] == manifest["reference_sha256"]
                and ref["file_bytes"] == manifest["reference_bytes"] == catalog["file_bytes"]
                and ref["format"] == catalog["format"] and ref["format_version"] == catalog["version"]
                and ref["byte_order"] == catalog["byte_order"]
                and ref["declared_assembly"] in (None, "GRCh38"),
                f"{operation}: reference identity differs from frozen hg38")

    contigs = {}
    whole_counts = Counter()
    for row in full["contigs"]:
        name = row["name"]
        require(name in catalog["sequences"] and name not in contigs, f"unexpected contig {name}")
        counts = bins(row["counts"], name)
        require(type(row["length"]) is int and row["length"] == catalog["sequences"][name]["length"]
                and counts.total() == row["length"], f"{name}: full summary length differs")
        contigs[name] = counts
        whole_counts.update(counts)
    require(contigs.keys() == catalog["sequences"].keys(), "full summary is missing contigs")
    require(whole_counts == bins(full["total"], "full total"), "full contig counts disagree with total")
    whole = scope(whole_counts, len(contigs))
    require(whole == UCSC_FASIZE, f"whole-image counts {whole} differ from UCSC faSize {UCSC_FASIZE}")

    require(tuple(row["name"] for row in atlas["chromosomes"]) == PRIMARY, "primary chromosome set or order differs")
    class_totals = {label: Counter() for label in CLASSES}
    primary_counts = Counter()
    for row in atlas["chromosomes"]:
        name = row["name"]
        require(type(row["length"]) is int and row["length"] == catalog["sequences"][name]["length"],
                f"{name}: primary length differs from catalog")
        require(set(row["classes"]) == set(CLASSES), f"{name}: unexpected annotation classes")
        across_classes = Counter()
        for label, values in row["classes"].items():
            counts = bins(values, f"{name}/{label}")
            across_classes.update(counts)
            class_totals[label].update(counts)
        require(across_classes == contigs[name], f"{name}: class sums differ from full summary")
        primary_counts.update(across_classes)
    require(set(atlas["total"]) == set(CLASSES), "unexpected classes in primary total")
    for label, expected in class_totals.items():
        require(expected == bins(atlas["total"][label], f"primary total/{label}"),
                f"{label}: primary total differs from its chromosomes")

    return {
        "schema_version": 1,
        "reference_sha256": manifest["reference_sha256"],
        "ucsc_published_md5": manifest["reference_published_md5"],
        "catalog_sha256": catalog_sha,
        "ucsc_build_log": BUILD_LOG,
        "ucsc_fasize_convention": FASIZE_SOURCE,
        "whole_image": whole,
        "whole_image_masked_n": whole_counts["masked", "n"],
        "primary_chromosomes": scope(primary_counts, len(PRIMARY)),
        "primary_masked_n": primary_counts["masked", "n"],
        "primary_mask_base_cells_reconciled": len(PRIMARY) * len(STRATA) * len(BASES),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("full_summary", type=Path)
    parser.add_argument("primary_atlas", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.full_summary, args.primary_atlas), indent=2, sort_keys=True))
