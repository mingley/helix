#!/usr/bin/env python3
"""Check chrM soft masking against case-preserved UCSC DNA and GenBank features."""
import argparse
from collections import Counter
import json
from pathlib import Path

from mitochondrial_control import ACCESSION, LENGTH, SOURCES, checked_partition, genbank_features


def intervals(covered):
    """Coalesce positions into sorted zero-based, half-open intervals."""
    runs = []
    for pos in sorted(covered):
        if runs and runs[-1][1] == pos:
            runs[-1][1] = pos + 1
        else:
            runs.append([pos, pos + 1])
    return runs


def reference_n_sites(dna, sequence, features):
    sites = []
    for i, base in enumerate(dna):
        if base.upper() != "N":
            continue
        annotations = []
        for kind, feature, covered in features:
            if kind == "misc_feature" and i in covered:
                notes = feature.findall("GBFeature_quals/GBQualifier[GBQualifier_name='note']/GBQualifier_value")
                annotations.append({
                    "genbank_location_1based": feature.findtext("GBFeature_location"),
                    "notes": [note.text for note in notes],
                })
        sites.append({
            "start0": i, "end0": i + 1, "ucsc_base": base, "genbank_base": sequence[i],
            "genbank_misc_features": annotations,
        })
    return sites


def control():
    sequence, features, _ = genbank_features(("D-loop", "misc_feature"))
    genbank_cds = set().union(*(covered for kind, _, covered in features if kind == "CDS"))
    dna, regions, _ = checked_partition(sequence, genbank_cds)
    masked = {i for i, base in enumerate(dna) if base.islower()}
    classes = {}
    for name, covered in regions.items():
        counts = Counter(dna[i] for i in covered)
        classes[name] = {
            "masked": {base: counts[base] for base in "acgtn"},
            "unmasked": {base: counts[base.upper()] for base in "acgtn"},
        }

    masked_by_kind = {kind: set() for kind in ("CDS", "rRNA", "tRNA")}
    feature_memberships = Counter()
    overlaps = []
    for kind, feature, covered in features:
        if kind not in masked_by_kind:
            continue
        hit = masked & covered
        if not hit:
            continue
        masked_by_kind[kind].update(hit)
        feature_memberships[kind] += len(hit)
        gene = feature.findtext("GBFeature_quals/GBQualifier[GBQualifier_name='gene']/GBQualifier_value")
        if not gene:
            raise ValueError(f"masked {kind} feature has no GenBank gene qualifier")
        overlaps.append({
            "gene": gene,
            "kind": kind,
            "genbank_location_1based": feature.findtext("GBFeature_location"),
            "masked_intervals0": intervals(hit),
        })
    masked_rna = masked_by_kind["rRNA"] | masked_by_kind["tRNA"]
    residual = masked - masked_rna - masked_by_kind["CDS"]
    in_d_loop = set()
    d_loop_hits = []
    for kind, feature, covered in features:
        if kind == "D-loop" and (hit := residual & covered):
            in_d_loop.update(hit)
            d_loop_hits.append({
                "genbank_location_1based": feature.findtext("GBFeature_location"),
                "masked_intervals0": intervals(hit),
            })
    return {
        "schema_version": 1,
        "accession": ACCESSION,
        "chromosome": "chrM",
        "length": LENGTH,
        "source_sha256": {name: digest for name, (_, digest) in SOURCES.items()},
        "masked_bases": len(masked),
        "masked_intervals0": intervals(masked),
        "classes": classes,
        "masked_genbank_features": overlaps,
        "masked_feature_memberships_by_kind": {kind: feature_memberships[kind] for kind in masked_by_kind},
        "masked_union_bases_by_kind": {kind: len(covered) for kind, covered in masked_by_kind.items()},
        "masked_genbank_rna_union_bases": len(masked_rna),
        "masked_outside_genbank_cds_and_rna": {
            "total_bases": len(residual),
            "in_d_loop_bases": len(in_d_loop),
            "d_loop_features": d_loop_hits,
            "outside_d_loop_intervals0": intervals(residual - in_d_loop),
        },
        "reference_n_sites": reference_n_sites(dna, sequence, features),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, help="check the chrM row of a Rust soft-mask atlas report")
    args = parser.parse_args()
    result = control()
    if args.report is not None:
        report = json.loads(args.report.read_bytes())
        if (report["schema_version"], report["operation"]) != (1, "atlas_composition_by_soft_mask"):
            raise ValueError("unexpected Rust report schema or operation")
        rows = [row for row in report["chromosomes"] if row["name"] == "chrM"]
        if len(rows) != 1 or rows[0]["length"] != LENGTH or rows[0]["classes"] != result["classes"]:
            raise ValueError("Rust chrM soft-mask counts differ from the retained UCSC sequence")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
