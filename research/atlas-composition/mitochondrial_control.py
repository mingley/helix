#!/usr/bin/env python3
"""Check the curated chrM partition against GenBank features and retained UCSC DNA."""
import argparse
from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
LENGTH = 16_569
ACCESSION = "NC_012920.1"
GENEPRED = "datasets/genome-atlas/ncbiRefSeqCurated.txt.gz"
GENBANK = "datasets/genome-atlas/mitochondrial-reference.xml"
UCSC_DNA = "research/context-reference/ucsc-api-chrM-0-16569.json"
SOURCES = {
    GENEPRED: (4_801_131, "67090326b606658bdabdfe5f8dcb97b7943ac802c7945399522be7174ec021ca"),
    GENBANK: (176_955, "a97065e862ae56799ca419378fea8acb4453c3cf3733bf88830f06573f7fe6f0"),
    UCSC_DNA: (16_716, "a1a66813ed710f745e032ba2ac4df685b50abdfee6e000402f12120bc26c2dec"),
}


def frozen_bytes(name):
    size, digest = SOURCES[name]
    with (ROOT / name).open("rb") as stream:
        data = stream.read(size + 1)
    if len(data) != size or hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"source differs from frozen bytes: {name}")
    return data


def positions(start, end):
    if not 0 <= start < end <= LENGTH:
        raise ValueError(f"invalid chrM interval [{start}, {end})")
    return range(start, end)


def curated_regions():
    spans, exons, coding = set(), set(), set()
    models = 0
    with gzip.GzipFile(fileobj=io.BytesIO(frozen_bytes(GENEPRED))) as table:
        for line in table:
            fields = line.rstrip(b"\r\n").split(b"\t")
            if fields[2] != b"chrM":
                continue
            start, end, cds_start, cds_end, exon_count = map(int, fields[4:9])
            starts = [int(value) for value in fields[9].rstrip(b",").split(b",")]
            ends = [int(value) for value in fields[10].rstrip(b",").split(b",")]
            if (len(fields) != 16 or len(starts) != exon_count or len(ends) != exon_count
                    or not start <= cds_start <= cds_end <= end):
                raise ValueError("inconsistent chrM genePred model")
            models += 1
            spans.update(positions(start, end))
            for left, right in zip(starts, ends, strict=True):
                if not start <= left < right <= end:
                    raise ValueError("exon exceeds transcript span")
                exons.update(positions(left, right))
                begin, stop = max(left, cds_start), min(right, cds_end)
                if begin < stop:
                    coding.update(positions(begin, stop))
    if models != 13 or not coding <= exons <= spans:
        raise ValueError("unexpected curated chrM models")
    return spans, exons, coding, models


def genbank_regions():
    record = ET.fromstring(frozen_bytes(GENBANK)).find("GBSeq")
    if record is None or record.findtext("GBSeq_accession-version") != ACCESSION:
        raise ValueError("unexpected GenBank accession")
    sequence = record.findtext("GBSeq_sequence", "").upper()
    regions = {kind: set() for kind in ("CDS", "rRNA", "tRNA")}
    counts = Counter()
    for feature in record.findall("GBSeq_feature-table/GBFeature"):
        kind = feature.findtext("GBFeature_key")
        if kind not in regions:
            continue
        counts[kind] += 1
        for interval in feature.findall("GBFeature_intervals/GBInterval"):
            if interval.findtext("GBInterval_accession") != ACCESSION:
                raise ValueError("GenBank feature refers to another accession")
            first = int(interval.findtext("GBInterval_from"))
            last = int(interval.findtext("GBInterval_to"))
            if interval.find("GBInterval_iscomp[@value='true']") is not None:
                first, last = last, first
            # GenBank range endpoints are 1-based, inclusive and strand-ordered.
            regions[kind].update(positions(first - 1, last))
    if counts != {"CDS": 13, "rRNA": 2, "tRNA": 22}:
        raise ValueError("unexpected GenBank feature counts")
    return sequence, regions, counts


def control():
    spans, exons, cds, models = curated_regions()
    seq, genbank, feature_counts = genbank_regions()
    api = json.loads(frozen_bytes(UCSC_DNA))
    if (api["genome"], api["chrom"], api["start"], api["end"]) != ("hg38", "chrM", 0, LENGTH):
        raise ValueError("UCSC DNA has different coordinates")
    if len(seq) != LENGTH or set(seq) - set("ACGTN") or seq != api["dna"].upper():
        raise ValueError("GenBank sequence differs from retained UCSC chrM")
    if cds != genbank["CDS"] or spans != cds or exons != cds:
        raise ValueError("GenBank CDS and UCSC coding, exonic and transcript coverage disagree")

    outside = set(range(LENGTH)) - spans
    regions = {
        "cds": cds,
        "exon_outside_cds": exons - cds,
        "transcript_span_outside_exon": spans - exons,
        "outside_transcript_span": outside,
    }
    classes = {}
    for name, covered in regions.items():
        bases = Counter(seq[index].lower() for index in covered)
        classes[name] = {base: bases[base] for base in "acgtn"}
    return {
        "schema_version": 1,
        "accession": ACCESSION,
        "chromosome": "chrM",
        "length": LENGTH,
        "source_sha256": {name: digest for name, (_, digest) in SOURCES.items()},
        "curated_models": models,
        "genbank_features": dict(feature_counts),
        "class_widths": {name: len(covered) for name, covered in regions.items()},
        "classes": classes,
        "outside_transcript_span_with_genbank_rna": {
            "rRNA": len(outside & genbank["rRNA"]),
            "tRNA": len(outside & genbank["tRNA"]),
            "union": len(outside & (genbank["rRNA"] | genbank["tRNA"])),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, help="check the chrM row of a Rust atlas-composition report")
    args = parser.parse_args()
    result = control()
    if args.report is not None:
        report = json.loads(args.report.read_text(encoding="utf-8"))
        if (report["schema_version"], report["operation"]) != (1, "atlas_composition"):
            raise ValueError("unexpected Rust report schema or operation")
        rows = [row for row in report["chromosomes"] if row["name"] == "chrM"]
        if len(rows) != 1 or rows[0]["length"] != LENGTH or rows[0]["classes"] != result["classes"]:
            raise ValueError("Rust chrM row differs from the independent control")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
