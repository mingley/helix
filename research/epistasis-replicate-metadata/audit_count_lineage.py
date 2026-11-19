#!/usr/bin/env python3
"""Check which input sequencing libraries the OL27 count tables share."""
import hashlib
import json
import urllib.request

SOURCE_REVISION = "fc2d3daabfcf90655147ba2f9d11680f3a8d1045"
SOURCE_BASE = (
    "https://raw.githubusercontent.com/ENCODE-DCC/geo-submission/"
    f"{SOURCE_REVISION}/experiments"
)
BARCODE_MAP = "ENCFF823YWM"
PLASMID = ("ENCSR915INT", "c5e4ef009653b55afb1bd26788ff4791477b04a4")
CONTEXTS = (
    ("A549", "ENCSR729IRS", "8469c5903b3219e94576878d0edbad97ac55c169", 1, "ENCFF018AMJ"),
    ("GM12878", "ENCSR645CIP", "04b4ad1226d5d32aca5eab004afb5b3d708834a5", 2, "ENCFF345ASG"),
    ("HepG2", "ENCSR847SNQ", "7feaa964c0a9335f5ce9c43f310b900079ffa3ae", 2, "ENCFF318XMJ"),
    ("K562", "ENCSR857JFG", "b36de81bf6d6e6cefaefad80b0d2dbdd720c858e", 2, "ENCFF970OLE"),
    ("SK-N-SH", "ENCSR623DJZ", "3c013ebec2501eafab9d9bf0b5d40332f70446fe", 2, "ENCFF996ECA"),
)
MAX_RESPONSE = 1024 * 1024


def read_experiment(accession, expected_blob):
    url = f"{SOURCE_BASE}/{accession}_modified.json"
    request = urllib.request.Request(url, headers={"User-Agent": "Helix-source-audit/1.0"})
    with urllib.request.urlopen(request, timeout=20) as response:
        content = response.read(MAX_RESPONSE + 1)
    if len(content) > MAX_RESPONSE:
        raise ValueError(f"{accession}: experiment JSON exceeds the read limit")
    git_blob = f"blob {len(content)}\0".encode() + content
    if hashlib.sha1(git_blob).hexdigest() != expected_blob:
        raise ValueError(f"{accession}: experiment JSON differs from the pinned Git blob")
    record = json.loads(content)
    if record["accession"] != accession:
        raise ValueError(f"{accession}: wrong experiment in source")
    return record, {
        "experiment": accession,
        "source_url": url,
        "source_git_blob_sha1": expected_blob,
        "run_type": "single-ended",
    }


def library_reads(record, analyte, files_per_library):
    accession = record["accession"]
    rows = {}
    for rep in record["replicates"]:
        number = rep["biological_replicate_number"]
        library = rep["library"]
        if (number in rows or rep["technical_replicate_number"] != 1
                or library["nucleic_acid_term_name"] != analyte):
            raise ValueError(f"{accession}: unexpected library for replicate {number}")
        rows[number] = {
            "biological_replicate": number,
            "technical_replicate": 1,
            "library": library["accession"],
            "fastq_accessions": [],
        }
    if set(rows) != set(range(1, 6)):
        raise ValueError(f"{accession}: expected biological replicate labels 1 through 5")

    reads = set()
    for file in record["files"]:
        if file["output_type"] != "reads":
            continue
        identity = file["accession"]
        rep = file["replicate"]
        number = rep["biological_replicate_number"]
        if ((file["file_format"], file["run_type"], file["status"])
                != ("fastq", "single-ended", "released")
                or rep["technical_replicate_number"] != 1 or number not in rows
                or identity in reads):
            raise ValueError(f"{accession}: unexpected FASTQ {identity}")
        reads.add(identity)
        rows[number]["fastq_accessions"].append(identity)

    libraries = [rows[index] for index in sorted(rows)]
    for row in libraries:
        row["fastq_accessions"].sort()
        if len(row["fastq_accessions"]) != files_per_library:
            raise ValueError(f"{accession}: unexpected number of FASTQs per library")
    return libraries, reads


def count_lineage():
    record, plasmid = read_experiment(*PLASMID)
    if (record["assay_term_name"], record["biosample_term_name"]) != (
            "pooled clone sequencing", "DNA cloning sample"):
        raise ValueError("expected the plasmid cloning experiment")
    plasmid["replicate_libraries"], plasmid_reads = library_reads(record, "DNA", 2)
    contexts = {}
    all_reads = plasmid_reads.copy()

    for cell, accession, blob, files_per_library, output in CONTEXTS:
        record, lineage = read_experiment(accession, blob)
        if (record["assay_term_name"], record["biosample_term_name"]) != ("MPRA", cell):
            raise ValueError(f"{accession}: wrong MPRA context in source")
        lineage["replicate_libraries"], rna_reads = library_reads(record, "RNA", files_per_library)
        if all_reads & rna_reads:
            raise ValueError(f"{accession}: local read accessions overlap another experiment")
        all_reads.update(rna_reads)

        matching = [file for file in record["files"] if file["accession"] == output]
        if len(matching) != 1:
            raise ValueError(f"{output}: expected a unique count table in {accession}")
        table = matching[0]
        if (table["file_format"], table["status"], table["output_type"]) != (
                "tsv", "released", "element quantifications"):
            raise ValueError(f"{output}: unexpected count table metadata")
        expected = {f"/files/{key}/" for key in plasmid_reads | rna_reads | {BARCODE_MAP}}
        if sorted(table["derived_from"]) != sorted(expected):
            raise ValueError(f"{output}: upstream inputs differ from the source libraries")
        lineage["element_quantifications"] = output
        lineage["shared_plasmid_fastq_count"] = len(plasmid_reads)
        contexts[cell] = lineage

    library_ids = [row["library"] for record in [plasmid, *contexts.values()]
                   for row in record["replicate_libraries"]]
    if len(set(library_ids)) != len(library_ids):
        raise ValueError("library accessions overlap across replicate labels or experiments")
    return {
        "source_revision": SOURCE_REVISION,
        "barcode_map": BARCODE_MAP,
        "plasmid": plasmid,
        "contexts": contexts,
    }


if __name__ == "__main__":
    print(json.dumps(count_lineage(), indent=2, sort_keys=True))
