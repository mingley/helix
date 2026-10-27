#!/usr/bin/env python3
"""Freeze the CC0 MaveDB HBB promoter MPRA with checked GRCh38 mapping."""
import argparse
import csv
import datetime
import hashlib
import io
import json
import math
import pathlib
import re
import urllib.request
import xml.etree.ElementTree as ET

URN = "urn:mavedb:00000018-a-1"
BASE = "https://api.mavedb.org/api/v1/score-sets/" + URN
SCOPE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6827219/fullTextXML"
MAX_BYTES = 2 * 1024 * 1024
EXPECTED_COLUMNS = [
    "accession", "hgvs_nt", "hgvs_splice", "hgvs_pro", "score", "p-value",
    "unique_tags", "dna_sequences", "rna_sequences"
]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url):
    with urllib.request.urlopen(url, timeout=90) as response:
        data = response.read(MAX_BYTES + 1)
        content_type = response.headers.get("Content-Type", "")
    if len(data) > MAX_BYTES:
        raise ValueError("response exceeds bounded acquisition size")
    return data, content_type


def validate_rows(raw, target):
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8")))
    if reader.fieldnames != EXPECTED_COLUMNS:
        raise ValueError("unexpected score CSV schema; review a new acquisition version")
    rows = list(reader)
    identifiers, snvs = set(), set()
    summary = dict(total_rows=len(rows), snv_rows=0, equality_rows=0,
                   other_representation_rows=0, finite_score_rows=0,
                   missing_score_rows=0, invalid_score_rows=0)
    for row in rows:
        identifier = row["accession"]
        if not identifier.startswith(URN + "#") or identifier in identifiers:
            raise ValueError("duplicate or invalid variant accession")
        identifiers.add(identifier)
        change = re.fullmatch(r"n\.(\d+)([ACGT])>([ACGT])", row["hgvs_nt"])
        equality = re.fullmatch(r"n\.(\d+)=", row["hgvs_nt"])
        if change:
            position, ref, alt = change.groups()
            position = int(position)
            if not 1 <= position <= len(target) or target[position - 1] != ref or ref == alt:
                raise ValueError("SNV reference/coordinate mismatch")
            key = (position, ref, alt)
            if key in snvs:
                raise ValueError("duplicate SNV representation")
            snvs.add(key)
            summary["snv_rows"] += 1
        elif equality:
            position = int(equality.group(1))
            if not 1 <= position <= len(target):
                raise ValueError("equality record outside source target")
            # Preserve literal submitted identity notation. It is not proof of
            # a neutral control or permission to infer a different mutation.
            summary["equality_rows"] += 1
        else:
            summary["other_representation_rows"] += 1
        score = row["score"].strip()
        if score in ("", "NA", "NaN", "nan", "null"):
            summary["missing_score_rows"] += 1
        else:
            try:
                if not math.isfinite(float(score)):
                    summary["invalid_score_rows"] += 1
                else:
                    summary["finite_score_rows"] += 1
            except ValueError:
                summary["invalid_score_rows"] += 1
    return summary


def load_reference(root):
    root = pathlib.Path(root).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    if manifest["assembly"] != "GRCh38" or manifest["gene"] != "HBB":
        raise ValueError("reference must be the checked GRCh38 HBB bundle")
    raw = (root / "genomic.fasta").read_bytes()
    if digest(raw) != manifest["files"]["genomic.fasta"]:
        raise ValueError("reference FASTA checksum mismatch")
    sequence = "".join(raw.decode().splitlines()[1:]).upper()
    if len(sequence) != manifest["sequence_end"] - manifest["sequence_start"]:
        raise ValueError("reference interval length mismatch")
    return manifest, sequence, digest(raw)


def validate_license(metadata):
    license_info = metadata["license"]
    if (metadata["urn"] != URN or license_info["shortName"] != "CC0"
            or license_info["link"] != "https://creativecommons.org/publicdomain/zero/1.0/"
            or license_info["version"] != "1.0"):
        raise ValueError("the score-set data license is not the reviewed CC0 1.0 grant")
    return license_info


def map_target(metadata, reference):
    target = metadata["targetGenes"][0]["targetSequence"]["sequence"].upper()
    taxonomy = metadata["targetGenes"][0]["targetSequence"]["taxonomy"]
    if len(target) != 187 or set(target) - set("ACGT") or taxonomy["code"] != 9606:
        raise ValueError("unexpected target length, residues or species")
    reference_manifest, genomic, reference_sha = load_reference(reference)
    reverse = target.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    matches = [(i, "+") for i in range(len(genomic)) if genomic.startswith(target, i)]
    matches += [(i, "-") for i in range(len(genomic)) if genomic.startswith(reverse, i)]
    if len(matches) != 1:
        raise ValueError("source target must map uniquely in one orientation")
    offset, strand = matches[0]
    start = reference_manifest["sequence_start"] + offset
    end = start + len(target)
    source_coords = metadata["extraMetadata"]
    if (source_coords != {"chr": "11", "end": 5227208, "start": 5227022, "reference": "hg38"}
            or (start, end, strand) != (5227021, 5227208, "+")):
        raise ValueError("target mapping differs from reviewed one-based-inclusive source coordinates")
    return target, start, end, strand, reference_sha


def freeze(output, reference):
    output = pathlib.Path(output)
    if output.exists():
        raise ValueError("output already exists; choose a new directory")
    metadata_raw, metadata_type = fetch(BASE)
    metadata = json.loads(metadata_raw)
    license_info = validate_license(metadata)
    policy_raw, _ = fetch(SCOPE_URL)
    policy = ET.fromstring(policy_raw)
    paragraphs = [
        " ".join("".join(p.itertext()).split()) for p in policy.findall(".//p")
    ]
    scope = next((p for p in paragraphs if "licenses for each score set" in p
                  and "downloaded file" in p), None)
    if scope is None:
        raise ValueError("score-set license scope could not be verified")
    target, start, end, strand, reference_sha = map_target(metadata, reference)
    # Retrieve scores only after checking the license and exact reference target.
    scores_raw, scores_type = fetch(BASE + "/scores")
    if "text/csv" not in scores_type:
        raise ValueError("score endpoint did not return CSV")
    summary = validate_rows(scores_raw, target)
    if summary["total_rows"] != metadata["numVariants"] or summary["total_rows"] != 623:
        raise ValueError("variant count differs from reviewed metadata")
    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    scope_document = dict(
        data_license=license_info,
        scope_source_url=SCOPE_URL,
        scope_source_sha256=digest(policy_raw),
        scope_source_doi="10.1186/s13059-019-1845-6",
        scope_excerpt=scope,
        interpretation="MaveDB's per-score-set Creative Commons license applies to the deposited/downloaded data, not merely catalog metadata.",
        producer_policy=metadata.get("dataUsagePolicy"),
        producer_publication=dict(doi="10.1038/s41467-019-11526-w", pmid="31395865",
                                  article_license="CC BY 4.0; distinct from this score-set CC0 grant"),
    )
    target_fasta = (f">{URN}|HBB_promoter|GRCh38:11:{start}-{end}|strand={strand}\n"
                    + "\n".join(target[i:i+80] for i in range(0, len(target), 80)) + "\n").encode()
    files = {
        "metadata.json": metadata_raw, "scores.csv": scores_raw,
        "target.fasta": target_fasta,
        "license-scope.json": (json.dumps(scope_document, indent=2) + "\n").encode(),
    }
    manifest = dict(
        schema_version=1, dataset_id="hbb-promoter-mpra", assay_urn=URN,
        source_url=BASE, score_url=BASE + "/scores", retrieved_at=now,
        assembly="GRCh38", chromosome="11",
        target_start=start, target_end=end, target_strand=strand,
        gene_transcription_strand="-", expected_rows=623,
        cell_line="HEL 92.1.7",
        score_name="score", score_units="log2 reporter-expression effect; linear-model coefficient",
        replicates=dict(cultures=3, library="Same constructed library across transfections",
                        independent_studies=1, separate_replicate_scores_in_this_csv=False),
        source_coordinate_convention="One-based inclusive; verified by exact sequence mapping",
        reference_fasta_sha256=reference_sha,
        source_target_sequence_sha256=digest(target.encode()),
        reference_5prime_utr_interval=[5227021, 5227071],
        upstream_interval=[5227071, 5227208],
        data_license=license_info,
        validation=summary,
        sources=[
            dict(file="metadata.json", url=BASE, content_type=metadata_type,
                 retrieved_at=now, sha256=digest(metadata_raw)),
            dict(file="scores.csv", url=BASE + "/scores", content_type=scores_type,
                 retrieved_at=now, sha256=digest(scores_raw)),
            dict(file="license-scope.json", url=SCOPE_URL, retrieved_at=now,
                 source_document_sha256=digest(policy_raw)),
        ],
        files={name: digest(raw) for name, raw in files.items()},
        analysis_scope="Acquisition, reference mapping and explicit missingness/representation accounting only; no fitted model, thresholds or outcome-pattern analysis",
    )
    output.mkdir(parents=True)
    for name, raw in files.items():
        (output / name).write_bytes(raw)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"mapping": [start, end, strand], "validation": summary}))


def verify(directory, reference):
    """Verify the retained snapshot without network requests or writes."""
    root = pathlib.Path(directory).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    if (manifest["schema_version"] != 1 or manifest["assay_urn"] != URN
            or manifest["assembly"] != "GRCh38" or manifest["chromosome"] != "11"
            or manifest["gene_transcription_strand"] != "-"):
        raise ValueError("unsupported assay manifest identity")
    required = {"metadata.json", "scores.csv", "target.fasta", "license-scope.json"}
    if set(manifest["files"]) != required:
        raise ValueError("manifest must pin exactly the reviewed input files")
    contents = {}
    for name, expected in manifest["files"].items():
        path = (root / name).resolve()
        if path.parent != root:
            raise ValueError("snapshot input escapes its directory")
        raw = path.read_bytes()
        if digest(raw) != expected:
            raise ValueError(f"checksum mismatch: {name}")
        contents[name] = raw
    metadata = json.loads(contents["metadata.json"])
    license_info = validate_license(metadata)
    scope = json.loads(contents["license-scope.json"])
    if (manifest["data_license"] != license_info or scope["data_license"] != license_info
            or scope["scope_source_url"] != SCOPE_URL
            or scope["scope_source_doi"] != "10.1186/s13059-019-1845-6"
            or "licenses for each score set" not in scope["scope_excerpt"]
            or "downloaded file" not in scope["scope_excerpt"]):
        raise ValueError("score-set license or its documented scope is inconsistent")
    target, start, end, strand, reference_sha = map_target(metadata, reference)
    if (manifest["target_start"], manifest["target_end"], manifest["target_strand"]) != (start, end, strand):
        raise ValueError("manifest target mapping is inconsistent")
    if (manifest["reference_fasta_sha256"] != reference_sha
            or manifest["source_target_sequence_sha256"] != digest(target.encode())):
        raise ValueError("reference or target commitment mismatch")
    fasta = contents["target.fasta"].decode().splitlines()
    if not fasta or not fasta[0].startswith(">") or "".join(fasta[1:]).upper() != target:
        raise ValueError("retained target FASTA differs from source")
    summary = validate_rows(contents["scores.csv"], target)
    if (summary != manifest["validation"] or summary["total_rows"] != manifest["expected_rows"]
            or summary["total_rows"] != metadata["numVariants"] or summary["total_rows"] != 623):
        raise ValueError("row accounting differs from source/manifest")
    print(json.dumps({"verified_offline": True, "mapping": [start, end, strand],
                      "validation": summary}))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?")
    parser.add_argument("--verify", metavar="DIRECTORY",
                        help="Check a frozen snapshot offline without writing files")
    parser.add_argument("--reference", default="datasets/hbb-reference")
    args = parser.parse_args()
    if bool(args.output) == bool(args.verify):
        parser.error("choose exactly one output directory or --verify DIRECTORY")
    if args.verify:
        verify(args.verify, args.reference)
    else:
        freeze(args.output, args.reference)
