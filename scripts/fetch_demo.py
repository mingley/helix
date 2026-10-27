#!/usr/bin/env python3
"""Refresh a public, non-patient-level HBB demonstration fixture into a NEW directory."""
import argparse
import datetime
import hashlib
import json
import pathlib
import time
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
QUERY = 'HBB[gene] AND "single nucleotide variant"[Type] AND ("pathogenic"[Clinical Significance] OR "benign"[Clinical Significance])'


def fetch(endpoint, **params):
    url = BASE + endpoint + "?" + urllib.parse.urlencode(params)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=90) as response:
                return response.read(), url
        except (TimeoutError, urllib.error.URLError):
            if attempt == 2:
                raise
            time.sleep(1 + attempt)
    raise RuntimeError("unreachable")


def consequence(names):
    names = " ".join(names).lower()
    for text, value in [
        ("splice", "splice"), ("stop", "stop"), ("missense", "missense"),
        ("synonymous", "synonymous"), ("intron", "intron"), ("utr", "utr")
    ]:
        if text in names:
            return value
    return "other"


def build(output, cached_source=None):
    output = pathlib.Path(output)
    if output.exists():
        raise ValueError("output already exists; preserve the frozen fixture and choose a new directory")
    if cached_source:
        raw = pathlib.Path(cached_source).read_bytes()
        result = json.loads(raw)["result"]
        ids = result["uids"]
        url = BASE + "esummary.fcgi?" + urllib.parse.urlencode(
            dict(db="clinvar", id=",".join(ids), retmode="json"))
    else:
        raw, _ = fetch("esearch.fcgi", db="clinvar", term=QUERY, retmax=100, retmode="json")
        ids = json.loads(raw)["esearchresult"]["idlist"]
        time.sleep(0.4)
        raw, url = fetch("esummary.fcgi", db="clinvar", id=",".join(ids), retmode="json")
        result = json.loads(raw)["result"]
    rows = []
    seen = set()
    for uid in ids:
        record = result[uid]
        classification = record["germline_classification"]["description"]
        labels = {"Benign": "benign", "Likely benign": "benign",
                  "Pathogenic": "pathogenic", "Likely pathogenic": "pathogenic"}
        if classification not in labels or len(record["variation_set"]) != 1:
            continue
        spdi = record["variation_set"][0].get("canonical_spdi", "").split(":")
        if len(spdi) != 4:
            continue
        accession, position, ref, alt = spdi
        if accession != "NC_000011.10" or len(ref) != 1 or len(alt) != 1:
            continue
        if ref not in "ACGT" or alt not in "ACGT" or ref == alt:
            continue
        key = (int(position), ref, alt)
        if key in seen:
            continue
        seen.add(key)
        rows.append(dict(id=record["accession_version"], position=int(position), ref=ref, alt=alt,
                         consequence=consequence(record["molecular_consequence_list"]),
                         label=labels[classification], classification=classification,
                         review_status=record["germline_classification"]["review_status"],
                         source_url=f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{uid}/"))
    rows.sort(key=lambda row: (row["position"], row["ref"], row["alt"]))
    if len(rows) < 10:
        raise ValueError("insufficient records")
    start = min(row["position"] for row in rows)
    end = max(row["position"] for row in rows) + 1
    time.sleep(0.4)
    fasta, fasta_url = fetch("efetch.fcgi", db="nuccore", id="NC_000011.10",
                             seq_start=start + 1, seq_stop=end, rettype="fasta", retmode="text")
    sequence = "".join(fasta.decode().splitlines()[1:]).upper()
    if len(sequence) != end - start:
        raise ValueError("reference length mismatch")
    for row in rows:
        if sequence[row["position"] - start] != row["ref"]:
            raise ValueError(f"reference mismatch: {row['id']}")
    output.mkdir(parents=True)
    (output / "source.json").write_bytes(raw)
    (output / "reference.fasta").write_bytes(fasta)
    vcf = [
        "##fileformat=VCFv4.3", "##reference=GRCh38",
        "##contig=<ID=11,length=135086622>",
        '##INFO=<ID=MC,Number=1,Type=String,Description="ClinVar molecular consequence group">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ]
    partitions = {"training": [], "heldout": []}
    for row in rows:
        vcf.append(f"11\t{row['position'] + 1}\t{row['id']}\t{row['ref']}\t{row['alt']}\t.\tPASS\tMC={row['consequence']}")
        # All alleles at a position stay together. Assignment does not inspect labels.
        group = hashlib.sha256(f"NC_000011.10:{row['position']}".encode()).digest()[0] % 2
        partition = "training" if group == 0 else "heldout"
        partitions[partition].append({key: row[key] for key in
                                     ("id", "label", "classification", "review_status", "source_url")})
    (output / "variants.vcf").write_text("\n".join(vcf) + "\n")
    for name, labels in partitions.items():
        if {row["label"] for row in labels} != {"benign", "pathogenic"}:
            raise ValueError(f"{name} needs both classes")
        (output / f"{name}.json").write_text(json.dumps(labels, indent=2) + "\n")
    files = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
             for path in sorted(output.iterdir())}
    manifest = dict(schema_version=1, assembly="GRCh38", gene="HBB",
                    reference_accession="NC_000011.10", reference_start=start, reference_end=end,
                    source_url=url, reference_url=fasta_url, query=QUERY,
                    retrieved_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    split="sha256(NC_000011.10:<zero-based position>)[0] % 2; 0=training; 1=heldout",
                    scope="Illustrative ClinVar classification benchmark; disjoint loci, NOT independent experimental cohorts.",
                    files=files)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({name: len(labels) for name, labels in partitions.items()}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    parser.add_argument("--cached-source", help="Existing unmodified NCBI esummary response")
    args = parser.parse_args()
    build(args.output, args.cached_source)
