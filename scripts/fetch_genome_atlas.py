#!/usr/bin/env python3
"""Freeze UCSC hg38 curated RefSeq annotation and build an offline context atlas."""
import argparse
import collections
import datetime
import gzip
import hashlib
import io
import json
import pathlib
import tempfile
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET

BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/"
MITO_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=gb&retmode=xml"
SOURCES = ["mitochondrial-reference.xml", "ncbiRefSeqCurated.txt.gz", "ncbiRefSeqCurated.sql",
           "ncbiRefSeqLink.txt.gz", "ncbiRefSeqLink.sql", "chromInfo.txt.gz", "chromInfo.sql"]
RNA_BIOTYPES = {"RNase_MRP_RNA", "RNase_P_RNA", "Y_RNA", "antisense_RNA", "lncRNA", "miRNA", "misc_RNA", "ncRNA", "rRNA", "scRNA", "snRNA", "snoRNA", "telomerase_RNA", "vault_RNA", "tRNA", "SRP_RNA"}
PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def fetch(name):
    url = MITO_URL if name == "mitochondrial-reference.xml" else BASE + name
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                data = response.read(35_000_001)
                if len(data) > 35_000_000:
                    raise ValueError("source exceeds bounded download budget")
                return data, dict(url=url, last_modified=response.headers.get("Last-Modified"),
                                  bytes=len(data), sha256=sha256(data))
        except (TimeoutError, urllib.error.URLError):
            if attempt == 2:
                raise
            time.sleep(1 + attempt)
    raise RuntimeError("unreachable")


def rows(data):
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as stream:
        for index, line in enumerate(stream):
            if index > 1_000_000:
                raise ValueError("source row budget exceeded")
            yield line.decode("utf-8").rstrip("\n").rstrip("\r").split("\t")


def merged(intervals):
    """Union of zero-based half-open intervals, merging adjacent ranges."""
    output = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError("invalid half-open interval")
        if output and start <= output[-1][1]:
            output[-1][1] = max(output[-1][1], end)
        else:
            output.append([start, end])
    return output


def covered(intervals):
    return sum(end - start for start, end in merged(intervals))


def normalize(raw):
    chromosome_sizes = {}
    for row in rows(raw["chromInfo.txt.gz"]):
        if row[0] in PRIMARY:
            chromosome_sizes[row[0]] = int(row[1])
    if set(chromosome_sizes) != set(PRIMARY):
        raise ValueError("missing primary chromosome lengths")
    links = {}
    protein_links = {}
    for row in rows(raw["ncbiRefSeqLink.txt.gz"]):
        if len(row) != 19:
            raise ValueError(f"unsupported RefSeqLink schema: {len(row)} columns")
        if row[0] in links:
            raise ValueError("duplicate RefSeq link accession")
        links[row[0]] = dict(gene_id=row[6], symbol=row[2], biotype=row[13])
        if row[5]:
            if row[5] in protein_links and protein_links[row[5]] != links[row[0]]:
                raise ValueError("ambiguous protein-to-gene mapping")
            protein_links[row[5]] = links[row[0]]
    mito = ET.fromstring(raw["mitochondrial-reference.xml"])
    if mito.findtext(".//GBSeq_accession-version") != "NC_012920.1":
        raise ValueError("mitochondrial accession mismatch")
    for feature in mito.findall(".//GBFeature"):
        if feature.findtext("GBFeature_key") != "CDS":
            continue
        qualifiers = collections.defaultdict(list)
        for qualifier in feature.findall("GBFeature_quals/GBQualifier"):
            qualifiers[qualifier.findtext("GBQualifier_name")].append(qualifier.findtext("GBQualifier_value"))
        protein = qualifiers["protein_id"][0]
        gene_ids = [value.removeprefix("GeneID:") for value in qualifiers["db_xref"] if value.startswith("GeneID:")]
        if len(gene_ids) != 1 or len(qualifiers["gene"]) != 1:
            raise ValueError("ambiguous mitochondrial gene metadata")
        metadata = dict(gene_id=gene_ids[0], symbol=qualifiers["gene"][0], biotype="protein_coding")
        if protein in protein_links and protein_links[protein] != metadata:
            raise ValueError("contradictory mitochondrial protein metadata")
        protein_links[protein] = metadata
    genes = {}
    excluded = collections.Counter()
    transcript_count = 0
    for row in rows(raw["ncbiRefSeqCurated.txt.gz"]):
        if len(row) != 16:
            raise ValueError(f"unsupported genePred schema: {len(row)} columns")
        accession, chrom, strand = row[1:4]
        if chrom not in chromosome_sizes:
            excluded[chrom] += 1
            continue
        link = links.get(accession, protein_links.get(accession))
        if link is None:
            raise ValueError(f"no gene metadata: {accession}")
        if not link["gene_id"].isdigit() or not link["symbol"] or not link["biotype"]:
            raise ValueError(f"missing gene identity or type: {accession}")
        start, end, cds_start, cds_end, exon_count = map(int, row[4:9])
        starts = [int(s) for s in row[9].strip(",").split(",")]
        ends = [int(s) for s in row[10].strip(",").split(",")]
        if strand not in ("+", "-") or len(starts) != exon_count or len(ends) != exon_count:
            raise ValueError(f"invalid transcript: {accession}")
        if not 0 <= start < end <= chromosome_sizes[chrom] or not start <= cds_start <= cds_end <= end:
            raise ValueError(f"out-of-bounds transcript: {accession}")
        exons = list(zip(starts, ends))
        if not all(start <= a < b <= end for a, b in exons) or any(exons[i][1] > exons[i+1][0] for i in range(len(exons)-1)):
            raise ValueError(f"invalid exon geometry: {accession}")
        key = (chrom, link["gene_id"], strand)
        if key not in genes:
            genes[key] = dict(gene_id=link["gene_id"], symbol=link["symbol"], biotype=link["biotype"],
                              chromosome=chrom, strand=strand, start=start, end=end,
                              transcript_count=0, transcript_intervals=[], exon_intervals=[], cds_intervals=[])
        gene = genes[key]
        if gene["symbol"] != link["symbol"] or gene["biotype"] != link["biotype"]:
            raise ValueError(f"inconsistent gene metadata: {key}")
        gene["start"] = min(gene["start"], start)
        gene["end"] = max(gene["end"], end)
        gene["transcript_count"] += 1
        gene["transcript_intervals"].append((start, end))
        gene["exon_intervals"].extend(exons)
        for a, b in exons:
            a, b = max(a, cds_start), min(b, cds_end)
            if a < b:
                gene["cds_intervals"].append((a, b))
        transcript_count += 1
    output = []
    for gene in genes.values():
        for field in ("transcript_intervals", "exon_intervals", "cds_intervals"):
            gene[field] = merged(gene[field])
        output.append(gene)
    rank = {name: i for i, name in enumerate(PRIMARY)}
    output.sort(key=lambda gene: (rank[gene["chromosome"]], gene["start"], gene["end"], gene["gene_id"], gene["strand"]))
    return dict(chromosomes=[dict(name=name, length=chromosome_sizes[name]) for name in PRIMARY],
                genes=output, retained_transcript_count=transcript_count,
                excluded_transcript_count=sum(excluded.values()), excluded_contig_count=len(excluded),
                excluded_contigs=dict(sorted(excluded.items())))


def summarize(data):
    genes = data["genes"]
    chromosomes = []
    for chrom in data["chromosomes"]:
        selected = [gene for gene in genes if gene["chromosome"] == chrom["name"]]
        chromosomes.append(dict(name=chrom["name"], length=chrom["length"], gene_loci=len(selected),
                                transcript_count=sum(gene["transcript_count"] for gene in selected),
                                transcript_span_union_bases=covered([i for gene in selected for i in gene["transcript_intervals"]]),
                                exon_union_bases=covered([i for gene in selected for i in gene["exon_intervals"]]),
                                cds_union_bases=covered([i for gene in selected for i in gene["cds_intervals"]])))
    hbb = [gene for gene in genes if gene["gene_id"] == "3043"]
    if len(hbb) != 1 or hbb[0]["symbol"] != "HBB" or hbb[0]["chromosome"] != "chr11" or hbb[0]["strand"] != "-":
        raise ValueError("HBB locus cannot be uniquely reconciled")
    hbb = hbb[0]
    def gap(gene):
        return max(hbb["start"] - gene["end"], gene["start"] - hbb["end"], 0)
    neighbors = sorted([gene for gene in genes if gene["chromosome"] == "chr11" and gene["gene_id"] != "3043"],
                       key=lambda gene: (gap(gene), gene["start"], gene["gene_id"]))[:12]
    neighbor_rows = [dict(gene_id=gene["gene_id"], symbol=gene["symbol"], biotype=gene["biotype"],
                          start=gene["start"], end=gene["end"], strand=gene["strand"], gap_bases=gap(gene),
                          genomic_side="overlapping" if gap(gene) == 0 else "lower_coordinate" if gene["end"] <= hbb["start"] else "higher_coordinate")
                     for gene in neighbors]
    totals = {key: sum(chrom[key] for chrom in chromosomes) for key in ("length", "transcript_span_union_bases", "exon_union_bases", "cds_union_bases")}
    return dict(schema_version=1, assembly="GRCh38/hg38", annotation="UCSC ncbiRefSeqCurated with ncbiRefSeqLink, frozen snapshot",
                scope="Primary chromosomes chr1-22, chrX, chrY, chrM only; curated RefSeq annotation is not a complete gene repertoire or functional map.",
                gene_loci=len(genes), distinct_gene_ids=len({gene["gene_id"] for gene in genes}),
                protein_coding_gene_loci=sum(gene["biotype"] == "protein_coding" for gene in genes),
                other_biotype_gene_loci=sum(gene["biotype"] != "protein_coding" for gene in genes),
                noncoding_rna_gene_loci=sum(gene["biotype"] in RNA_BIOTYPES for gene in genes),
                pseudogene_loci=sum("pseudogene" in gene["biotype"] for gene in genes),
                unclassified_biotype_gene_loci=sum(gene["biotype"] != "protein_coding" and gene["biotype"] not in RNA_BIOTYPES and "pseudogene" not in gene["biotype"] for gene in genes),
                gene_biotype_counts=dict(sorted(collections.Counter(gene["biotype"] for gene in genes).items())),
                retained_transcript_count=data["retained_transcript_count"],
                excluded_transcript_count=data["excluded_transcript_count"], excluded_contig_count=data["excluded_contig_count"],
                totals=totals, chromosomes=chromosomes, hbb=hbb, hbb_neighbors=neighbor_rows,
                interpretation=["Intervals are zero-based, half-open; coverage is the length of merged unions, never a sum over overlapping isoforms.",
                                "Gene loci group NCBI Gene ID, chromosome and strand; genes in X/Y pseudoautosomal regions can occupy more than one locus.",
                                "Transcript spans include introns; exon coverage includes coding and untranslated exons; CDS coverage is coding parts of annotated exons.",
                                "Other biotypes include noncoding RNAs and pseudogenes; the exact biotype counts are reported separately.",
                                "Sequence outside these annotations is unannotated by this dataset, not demonstrated nonfunctional.",
                                "Physical proximity and gene annotation alone do not establish regulatory interactions or disease mechanisms.",
                                "The source table includes 13 mitochondrial coding models, represented by protein accessions; their GeneIDs come from NC_012920.1 CDS qualifiers. It omits mitochondrial RNA genes."])


def catalog_text(data):
    metadata = {key: value for key, value in data.items() if key != "genes"}
    prefix = json.dumps(metadata, separators=(",", ":"))[:-1]
    records = ",\n".join(json.dumps(gene, separators=(",", ":")) for gene in data["genes"])
    return prefix + ',"genes":[\n' + records + '\n]}\n'


def build(output):
    output = pathlib.Path(output)
    if output.exists():
        raise ValueError("output exists; preserve snapshots and choose a new directory")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="genome-atlas-", dir=output.parent) as temp:
        temp = pathlib.Path(temp)
        raw, metadata = {}, []
        for name in SOURCES:
            payload, info = fetch(name)
            raw[name] = payload
            info["file"] = name
            metadata.append(info)
            (temp / name).write_bytes(payload)
        data = normalize(raw)
        report = summarize(data)
        (temp / "genes.json").write_text(catalog_text(data))
        (temp / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
        manifest = dict(schema_version=1, assembly="GRCh38/hg38", annotation="ncbiRefSeqCurated", coordinate_system="zero-based half-open",
                        retrieved_at=datetime.datetime.now(datetime.timezone.utc).isoformat(), sources=metadata,
                        primary_chromosomes=PRIMARY, files={path.name: sha256(path.read_bytes()) for path in sorted(temp.iterdir())})
        (temp / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        temp.rename(output)
    print(json.dumps({key: report[key] for key in ("gene_loci", "distinct_gene_ids", "protein_coding_gene_loci", "other_biotype_gene_loci", "retained_transcript_count", "totals")}, indent=2))


def verify(directory):
    root = pathlib.Path(directory)
    manifest = json.loads((root / "manifest.json").read_text())
    for name, expected in manifest["files"].items():
        if pathlib.Path(name).name != name:
            raise ValueError("unsupported manifest path")
        path = (root / name).resolve()
        if not path.is_relative_to(root.resolve()):
            raise ValueError("source path escapes snapshot")
        if sha256(path.read_bytes()) != expected:
            raise ValueError(f"checksum mismatch: {name}")
    raw = {name: (root / name).read_bytes() for name in SOURCES}
    data = normalize(raw)
    if data != json.loads((root / "genes.json").read_text()):
        raise ValueError("normalized catalog differs from raw sources")
    if summarize(data) != json.loads((root / "summary.json").read_text()):
        raise ValueError("summary differs from recomputed catalog")
    print("Atlas source hashes, raw-source normalization and summary verified offline.")


def self_test():
    assert merged([(4, 7), (1, 5), (7, 9), (20, 21), (1, 3)]) == [[1, 9], [20, 21]]
    assert covered([(0, 10), (2, 5), (9, 12)]) == 12
    assert covered([]) == 0
    try:
        merged([(10, 10)])
    except ValueError:
        pass
    else:
        raise AssertionError("invalid interval was accepted")
    print("Atlas importer interval tests passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--verify", metavar="DIRECTORY", help="Recompute a frozen atlas from raw sources offline")
    args = parser.parse_args()
    if sum(bool(value) for value in (args.output, args.self_test, args.verify)) != 1:
        parser.error("choose exactly one output directory, --self-test, or --verify")
    if args.self_test:
        self_test()
    elif args.verify:
        verify(args.verify)
    elif args.output:
        build(args.output)
    else:
        parser.error("provide an output directory or --self-test")
