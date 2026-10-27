#!/usr/bin/env python3
"""Freeze and reconcile a public GRCh38 HBB reference bundle in a new directory.

Raw NCBI and UCSC responses are retained unchanged. All derived coordinates are
zero-based half-open. A successful build requires exact transcript and protein
agreement across the genomic annotation and versioned RefSeq sequences.
"""
import argparse
import datetime
import hashlib
import itertools
import json
import pathlib
import re
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request


ACCESSION = "NC_000011.10"
TRANSCRIPT = "NM_000518.5"
PROTEIN = "NP_000509.1"
SEQUENCE_START = 5_225_000
SEQUENCE_END = 5_227_500
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
SOURCES = {
    "annotation.json": (
        "https://api.genome.ucsc.edu/getData/track?"
        "genome=hg38;track=ncbiRefSeq;chrom=chr11;start=5225000;end=5227500"
    ),
    "genomic.fasta": EUTILS + urllib.parse.urlencode(
        dict(db="nuccore", id=ACCESSION, seq_start=SEQUENCE_START + 1,
             seq_stop=SEQUENCE_END, rettype="fasta", retmode="text")),
    "transcript.fasta": EUTILS + urllib.parse.urlencode(
        dict(db="nuccore", id=TRANSCRIPT, rettype="fasta", retmode="text")),
    "protein.fasta": EUTILS + urllib.parse.urlencode(
        dict(db="protein", id=PROTEIN, rettype="fasta", retmode="text")),
    "transcript.gb": EUTILS + urllib.parse.urlencode(
        dict(db="nuccore", id=TRANSCRIPT, rettype="gb", retmode="text")),
}


def sha256(content):
    return hashlib.sha256(content).hexdigest()


def download(url):
    """Use normal verified HTTPS and inherited proxy/CA configuration."""
    request = urllib.request.Request(
        url, headers={"User-Agent": "Helix-public-genomics-research/0.1"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                content = response.read()
            time.sleep(0.4)  # NCBI's unauthenticated E-utilities rate limit.
            return content
        except urllib.error.HTTPError as error:
            if error.code not in {429, 500, 502, 503, 504} or attempt == 2:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == 2:
                raise
        time.sleep(1 + attempt)
    raise RuntimeError("unreachable")


def fasta(content, expected_id):
    lines = content.decode("ascii").splitlines()
    if not lines or not lines[0].startswith(">" + expected_id):
        raise ValueError(f"FASTA must identify {expected_id}")
    if any(line.startswith(">") for line in lines[1:]):
        raise ValueError("multiple FASTA sequences are not supported")
    sequence = "".join(lines[1:]).upper()
    if not sequence or any(char.isspace() for char in sequence):
        raise ValueError("invalid FASTA sequence")
    return sequence


def reverse_complement(sequence):
    if set(sequence) - set("ACGT"):
        raise ValueError("reference must contain unambiguous DNA")
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def translate(sequence):
    """NCBI standard genetic code (table 1), including terminal stops."""
    if len(sequence) % 3 or set(sequence) - set("ACGT"):
        raise ValueError("CDS must have complete, unambiguous codons")
    alphabet = "TCAG"
    residues = (
        "FFLLSSSSYY**CC*W"
        "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG"
    )
    code = dict(zip(
        ("".join(codon) for codon in itertools.product(alphabet, repeat=3)),
        residues,
    ))
    return "".join(code[sequence[index:index + 3]]
                   for index in range(0, len(sequence), 3))


def reconcile(contents):
    annotation = json.loads(contents["annotation.json"])
    candidates = [row for row in annotation.get("ncbiRefSeq", [])
                  if row.get("name") == TRANSCRIPT and row.get("name2") == "HBB"
                  and row.get("chrom") == "chr11"]
    if len(candidates) != 1:
        raise ValueError("need exactly one primary chr11 HBB RefSeq transcript")
    row = candidates[0]
    if annotation.get("genome") != "hg38" or row["strand"] != "-":
        raise ValueError("annotation must identify GRCh38 negative-strand HBB")
    if row["cdsStartStat"] != "cmpl" or row["cdsEndStat"] != "cmpl":
        raise ValueError("CDS must have complete start and stop annotations")
    starts = [int(value) for value in row["exonStarts"].rstrip(",").split(",")]
    ends = [int(value) for value in row["exonEnds"].rstrip(",").split(",")]
    if len(starts) != len(ends) or len(starts) != row["exonCount"]:
        raise ValueError("exon count mismatch")
    exons = list(zip(starts, ends))
    if any(start >= end for start, end in exons):
        raise ValueError("empty or reversed exon")
    if any(left[1] > right[0] for left, right in zip(exons, exons[1:])):
        raise ValueError("exons must be ascending and non-overlapping")
    if not (SEQUENCE_START <= row["txStart"] <= row["cdsStart"]
            < row["cdsEnd"] <= row["txEnd"] <= SEQUENCE_END):
        raise ValueError("annotation lies outside reference interval")
    if starts[0] != row["txStart"] or ends[-1] != row["txEnd"]:
        raise ValueError("transcript bounds do not match terminal exons")
    genomic = fasta(contents["genomic.fasta"],
                    f"{ACCESSION}:{SEQUENCE_START + 1}-{SEQUENCE_END}")
    transcript = fasta(contents["transcript.fasta"], TRANSCRIPT + " ")
    protein = fasta(contents["protein.fasta"], PROTEIN + " ")
    if len(genomic) != SEQUENCE_END - SEQUENCE_START:
        raise ValueError("genomic reference interval length mismatch")

    def subsequence(start, end):
        if not SEQUENCE_START <= start <= end <= SEQUENCE_END:
            raise ValueError("sequence extraction outside reference")
        return genomic[start - SEQUENCE_START:end - SEQUENCE_START]

    spliced = reverse_complement("".join(subsequence(start, end)
                                         for start, end in exons))
    if spliced != transcript:
        raise ValueError("spliced genomic exons disagree with versioned RefSeq mRNA")
    cds_segments = [(max(start, row["cdsStart"]), min(end, row["cdsEnd"]))
                    for start, end in exons]
    cds = reverse_complement("".join(subsequence(start, end)
                                    for start, end in cds_segments if start < end))
    translated = translate(cds)
    if not cds.startswith("ATG") or translated != protein + "*":
        raise ValueError("spliced CDS disagrees with versioned protein or terminal stop")
    genbank = contents["transcript.gb"].decode("ascii")
    cds_match = re.search(r"^     CDS\s+(\d+)\.\.(\d+)\s*$",
                          genbank, re.MULTILINE)
    if not cds_match or f'/protein_id="{PROTEIN}"' not in genbank:
        raise ValueError("RefSeq GenBank CDS/protein metadata unavailable")
    cds_tx_start, cds_tx_end = (int(cds_match[1]) - 1, int(cds_match[2]))
    if transcript[cds_tx_start:cds_tx_end] != cds:
        raise ValueError("GenBank CDS interval disagrees with genomic model")
    transcript_model = dict(
        accession=TRANSCRIPT, protein_accession=PROTEIN, strand=row["strand"],
        tx_start=row["txStart"], tx_end=row["txEnd"],
        cds_start=row["cdsStart"], cds_end=row["cdsEnd"],
        exons=[list(exon) for exon in exons],
    )
    validation = dict(
        genomic_bases=len(genomic), transcript_bases=len(transcript),
        cds_bases=len(cds), protein_residues=len(protein),
        transcript_matches_spliced_genome=True,
        protein_matches_spliced_cds=True,
        transcript_cds_start=cds_tx_start, transcript_cds_end=cds_tx_end,
        terminal_stop_included_in_cds=True,
    )
    return transcript_model, validation


def build(output, fetch=download):
    output = pathlib.Path(output)
    if output.exists():
        raise ValueError("output already exists; use a new reference bundle directory")
    contents = {name: fetch(url) for name, url in SOURCES.items()}
    transcript_model, validation = reconcile(contents)
    retrieved_at = datetime.datetime.now(datetime.timezone.utc).isoformat()
    manifest = dict(
        schema_version=1, assembly="GRCh38", gene="HBB", chromosome="11",
        sequence_start=SEQUENCE_START, sequence_end=SEQUENCE_END,
        transcript=transcript_model,
        sources=[dict(file=name, url=url, retrieved_at=retrieved_at,
                      sha256=sha256(contents[name]))
                 for name, url in SOURCES.items()],
        files={name: sha256(content) for name, content in sorted(contents.items())},
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".hbb-reference-", dir=output.parent) as temporary:
        temporary = pathlib.Path(temporary)
        for name, content in contents.items():
            (temporary / name).write_bytes(content)
        (temporary / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n")
        if output.exists():
            raise ValueError("output appeared during download; refusing to overwrite")
        temporary.rename(output)
    print(json.dumps(validation, sort_keys=True))
    return manifest, validation


def verify(directory):
    """Offline verification of a frozen bundle, including sequence reconciliation."""
    directory = pathlib.Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text())
    if set(manifest["files"]) != set(SOURCES):
        raise ValueError("unexpected bundle files")
    contents = {}
    root = directory.resolve()
    for name, expected in manifest["files"].items():
        path = (directory / name).resolve()
        if path.parent != root:
            raise ValueError("bundle path escapes directory")
        contents[name] = path.read_bytes()
        if sha256(contents[name]) != expected:
            raise ValueError(f"checksum mismatch: {name}")
    transcript_model, validation = reconcile(contents)
    if (manifest["schema_version"] != 1 or manifest["assembly"] != "GRCh38"
            or manifest["gene"] != "HBB" or manifest["chromosome"] != "11"
            or manifest["sequence_start"] != SEQUENCE_START
            or manifest["sequence_end"] != SEQUENCE_END
            or manifest["transcript"] != transcript_model):
        raise ValueError("manifest disagrees with reconciled reference")
    for source in manifest["sources"]:
        if (source["url"] != SOURCES[source["file"]]
                or source["sha256"] != manifest["files"][source["file"]]):
            raise ValueError("source provenance mismatch")
    if {source["file"] for source in manifest["sources"]} != set(SOURCES):
        raise ValueError("source inventory mismatch")
    return validation


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", help="New bundle directory, or existing one with --verify")
    parser.add_argument("--verify", action="store_true", help="Verify a frozen bundle offline")
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(verify(args.output), sort_keys=True))
    else:
        build(args.output)
