#!/usr/bin/env python3
"""Compare audited S18 constructs with released OL27 design references."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import urllib.request

import inventory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
S18 = ROOT / "research/epistasis-source-audit/siraj2026-S18.xlsx"
S18_SHA256 = "92929b4cef82797237c4c2be770feb744daa97255ac6dbfe3d53d5a34c996fbb"
S18_READER = ROOT / "research/epistasis-feature-collisions-independent/verify.py"
S18_READER_SHA256 = "f07ecefaf4231d1d03196c899ad27ddc42deddc294ac1ebfeb39a88738646233"
LOOKUP_READER = ROOT / "research/epistasis-replicate-independent/verify.py"
LOOKUP_READER_SHA256 = "c37987951b395c61cd64980d538fd2327929ff41820444eac0d645da41adb768"
ENCODE = "ENCFF715XTT.fasta.gz"
IGVF = "IGVFFI9985FTZH.tsv.gz"
LOOKUP = "ENCFF823YWM.tsv.gz"
DESIGN_FIELDS = (
    "name", "sequence", "category", "class", "source", "ref", "chr", "start",
    "end", "strand", "variant_class", "variant_pos", "SPDI", "allele", "info",
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path, algorithm):
    result = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            result.update(chunk)
    return result.hexdigest()


def load_reader(name, path, expected_sha):
    require(inventory.sha256(path) == expected_sha, f"reader hash differs: {path.name}")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_specs():
    evidence = json.loads((HERE / "evidence.json").read_text(encoding="utf-8"))
    metadata = json.loads((HERE / "ol27-source-metadata.json").read_text(encoding="utf-8"))
    encode = metadata["encode_design"]["metadata"]
    igvf = metadata["igvf_reference"]["metadata"]
    require((encode["accession"], encode["status"], encode["file_format"], encode["output_type"])
            == ("ENCFF715XTT", "released", "fasta", "elements reference"),
            "unexpected ENCODE source metadata")
    require((igvf["accession"], igvf["status"], igvf["file_format"], igvf["content_type"])
            == ("IGVFFI9985FTZH", "released", "tsv", "MPRA sequence designs"),
            "unexpected IGVF source metadata")
    require(igvf["controlled_access"] is False, "IGVF design is not public")
    specs = {
        ENCODE: ("https://www.encodeproject.org" + encode["href"], encode["file_size"]),
        IGVF: ("https://api.data.igvf.org" + igvf["href"], igvf["file_size"]),
        LOOKUP: (evidence["lookup_table"]["url"], evidence["lookup_table"]["gzip_bytes"]),
    }
    return specs, {ENCODE: encode, IGVF: igvf}


def download(directory, specs):
    for filename, (url, size) in specs.items():
        request = urllib.request.Request(url, headers={"User-Agent": "Helix-source-audit/1.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            with (directory / filename).open("xb") as out:
                total = 0
                for chunk in iter(lambda: response.read(1 << 20), b""):
                    total += len(chunk)
                    require(total <= size, f"download exceeds publisher size: {filename}")
                    out.write(chunk)
                require(total == size, f"download size differs: {filename}")


def verify_reference(path, metadata):
    require(path.stat().st_size == metadata["file_size"], f"size differs: {path.name}")
    require(digest(path, "md5") == metadata["md5sum"], f"download MD5 differs: {path.name}")
    content = hashlib.md5()
    with gzip.open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            content.update(chunk)
    require(content.hexdigest() == metadata["content_md5sum"],
            f"uncompressed content MD5 differs: {path.name}")
    return inventory.sha256(path)


def fasta_records(stream):
    name = None
    sequence = []
    for line in stream:
        line = line.rstrip("\r\n")
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                require(bool(sequence), "empty FASTA record")
                yield name, "".join(sequence)
            header = line[1:].split()
            require(bool(header), "empty FASTA header")
            name = header[0]
            sequence = []
        else:
            require(name is not None, "FASTA sequence occurs before a header")
            require(not any(ch.isspace() for ch in line), "unexpected space in FASTA sequence")
            sequence.append(line.upper())
    if name is not None:
        require(bool(sequence), "empty final FASTA record")
        yield name, "".join(sequence)


def fasta_aliases(name):
    yield name, "literal_fasta_header"
    if name.startswith("(") and name.endswith(")"):
        members = name[1:-1].split(";")
        require(len(members) >= 2 and len(set(members)) == len(members) and all(members),
                "malformed compound FASTA header")
        for alias in members:
            yield alias, "compound_fasta_member"


def tsv_records(stream):
    rows = csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE, strict=True)
    require(tuple(next(rows, ())) == DESIGN_FIELDS, "IGVF design TSV header differs")
    for row in rows:
        require(len(row) == len(DESIGN_FIELDS) and bool(row[0]) and bool(row[1]),
                "invalid IGVF design row width or identity")
        yield row[0], row[1].upper()


def select_aliases(records, expected, aliases):
    selected = defaultdict(lambda: defaultdict(Counter))
    for name, sequence in records:
        for alias, kind in aliases(name):
            if alias in expected:
                selected[alias][sequence][kind] += 1
    return selected


def evaluate_alias(sequences, expected):
    if not sequences:
        status = "missing"
    elif len(sequences) > 1:
        status = "conflicting_sequences"
    elif expected not in sequences:
        status = "sequence_mismatch"
    else:
        status = "exact_sequence"
    return {
        "status": status,
        "sequence_hashes": [
            {"sha256": hashlib.sha256(seq.encode("ascii")).hexdigest(),
             "length": len(seq), "records_by_alias_type": dict(kinds)}
            for seq, kinds in sorted(sequences.items())
        ],
    }


def validated_candidates(lookup):
    require(inventory.sha256(S18) == S18_SHA256, "S18 hash differs")
    require(inventory.sha256(inventory.DESIGN) == inventory.DESIGN_SHA, "design audit hash differs")
    require(inventory.sha256(inventory.FEATURES) == inventory.FEATURES_SHA, "feature audit hash differs")
    require(inventory.sha256(lookup) == inventory.LOOKUP_SHA, "lookup gzip hash differs")
    source_reader = load_reader("s18_sequence_reader", S18_READER, S18_READER_SHA256)
    lookup_reader = load_reader("ol27_lookup_reader", LOOKUP_READER, LOOKUP_READER_SHA256)
    prior = json.loads(inventory.DESIGN.read_text(encoding="utf-8"))
    groups = prior["eligible_groups"]
    classes = inventory.feature_classes(inventory.FEATURES, groups)
    names = {inventory.key_for(group): inventory.candidate_names(group) for group in groups}
    all_names = {name for entry in names.values() for name in entry.values()}
    seen, strict, totals, lookup_md5, lookup_bytes = lookup_reader.scan_lookup(lookup, all_names)
    require((lookup_md5, lookup_bytes) == (inventory.LOOKUP_MD5, inventory.LOOKUP_UNCOMPRESSED_BYTES),
            "uncompressed lookup bytes differ")
    original = json.loads((HERE / "ol27-name-overlap.json").read_text(encoding="utf-8"))
    require(totals == original["lookup_row_counts"], "lookup row totals differ")
    require(inventory.overlap(groups, seen, classes) == original["seen_under_any_flag"] and
            inventory.overlap(groups, strict, classes) == original["seen_on_a_flag_0_single_name_row"],
            "lookup cross-tabs differ")

    row_count, accepted, rejected = source_reader.collect_design(S18)
    require(row_count == prior["data_rows"] and
            {key: item[0] for key, item in accepted.items()} ==
            {inventory.key_for(group): group for group in groups} and
            rejected == {tuple(group["key"]): sorted(group["rows"])
                         for group in prior["ineligible_groups"]}, "S18 group replay differs")
    candidates = {}
    for key, labels in names.items():
        if classes[key] != "no_collision" and set(labels.values()) <= strict:
            candidates[key] = {
                label: (alias, accepted[key][1][label]) for label, alias in labels.items()
            }
    require(len(candidates) == 12, "unexpected number of candidate collision windows")
    return candidates, classes


def audit(directory, metadata):
    compressed = {name: verify_reference(directory / name, info)
                  for name, info in metadata.items()}
    windows, classes = validated_candidates(directory / LOOKUP)
    aliases = {alias for entries in windows.values() for alias, _ in entries.values()}
    with gzip.open(directory / ENCODE, "rt", encoding="utf-8") as stream:
        fasta = select_aliases(fasta_records(stream), aliases, fasta_aliases)

    try:
        with gzip.open(directory / IGVF, "rt", encoding="utf-8", newline="") as stream:
            tsv = select_aliases(tsv_records(stream), aliases,
                                 lambda name: [(name, "literal_igvf_name")])
        tsv_failure = None
    except (csv.Error, ValueError) as error:
        tsv = {}
        tsv_failure = str(error)

    details = []
    for key, candidates in sorted(windows.items()):
        constructs = {}
        for label, (alias, sequence) in candidates.items():
            constructs[label] = {
                "alias": alias,
                "s18_sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "encode": evaluate_alias(fasta.get(alias, {}), sequence),
                "igvf": (None if tsv_failure else evaluate_alias(tsv.get(alias, {}), sequence)),
            }
        details.append({"window": key, "feature_class": classes[key], "constructs": constructs})

    def all_exact(window, reference):
        return all(entry[reference] and entry[reference]["status"] == "exact_sequence"
                   for entry in window["constructs"].values())

    return {
        "s18_sha256": S18_SHA256,
        "lookup_gzip_sha256": inventory.LOOKUP_SHA,
        "reference_gzip_sha256": compressed,
        "candidate_collision_windows": len(windows),
        "candidate_construct_aliases": len(aliases),
        "encode_four_exact_sequences": sum(all_exact(w, "encode") for w in details),
        "igvf_four_exact_sequences": None if tsv_failure else sum(all_exact(w, "igvf") for w in details),
        "igvf_schema_failure": tsv_failure,
        "windows": details,
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    inputs = cli.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--inputs", type=Path, help="directory with the three bound gzip inputs")
    inputs.add_argument("--download", action="store_true", help="replay using a temporary directory")
    args = cli.parse_args()
    specs, metadata = source_specs()
    if args.download:
        with tempfile.TemporaryDirectory(prefix="helix-ol27-") as temp:
            directory = Path(temp)
            download(directory, specs)
            result = audit(directory, metadata)
    else:
        result = audit(args.inputs, metadata)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
