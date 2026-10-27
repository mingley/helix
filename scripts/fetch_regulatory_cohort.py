#!/usr/bin/env python3
"""Metadata-first acquisition of a fixed, licensed Kircher2019 promoter cohort."""
import argparse
import collections
import concurrent.futures
import csv
import datetime
import hashlib
import io
import json
import math
import pathlib
import re
import shutil
import urllib.request

API = "https://api.mavedb.org/api/v1/score-sets/"
DOI = "10.1038/s41467-019-11526-w"
PMID = "31395865"
HBB = "urn:mavedb:00000018-a-1"
HEADERS = ["accession", "hgvs_nt", "hgvs_splice", "hgvs_pro", "score", "p-value", "unique_tags", "dna_sequences", "rna_sequences"]
PRIMARY = ["00000015-a-1", "00000016-a-1", "00000017-a-1", "00000019-a-1", "00000020-a-1", "00000023-a-1", "00000024-a-1", "00000027-a-1", "00000031-a-1"]
VALIDATION = ["00000023-a-2", "00000027-b-1", "00000031-b-1", "00000031-c-1", "00000031-d-1"]
COHORT = sorted(PRIMARY + VALIDATION + ["00000018-a-1"])


def digest(data):
    return hashlib.sha256(data).hexdigest()


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def fetch(url, body=None):
    req = urllib.request.Request(url, data=None if body is None else json.dumps(body).encode(),
                                 headers={} if body is None else {"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as response:
        data = response.read(4_000_001)
        if len(data) > 4_000_000:
            raise ValueError("source response exceeds bounded acquisition size")
        return data, dict(url=url, method="GET" if body is None else "POST", retrieved_at=now(),
                          content_type=response.headers.get("Content-Type", ""), bytes=len(data), sha256=digest(data))


def reverse(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def numeric(raw):
    value = raw.strip()
    if not value or value.upper() in ("NA", "N/A", "NULL", "."):
        return dict(raw=raw, value=None, status="missing")
    try:
        number = float(value)
    except ValueError:
        return dict(raw=raw, value=None, status="non_numeric")
    if not math.isfinite(number):
        return dict(raw=raw, value=None, status="non_finite")
    return dict(raw=raw, value=number, status="finite")


def classify(hgvs, target):
    match = re.fullmatch(r"n\.([1-9][0-9]*)([ACGT])>([ACGT])", hgvs)
    if match:
        position, ref, alt = match.groups()
        position = int(position)
        if not 1 <= position <= len(target) or target[position-1] != ref or ref == alt:
            raise ValueError("source SNV coordinate/reference mismatch")
        return "single_nucleotide_substitution", (position, ref, alt)
    equal = re.fullmatch(r"n\.([1-9][0-9]*)=", hgvs)
    if equal:
        if not 1 <= int(equal.group(1)) <= len(target):
            raise ValueError("source equality outside target")
        return "source_equality", None
    if not hgvs.startswith("n."):
        return "unsupported_coordinate_system", None
    if "delins" in hgvs or ";" in hgvs or "[" in hgvs or "_" in hgvs:
        return "unsupported_complex_allele", None
    if "del" in hgvs:
        return "unsupported_deletion", None
    if "ins" in hgvs or "dup" in hgvs:
        return "unsupported_insertion_or_duplication", None
    return "unsupported_or_malformed_description", None


def metadata_checks(metadata):
    urn = metadata["urn"]
    if urn.split(":")[-1] not in COHORT or "promoter" not in metadata["title"].lower():
        raise ValueError("score set outside fixed promoter accession inventory")
    license_info = metadata["license"]
    if (license_info["shortName"], license_info["version"], license_info["link"]) != (
            "CC0", "1.0", "https://creativecommons.org/publicdomain/zero/1.0/"):
        raise ValueError("source score-set data permission is not CC0 1.0")
    targets = metadata["targetGenes"]
    if len(targets) != 1 or targets[0]["category"] != "regulatory":
        raise ValueError("target must be one regulatory fragment")
    target = targets[0]["targetSequence"]
    sequence = target["sequence"]
    if target["taxonomy"]["code"] != 9606 or target["sequenceType"] != "dna" or set(sequence) - set("ACGT") or not 150 <= len(sequence) <= 1000:
        raise ValueError("unsupported target organism/sequence/length")
    method = metadata["methodText"]
    if "Log2 variant expression effect" not in method or "linear model" not in method:
        raise ValueError("source endpoint is not the declared log2 regression coefficient")
    if not any(p.get("identifier") == PMID for p in metadata["experiment"]["primaryPublicationIdentifiers"]):
        raise ValueError("metadata does not identify the fixed primary study")
    coords = metadata["extraMetadata"]
    if coords["reference"] != "hg38" or coords["end"] - coords["start"] + 1 != len(sequence):
        raise ValueError("unsupported genomic coordinate declaration")
    return sequence


def cell_context(metadata):
    text = metadata["shortDescription"]
    names = [name for name in ("HepG2", "HeLa", "HEL 92.1.7", "HEK293T", "K562", "SF7996") if name in text]
    if len(names) != 1:
        raise ValueError("cell context must identify exactly one declared human cell line")
    return names[0]


def map_target(metadata, reference):
    target = metadata_checks(metadata)
    coords = metadata["extraMetadata"]
    if reference["genome"] != "hg38" or reference["chrom"] != "chr" + coords["chr"]:
        raise ValueError("reference assembly/chromosome mismatch")
    sequence = reference["dna"].upper()
    if len(sequence) != reference["end"] - reference["start"]:
        raise ValueError("reference response length mismatch")
    matches = [(reference["start"]+i, strand) for i in range(len(sequence))
               for fragment, strand in ((target, "+"), (reverse(target), "-")) if sequence.startswith(fragment, i)]
    if len(matches) != 1:
        raise ValueError("fragment does not match uniquely within bounded reference window")
    start, strand = matches[0]
    end = start + len(target)
    if (start, end) != (coords["start"]-1, coords["end"]):
        raise ValueError("source coordinates do not reconcile as one-based inclusive")
    return dict(assembly="GRCh38", chromosome=coords["chr"], start=start, end=end,
                target_reference_alignment_strand=strand,
                source_coordinate_convention="one-based inclusive, resolved by exact full-target sequence mapping",
                genome_wide_uniqueness_established=False)


def metadata_phase(output):
    output = pathlib.Path(output)
    if output.exists():
        raise ValueError("output exists; choose a new snapshot directory")
    payloads, sources, assays = {}, [], []
    query = dict(publicationIdentifiers=[PMID], published=True, limit=100, offset=0)
    raw, info = fetch(API + "search", query)
    search = json.loads(raw)
    if search["numScoreSets"] != len(search["scoreSets"]):
        raise ValueError("metadata discovery is paginated; expand protocol before selecting")
    payloads["discovery-search.json"] = raw
    info.update(file="discovery-search.json", request_body=query, role="metadata_only_cohort_discovery")
    sources.append(info)
    references = {}
    for short in COHORT:
        urn = "urn:mavedb:" + short
        raw, info = fetch(API + urn)
        metadata = json.loads(raw)
        target = metadata_checks(metadata)
        metadata_file = short + ".metadata.json"
        payloads[metadata_file] = raw
        info.update(file=metadata_file, role="metadata_only_accession_selection")
        sources.append(info)
        coords = metadata["extraMetadata"]
        key = (coords["chr"], coords["start"], coords["end"])
        reference_file = f"reference-{key[0]}-{key[1]}-{key[2]}.json"
        if key not in references:
            url = ("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr" + key[0]
                   + f";start={key[1]-101};end={key[2]+100}")
            reference_raw, reference_info = fetch(url)
            payloads[reference_file] = reference_raw
            reference_info.update(file=reference_file, role="bounded_reference_sequence")
            sources.append(reference_info)
            references[key] = json.loads(reference_raw)
        mapping = map_target(metadata, references[key])
        element = metadata["targetGenes"][0]["name"].split(" promoter")[0]
        role = "primary" if short in PRIMARY else "previously_exposed_reference" if urn == HBB else "validation_only"
        assays.append(dict(assay_id=urn, element_id=element, role=role, title=metadata["title"],
                           cell_context=cell_context(metadata), assay_context=metadata["shortDescription"], target_sequence=target, target_length=len(target),
                           metadata_file=metadata_file, reference_file=reference_file,
                           expected_source_rows=metadata["numVariants"], mapping=mapping))
    primary = [a for a in assays if a["role"] == "primary"]
    for i, a in enumerate(primary):
        for b in primary[i+1:]:
            ma, mb = a["mapping"], b["mapping"]
            if ma["chromosome"] == mb["chromosome"] and max(ma["start"],mb["start"]) < min(ma["end"],mb["end"]):
                raise ValueError("overlapping primary promoter loci violate inclusion rule")
            if a["target_sequence"] in (b["target_sequence"], reverse(b["target_sequence"])):
                raise ValueError("duplicate primary target sequence")
    selection = dict(schema_version=1, phase="metadata_only_no_score_endpoints_requested", study_doi=DOI,
                     selection_rule="All nine distinct eligible non-HBB human promoters from the publication; lexicographically first accession per locus is primary. Five other same-locus contexts/replicates are validation-only. Existing HBB is pre-exposed reference-only.",
                     primary_elements=9, new_score_sets=14, all_score_sets=15, assays=assays,
                     discovery_total_score_sets=search["numScoreSets"], source_roles=sources,
                     outcome_access_before_freeze="No new score CSV endpoints were requested in this acquisition phase; existing HBB outcomes were previously exposed.")
    payloads["metadata-selection.json"] = (json.dumps(selection,indent=2)+"\n").encode()
    metadata_manifest = dict(schema_version=1, phase=selection["phase"], retained_at=now(),
                             files={name:digest(data) for name,data in payloads.items()})
    output.mkdir(parents=True)
    for name,data in payloads.items():
        (output/name).write_bytes(data)
    (output/"metadata-manifest.json").write_text(json.dumps(metadata_manifest,indent=2)+"\n")
    print(json.dumps({"primary_elements":9,"new_score_sets":14,"all_score_sets":15,"phase":selection["phase"]}))


def validate_scores(raw, assay):
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8")))
    if reader.fieldnames != HEADERS:
        raise ValueError("unsupported score CSV header")
    count = 0
    ids, representations, snvs = set(), set(), set()
    categories, statuses = collections.Counter(), collections.Counter()
    for row in reader:
        count += 1
        if count > assay["expected_source_rows"]:
            raise ValueError("source row count exceeds declared metadata")
        identifier, hgvs = row["accession"], row["hgvs_nt"]
        if not identifier.startswith(assay["assay_id"]+"#") or identifier in ids or hgvs in representations:
            raise ValueError("duplicate or invalid source accession/representation")
        ids.add(identifier);representations.add(hgvs)
        category, snv = classify(hgvs,assay["target_sequence"])
        if snv:
            if snv in snvs:
                raise ValueError("duplicate normalized target SNV")
            snvs.add(snv)
        categories[category] += 1
        statuses[numeric(row["score"])["status"]] += 1
    if count != assay["expected_source_rows"]:
        raise ValueError("source rows differ from metadata")
    return dict(all_source_rows=count, representation_counts=dict(sorted(categories.items())),
                score_status_counts=dict(sorted(statuses.items())), reference_validated_snv_rows=len(snvs),
                excluded_representation_rows=count-len(snvs),
                unique_target_positions_with_snvs=len({v[0] for v in snvs}))


def acquire(output, protocol, expected_protocol_sha256):
    output = pathlib.Path(output)
    if (output/"manifest.json").exists():
        raise ValueError("cohort already acquired; preserve it and choose a new metadata snapshot")
    phase = json.loads((output/"metadata-manifest.json").read_text())
    if phase["phase"] != "metadata_only_no_score_endpoints_requested":
        raise ValueError("missing metadata-first phase")
    for name, expected in phase["files"].items():
        if digest((output/name).read_bytes()) != expected:
            raise ValueError("metadata source changed before acquisition")
    selection = json.loads((output/"metadata-selection.json").read_text())
    protocol_raw = pathlib.Path(protocol).read_bytes()
    if digest(protocol_raw) != expected_protocol_sha256:
        raise ValueError("protocol changed after explicit acquisition authorization")
    frozen = json.loads(protocol_raw)
    roles = {a["assay_id"]:a["role"] for a in selection["assays"]}
    declared = {item["assay_id"]:item["role"] for item in frozen["cohort"]["assay_roles"]}
    if declared != roles:
        raise ValueError("frozen protocol must declare the exact metadata-selected accession role map")
    if frozen.get("status") != "frozen" or not frozen.get("frozen_before_new_outcome_acquisition"):
        raise ValueError("protocol does not freeze roles before new outcomes")
    # Permission scope is reviewed before any score request, not inferred from
    # a platform software license. Reuse the already frozen public scope proof.
    license_root = pathlib.Path("datasets/hbb-promoter-assay")
    license_raw = (license_root / "license-scope.json").read_bytes()
    license_manifest = json.loads((license_root / "manifest.json").read_text())
    if digest(license_raw) != license_manifest["files"]["license-scope.json"]:
        raise ValueError("frozen data-permission scope proof checksum mismatch")
    scope = json.loads(license_raw)
    if scope["scope_source_doi"] != "10.1186/s13059-019-1845-6" or "downloaded file" not in scope["scope_excerpt"]:
        raise ValueError("score-set data license scope is unverified")
    # This commitment is persisted BEFORE any new score endpoint is requested.
    gate = dict(schema_version=1, protocol_sha256=digest(protocol_raw), metadata_manifest_sha256=digest((output/"metadata-manifest.json").read_bytes()),
                frozen_before_new_outcome_acquisition=True, gate_written_at=now(), assay_roles=roles)
    gate_raw=(json.dumps(gate,indent=2)+"\n").encode()
    for name,data in [("frozen-protocol.json",protocol_raw),("outcome-acquisition-gate.json",gate_raw)]:
        with (output/name).open("xb") as stream:stream.write(data)
    sources = list(selection["source_roles"])
    files = dict(phase["files"])
    files["metadata-manifest.json"] = digest((output/"metadata-manifest.json").read_bytes())
    files["frozen-protocol.json"] = digest(protocol_raw)
    files["outcome-acquisition-gate.json"] = digest(gate_raw)
    reports=[]
    assays=[]
    for original in selection["assays"]:
        assay=dict(original)
        short=assay["assay_id"].split(":")[-1]
        name=short+".scores.csv"
        if assay["assay_id"] == HBB:
            raw=pathlib.Path("datasets/hbb-promoter-assay/scores.csv").read_bytes()
            info=dict(url=API+HBB+"/scores",method="retained_existing_snapshot",retrieved_at=None,
                      content_type="text/csv",bytes=len(raw),sha256=digest(raw))
        else:
            raw,info=fetch(API+assay["assay_id"]+"/scores")
            if "text/csv" not in info["content_type"]:
                raise ValueError("score endpoint did not return CSV")
        report=validate_scores(raw,assay)
        with (output/name).open("xb") as stream:stream.write(raw)
        files[name]=digest(raw)
        info.update(file=name,role="previously_exposed_reference_scores" if assay["assay_id"]==HBB else "post_protocol_measured_outcomes")
        sources.append(info)
        assay["scores_file"]=name
        assays.append(assay)
        reports.append(dict(assay_id=assay["assay_id"],element_id=assay["element_id"],role=assay["role"],**report))
    with (output/"license-scope.json").open("xb") as stream:stream.write(license_raw)
    files["license-scope.json"]=digest(license_raw)
    manifest=dict(schema_version=1,cohort_id="kircher2019-human-promoters-v1",study_doi=DOI,
                  assay_count=len(assays),primary_element_count=9,score_unit="log2_reporter_expression_coefficient",
                  endpoint="Source-fitted additive linear-model coefficient in plasmid MPRA; not endogenous expression or clinical risk",
                  independent_study_count=1,metadata_manifest_sha256=files["metadata-manifest.json"],
                  frozen_protocol_sha256=digest(protocol_raw),outcome_acquisition_gate_sha256=digest(gate_raw),
                  acquired_at=now(),assays=assays,row_accounting=reports,sources=sources,files=files,
                  limitations=["New elements are held out within one published study, not independent-study replication", "Cell context and element identity are confounded in cross-element comparisons", "Validation-only replicate/time/cell/siRNA score sets share loci with primary elements and must not enter training as independent elements", "Deposited target alignment strand does not establish physical reporter construct orientation", "Equality/unsupported rows are preserved; no neutrality or deletion is inferred from equality notation"])
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({"acquired_assays":len(assays),"primary_elements":9,"rows":sum(x["all_source_rows"] for x in reports)}))


def verify(output):
    output=pathlib.Path(output).resolve()
    manifest=json.loads((output/"manifest.json").read_text())
    for name,expected in manifest["files"].items():
        path=(output/name).resolve()
        if path.parent!=output or digest(path.read_bytes())!=expected:
            raise ValueError("source path/checksum mismatch")
    gate=json.loads((output/"outcome-acquisition-gate.json").read_text())
    protocol_raw=(output/"frozen-protocol.json").read_bytes()
    roles={a["assay_id"]:a["role"] for a in manifest["assays"]}
    declared={item["assay_id"]:item["role"] for item in json.loads(protocol_raw)["cohort"]["assay_roles"]}
    if gate["protocol_sha256"]!=digest(protocol_raw) or declared!=roles:
        raise ValueError("cohort role/protocol commitment mismatch")
    for assay,account in zip(manifest["assays"],manifest["row_accounting"]):
        metadata=json.loads((output/assay["metadata_file"]).read_text())
        reference=json.loads((output/assay["reference_file"]).read_text())
        if map_target(metadata,reference)!=assay["mapping"]:
            raise ValueError("stored reference mapping differs from exact source match")
        report=validate_scores((output/assay["scores_file"]).read_bytes(),assay)
        if report!={k:v for k,v in account.items() if k not in ("assay_id","element_id","role")}:
            raise ValueError("row accounting differs from raw CSV")
    print(json.dumps({"verified_offline":True,"assays":len(manifest["assays"]),"primary_elements":9}))


def self_test():
    assert classify("n.2C>T","ACGT")==('single_nucleotide_substitution',(2,'C','T'))
    assert classify("n.2=","ACGT")==('source_equality',None)
    for s in ('n.02C>T','c.2C>T','n.[1A>C;2C>T]','n.2del'):
        assert classify(s,'ACGT')[1] is None
    for s in ('n.0A>C','n.5A>C','n.2A>T','n.2C>C'):
        try: classify(s,'ACGT')
        except ValueError: pass
        else:
            if s!='n.0A>C':raise AssertionError(s)
    assert numeric('NaN')['status']=='non_finite' and numeric('NA')['status']=='missing'
    print('Strict representation and missingness checks passed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    for command in ('metadata','acquire','verify'):
        p=sub.add_parser(command);p.add_argument('output')
        if command=='acquire':
            p.add_argument('--protocol',required=True)
            p.add_argument('--expected-protocol-sha256',required=True)
    sub.add_parser('self-test')
    args=parser.parse_args()
    if args.command=='metadata':metadata_phase(args.output)
    elif args.command=='acquire':acquire(args.output,args.protocol,args.expected_protocol_sha256)
    elif args.command=='verify':verify(args.output)
    else:self_test()
