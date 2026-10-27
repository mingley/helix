#!/usr/bin/env python3
"""Freeze a small, prior-biology JASPAR motif set without reading assay outcomes."""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import urllib.request

BASE = "https://jaspar.elixir.no"
# Accession/version selection is fixed before MPRA positional pattern inspection.
PROFILES = (
    ("MA0035.4", "GATA1", 11, True),
    ("MA0493.3", "KLF1", 8, True),
    ("MA0060.3", "NFYA", 11, True),
    ("MA0108.3", "TBP", 7, False),
)
BACKGROUND = {base: 0.25 for base in "ACGT"}
SCORING = {
    "id": "jaspar-overlap-max-delta-v1",
    "pseudocount_total_per_column": 1.0,
    "pseudocount_allocation": "background_proportions",
    "background": BACKGROUND,
    "log_base": 2,
    "probability_formula": "(count[base,column]+0.25)/(sum_base count[base,column]+1.0)",
    "window_selection": "all complete motif-width windows within the supplied target that overlap the SNV position",
    "orientations": ["forward", "reverse_complement"],
    "comparison": "reference and alternate use exactly the same window-coordinate/orientation candidate set",
    "features": ["reference_max_log2_odds", "alternate_max_log2_odds", "delta_max_log2_odds"],
    "delta_formula": "alternate_max_log2_odds-reference_max_log2_odds",
    "maxima_note": "The maximizing coordinate/orientation can differ between alleles; retain both argmax locations when reporting a candidate.",
    "threshold": None,
    "ambiguous_bases": "reject",
    "snv_policy": "one unambiguous reference-matching substitution; no indels or multiallelic haplotypes",
    "edge_policy": "use only complete windows; return an explicit missing feature if no candidate window exists",
    "units": "bits of motif log-odds similarity, not binding energy, occupancy or expression",
    "maximum_motif_width": 11,
    "maximum_context_radius": 10,
    "strict_disjoint_support_center_guard": 20,
    "guard_note": "For maximum width11, each SNV feature can use sequence from position-10 through position+10. Strictly disjoint sequence supports require training/test centers more than20bp apart.",
}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def get(url):
    with urllib.request.urlopen(url, timeout=90) as response:
        payload = response.read()
    return payload, {
        "url": url,
        "retrieved_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "sha256": sha256(payload),
    }


def validate_profile(profile, accession, name, width, require_human):
    if profile.get("matrix_id") != accession or profile.get("name") != name:
        raise ValueError(f"JASPAR identity mismatch for {accession}")
    if profile.get("collection") != "CORE" or profile.get("tax_group") != "vertebrates":
        raise ValueError(f"unexpected collection/taxonomic group for {accession}")
    if profile.get("version") != int(accession.split(".")[1]):
        raise ValueError(f"version mismatch for {accession}")
    species = profile.get("species")
    if not isinstance(species, list):
        raise ValueError(f"missing species metadata array for {accession}")
    if require_human and not any(str(s.get("tax_id")) == "9606" for s in species):
        raise ValueError(f"explicit human taxon missing for {accession}")
    pfm = profile.get("pfm")
    if not isinstance(pfm, dict) or set(pfm) != set("ACGT"):
        raise ValueError(f"invalid PFM alphabet for {accession}")
    if any(len(pfm[base]) != width for base in "ACGT"):
        raise ValueError(f"unexpected PFM width for {accession}")
    for column in range(width):
        values = [pfm[base][column] for base in "ACGT"]
        if any(not isinstance(value, (int, float)) or isinstance(value, bool)
               or not math.isfinite(value) or value < 0 for value in values):
            raise ValueError(f"invalid PFM count for {accession}")
        if sum(values) <= 0:
            raise ValueError(f"empty PFM column for {accession}")


def write_new(path, data):
    with path.open("xb") as handle:
        handle.write(data)


def pretty(data):
    return (json.dumps(data, indent=2, ensure_ascii=False) + "\n").encode()


def freeze(output):
    output = Path(output)
    if output.exists():
        raise ValueError("output already exists; preserve the frozen motif set and choose a new directory")
    payloads = {}
    sources = []
    motifs = []
    for accession, name, width, require_human in PROFILES:
        url = f"{BASE}/api/v1/matrix/{accession}/"
        raw, provenance = get(url)
        profile = json.loads(raw)
        validate_profile(profile, accession, name, width, require_human)
        filename = f"{accession}.json"
        payloads[filename] = raw
        sources.append(dict(provenance, file=filename))
        species = profile["species"]
        motifs.append({
            "matrix_id": accession,
            "name": name,
            "width": width,
            "collection": profile["collection"],
            "tax_group": profile["tax_group"],
            "species": species,
            "species_caveat": None if species else "JASPAR does not specify a species for this profile; linked human TF protein identity does not establish the matrix training species.",
            "experimental_type": profile.get("type"),
            "uniprot_ids": profile.get("uniprot_ids", []),
            "pubmed_ids": profile.get("pubmed_ids", []),
            "pfm": profile["pfm"],
            "source_file": filename,
            "source_url": url,
        })
    about, provenance = get(f"{BASE}/about/")
    about_text = about.decode("utf-8")
    if "Creative Commons Attribution 4.0 International License" not in about_text:
        raise ValueError("JASPAR official license text was not verified")
    if "creativecommons.org/licenses/by/4.0" not in about_text:
        raise ValueError("JASPAR official CC BY4.0 license link was not verified")
    payloads["jaspar-about.html"] = about
    sources.append(dict(provenance, file="jaspar-about.html"))
    payloads["motifs.json"] = pretty({"schema_version": 1, "motifs": motifs, "scoring": SCORING})
    manifest = {
        "schema_version": 1,
        "name": "Prior-biology HBB promoter motif set",
        "motifs_file": "motifs.json",
        "selection": {
            "fixed_accessions": [profile[0] for profile in PROFILES],
            "reason": "Erythroid GATA1/KLF1 and canonical promoter CCAAT/TATA factors; selected from prior biology before positional MPRA score inspection.",
            "species_policy": "GATA1/KLF1/NFYA must explicitly include Homo sapiens taxid9606; TBP missing species metadata is preserved and never silently assigned human provenance.",
            "release_policy": "Exact accessions/versions are frozen. No claim is made that these are all from one named JASPAR release or are currently the latest profiles.",
            "outcome_independence": "The importer never reads assay files, scores, positions or model outcomes.",
        },
        "terms": {
            "license": "CC BY 4.0",
            "license_url": "https://creativecommons.org/licenses/by/4.0/",
            "license_evidence_file": "jaspar-about.html",
            "license_evidence_url": f"{BASE}/about/",
            "attribution": "JASPAR collaborative database; cite each versioned profile and the database publication when reusing these matrices.",
            "scope": "JASPAR database content as explicitly licensed on the official About page. Derived motif probabilities/log-odds are transformations of these PFMs, not new experimental measurements.",
        },
        "sources": sources,
        "files": {name: sha256(data) for name, data in sorted(payloads.items())},
    }
    output.mkdir(parents=True)
    for filename, data in sorted(payloads.items()):
        write_new(output / filename, data)
    write_new(output / "manifest.json", pretty(manifest))
    validate(output)
    print(json.dumps({"profiles": [m["matrix_id"] for m in motifs], "scoring": SCORING["id"]}))


def validate(directory):
    directory = Path(directory).resolve(strict=True)
    manifest = json.loads((directory / "manifest.json").read_text())
    if manifest.get("schema_version") != 1 or manifest.get("motifs_file") != "motifs.json":
        raise ValueError("unsupported motif manifest")
    expected_files = {f"{p[0]}.json" for p in PROFILES} | {"motifs.json", "jaspar-about.html"}
    if set(manifest["files"]) != expected_files:
        raise ValueError("motif manifest file inventory mismatch")
    for filename, expected in manifest["files"].items():
        path = (directory / filename).resolve(strict=True)
        if path.parent != directory or sha256(path.read_bytes()) != expected:
            raise ValueError(f"source path/hash mismatch: {filename}")
    bundle = json.loads((directory / "motifs.json").read_text())
    if bundle.get("schema_version") != 1 or bundle.get("scoring") != SCORING:
        raise ValueError("motif scoring contract changed")
    if [m.get("matrix_id") for m in bundle["motifs"]] != [p[0] for p in PROFILES]:
        raise ValueError("motif selection/order mismatch")
    for frozen, specification in zip(bundle["motifs"], PROFILES):
        accession, name, width, require_human = specification
        raw = json.loads((directory / f"{accession}.json").read_text())
        validate_profile(raw, accession, name, width, require_human)
        if (frozen["pfm"] != raw["pfm"] or frozen["species"] != raw["species"]
                or frozen["name"] != raw["name"] or frozen["width"] != width
                or frozen["source_file"] != f"{accession}.json"):
            raise ValueError(f"normalized motif disagrees with source: {accession}")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?")
    parser.add_argument("--validate-existing", metavar="DIRECTORY")
    args = parser.parse_args()
    if args.validate_existing:
        if args.output:
            parser.error("choose an output directory or --validate-existing, not both")
        validate(args.validate_existing)
        print("Motif source hashes, identities, counts and scoring contract verified.")
    elif args.output:
        freeze(args.output)
    else:
        parser.error("an output directory or --validate-existing is required")
