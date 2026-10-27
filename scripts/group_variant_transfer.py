#!/usr/bin/env python3
"""Outcome-free full-population components, followed by prior-exposure exclusion.

Prepare and commit parameters before running. Every new reference-valid pair is
grouped before exposure exclusions. No outcome file is accepted by this tool.
"""
from __future__ import annotations

import argparse
from array import array
from bisect import bisect_right
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys

LENGTH = 200
MISMATCHES = 20
BIN_WIDTH = 1_000_000
GUARD = 10_000
CANDIDATE_BUDGET = 1_000_000_000
HAMMING_MASK = int("01" * LENGTH, 2)
PRIMARY = re.compile(r"chr(?:[1-9]|1\d|2[0-2]|X|Y)\Z")
HEX = re.compile(r"[0-9a-f]{64}\Z")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
BASE4 = str.maketrans("ACGT", "0123")
EXPLICIT_OLD_FILES = ["agarwal-cores.jsonl.gz", "kircher-deletion-windows200.jsonl.gz",
                      "hbb-auxiliary-windows200.jsonl.gz"]
FACTORIZED_OLD_FILE = "kircher-reference-windows200.jsonl.gz"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def strict_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(path):
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle, object_pairs_hook=strict_object,
                         parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def encode(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def read_jsonl(path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            yield json.loads(line, object_pairs_hook=strict_object,
                             parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def write_jsonl(path, rows):
    with Path(path).open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            for row in rows:
                handle.write((json.dumps(row, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode())


def label(value):
    return isinstance(value, str) and value and value.strip() == value and "\n" not in value


def valid_dna(sequence):
    return isinstance(sequence, str) and len(sequence) == LENGTH and not set(sequence) - set("ACGT")


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def packed(sequence):
    return int(sequence.translate(BASE4), 4)


def hamming(left, right):
    different = left ^ right
    return ((different | (different >> 1)) & HAMMING_MASK).bit_count()


def segment_parameters(count):
    return [(2 * (LENGTH - (j + 1) * LENGTH // count),
             (1 << (2 * ((j + 1) * LENGTH // count - j * LENGTH // count))) - 1)
            for j in range(count)]


SEGMENTS21 = segment_parameters(21)
SEGMENTS22 = segment_parameters(22)


def seeds(value, segments):
    # Every segment is at most10nt, so the segment number occupies bits >=20.
    return [(j << 20) | ((value >> shift) & mask) for j, (shift, mask) in enumerate(segments)]


class UnionFind:
    def __init__(self, size):
        self.parents = list(range(size))
        self.sizes = [1] * size

    def find(self, node):
        while node != self.parents[node]:
            self.parents[node] = self.parents[self.parents[node]]
            node = self.parents[node]
        return node

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left == right:
            return False
        if (self.sizes[left], -left) < (self.sizes[right], -right):
            left, right = right, left
        self.parents[right] = left
        self.sizes[left] += self.sizes[right]
        return True


class Budget:
    def __init__(self, limit):
        require(type(limit) is int and limit >= 0, "candidate budget must be a nonnegative integer")
        self.limit = limit
        self.counts = Counter()

    def add(self, stage, count=1):
        self.counts[stage] += count
        if sum(self.counts.values()) > self.limit:
            raise ValueError("homology candidate budget exceeded; no candidate pairs were dropped")


def validate_new(rows):
    rows = sorted(rows, key=lambda row: row["id"])
    require(bool(rows) and len(rows) < 2**32, "invalid new pair count")
    require(len({row["id"] for row in rows}) == len(rows), "duplicate new variant ID")
    for row in rows:
        require(set(row) <= {"id", "source_row_1based", "reasons", "source_hg19",
                "whole_insert_chain_candidates", "hg38", "sequence_ref", "sequence_alt",
                "gc_preserving", "parent_id"}, "unexpected input fields; outcomes are not accepted")
        require(label(row.get("id")), "invalid new variant ID")
        require(valid_dna(row.get("sequence_ref")) and valid_dna(row.get("sequence_alt")),
                "expected exact uppercase200nt ACGT allele cores")
        require([i for i, (a, b) in enumerate(zip(row["sequence_ref"], row["sequence_alt"])) if a != b] == [99],
                "new pair is not the declared centered single SNV")
        locus = row.get("hg38", {})
        require(PRIMARY.fullmatch(locus.get("chromosome", "")) and type(locus.get("start")) is int
                and type(locus.get("end")) is int and locus["start"] >= 0
                and locus["end"] - locus["start"] == LENGTH and locus.get("strand") in {"+", "-"},
                "invalid GRCh38 zero-based half-open200nt interval")
        require(type(row.get("gc_preserving")) is bool and row["gc_preserving"] == (
            (row["sequence_ref"][99] in "GC") == (row["sequence_alt"][99] in "GC")),
            "GC metadata does not match alleles")
        require("parent_id" not in row or label(row["parent_id"]), "invalid source parent ID")
    return rows


def connect_new(rows, budget, progress=None):
    """Exact connected components of every frozen new-to-new relation."""
    uf = UnionFind(len(rows))
    first_bin, first_parent, first_exact = {}, {}, {}
    encodings = []
    physical_unions = exact_unions = 0
    for index, row in enumerate(rows):
        parent = row.get("parent_id", row["id"])
        if parent in first_parent:
            physical_unions += uf.union(index, first_parent[parent])
        first_parent.setdefault(parent, index)
        locus = row["hg38"]
        key = locus["chromosome"], locus["start"] // BIN_WIDTH
        if key in first_bin:
            physical_unions += uf.union(index, first_bin[key])
        first_bin.setdefault(key, index)
        values = tuple(packed(s) for s in [row["sequence_ref"], row["sequence_alt"],
                                          rc(row["sequence_ref"]), rc(row["sequence_alt"])])
        encodings.append(values)
        for forward, reverse in [(values[0], values[2]), (values[1], values[3])]:
            canonical = min(forward, reverse)
            if canonical in first_exact:
                exact_unions += uf.union(index, first_exact[canonical])
            first_exact.setdefault(canonical, index)
    # A start-bin relation does not cover every interval crossing a boundary.
    # Join actual positive overlaps separately; touching endpoints do not overlap.
    interval_groups = defaultdict(list)
    for index, row in enumerate(rows):
        locus = row["hg38"]
        interval_groups[locus["chromosome"]].append((locus["start"], locus["end"], index))
    overlap_unions = 0
    for spans in interval_groups.values():
        maximum_end, representative = -1, None
        for start, end, index in sorted(spans):
            if start < maximum_end:
                overlap_unions += uf.union(index, representative)
            if end > maximum_end:
                maximum_end, representative = end, index
    del first_bin, first_parent, first_exact, interval_groups
    initial_roots = [uf.find(i) for i in range(len(rows))]
    initial_components = len(set(initial_roots))
    index = {}
    comparisons = joins = bucket_hits = 0
    witnesses = []
    for node, values in enumerate(encodings):
        peers = set()
        query_keys = set(seeds(values[0], SEGMENTS21)) | set(seeds(values[1], SEGMENTS21))
        root = uf.find(node)
        for key in query_keys:
            for initial_root, bucket in index.get(key, {}).items():
                # These complete buckets are already connected. Their edges
                # cannot alter components, so skipping them is exact.
                if uf.find(initial_root) != root:
                    bucket_hits += len(bucket)
                    peers.update(bucket)
        budget.add("new_new_distinct_variant_candidates", len(peers))
        for peer in sorted(peers):
            if uf.find(peer) == uf.find(node):
                continue
            comparisons += 1
            distance = min(hamming(query, target) for query in values[:2] for target in encodings[peer])
            if distance <= MISMATCHES:
                uf.union(node, peer)
                joins += 1
                witnesses.append({"left": rows[node]["id"], "right": rows[peer]["id"],
                                  "minimum_mismatches": distance})
        insertion_keys = set()
        for value in set(values):
            insertion_keys.update(seeds(value, SEGMENTS21))
        for key in insertion_keys:
            buckets = index.setdefault(key, {})
            if initial_roots[node] not in buckets:
                buckets[initial_roots[node]] = array("I")
            buckets[initial_roots[node]].append(node)
        if progress and node and node % 10000 == 0:
            progress({"phase": "new_new", "processed": node, "total": len(rows),
                      "candidates": dict(budget.counts), "homology_union_edges": joins})
    return uf, {"initial_physical_exact_components": initial_components,
                "physical_parent_start_bin_union_edges": physical_unions,
                "positive_interval_overlap_union_edges": overlap_unions,
                "exact_RC_union_edges": exact_unions, "cross_component_seed_bucket_entries": bucket_hits,
                "cross_component_hamming_checks": comparisons, "homology_union_edges": joins,
                "homology_union_witnesses": witnesses}


class GuardIndex:
    def __init__(self, rows):
        grouped = defaultdict(list)
        for row in rows:
            require(row.get("assembly") == "GRCh38" and row.get("coordinate_convention") == "zero_based_half_open",
                    "prior genomic guard uses unbound build or coordinate convention")
            require(PRIMARY.fullmatch(row.get("chromosome", "")) and type(row.get("start")) is int
                    and type(row.get("end")) is int and 0 <= row["start"] < row["end"],
                    "invalid prior interval")
            grouped[row["chromosome"]].append((row["start"], row["end"], row["id"]))
        self.data = {}
        for chromosome, spans in grouped.items():
            spans.sort()
            maximum, maxima = -1, []
            for _, end, _ in spans:
                maximum = max(maximum, end)
                maxima.append(maximum)
            self.data[chromosome] = spans, [s[0] for s in spans], maxima

    def witness(self, interval):
        if interval["chromosome"] not in self.data:
            return None
        spans, starts, maxima = self.data[interval["chromosome"]]
        position = bisect_right(starts, interval["end"] + GUARD) - 1
        while position >= 0 and maxima[position] >= interval["start"] - GUARD:
            start, end, identity = spans[position]
            if end >= interval["start"] - GUARD:
                gap = max(0, start - interval["end"], interval["start"] - end)
                return {"old_id": identity, "chromosome": interval["chromosome"],
                        "old_start": start, "old_end": end, "gap_bp": gap}
            position -= 1
        return None


def validate_old(row):
    require(label(row.get("id")) and valid_dna(row.get("sequence")), "invalid prior200nt DNA record")
    sequence = row["sequence"]
    require(row.get("sequence_sha256") == hashlib.sha256(sequence.encode()).hexdigest()
            and row.get("canonical_rc_sha256") == hashlib.sha256(min(sequence, rc(sequence)).encode()).hexdigest(),
            "prior sequence checksum differs")
    return row


def old_index(records, segments):
    index = defaultdict(lambda: array("I"))
    values = []
    for node, record in enumerate(records):
        value = packed(record["sequence"])
        values.append(value)
        for key in seeds(value, segments):
            index[key].append(node)
    return index, values


def old_match(query_sequence, record, reference_distance):
    """Resolve a factorized observed single-SNV family without false negatives."""
    if reference_distance <= MISMATCHES:
        return {"mismatches": reference_distance, "old_DNA_kind": "reference_or_explicit_edited_window"}
    if reference_distance != MISMATCHES + 1 or "observed_snv_alts" not in record:
        return None
    for position, alternatives in record["observed_snv_alts"]:
        require(type(position) is int and 0 <= position < LENGTH and isinstance(alternatives, str)
                and alternatives and set(alternatives) <= set("ACGT")
                and record["sequence"][position] not in alternatives,
                "invalid factorized observed-SNV inventory")
        if query_sequence[position] != record["sequence"][position] and query_sequence[position] in alternatives:
            return {"mismatches": MISMATCHES, "old_DNA_kind": "observed_single_SNV_window",
                    "old_window_position_0based": position, "observed_old_alternate": query_sequence[position]}
    return None


def screen_old(rows, uf, explicit_old, factorized_old, guard_sets, budget, progress=None):
    """Every retained component is exhaustively negative; excluded ones need one witness."""
    reasons, witnesses, row_coordinate_hits = defaultdict(set), defaultdict(list), {}
    for node, row in enumerate(rows):
        for collection, guard_index in guard_sets.items():
            witness = guard_index.witness(row["hg38"])
            if witness is not None:
                reason = "within_inclusive_10000bp_of_" + collection
                root = uf.find(node)
                reasons[root].add(reason)
                row_coordinate_hits.setdefault(row["id"], []).append(dict(witness, reason=reason))
                if not any(w.get("reason") == reason for w in witnesses[root]):
                    witnesses[root].append(dict(witness, new_id=row["id"], reason=reason))
    checked_rows = Counter()
    for stage, records, segments in [("old_explicit", explicit_old, SEGMENTS21),
                                     ("old_factorized_SNV", factorized_old, SEGMENTS22)]:
        index, values = old_index(records, segments)
        for node, row in enumerate(rows):
            root = uf.find(node)
            if reasons[root]:
                continue
            checked_rows[stage] += 1
            found = False
            for allele in ["sequence_ref", "sequence_alt"]:
                for orientation, query_sequence in [("+", row[allele]), ("-", rc(row[allele]))]:
                    query = packed(query_sequence)
                    peers = set()
                    for key in seeds(query, segments):
                        peers.update(index.get(key, ()))
                    budget.add(stage + "_oriented_allele_candidates", len(peers))
                    for peer in sorted(peers):
                        distance = hamming(query, values[peer])
                        match = old_match(query_sequence, records[peer], distance)
                        if match:
                            reason = "full200nt_ungapped_ge90pct_prior_DNA_match"
                            reasons[root].add(reason)
                            witnesses[root].append(dict(match, reason=reason, new_id=row["id"],
                                new_allele=allele, query_strand=orientation,
                                old_id=records[peer]["id"], old_inventory_stage=stage))
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if progress and node and node % 10000 == 0:
                progress({"phase": stage, "processed": node, "total": len(rows),
                          "candidates": dict(budget.counts)})
        del index, values
    return reasons, witnesses, row_coordinate_hits, dict(checked_rows)


def group(rows, explicit_old, factorized_old, guard_sets, budget, progress=None):
    rows = validate_new(rows)
    uf, new_stats = connect_new(rows, budget, progress)
    reasons, witnesses, row_hits, checked_rows = screen_old(
        rows, uf, explicit_old, factorized_old, guard_sets, budget, progress)
    members = defaultdict(list)
    for node, row in enumerate(rows):
        members[uf.find(node)].append(row["id"])
    components, assignments = [], {}
    for root, identities in members.items():
        identities.sort()
        identity = hashlib.sha256("\n".join(identities).encode()).hexdigest()
        component = {"component": identity, "variant_ids": identities, "size": len(identities),
                     "retained_after_prior_exposure": not reasons[root],
                     "exclusion_reasons": sorted(reasons[root]), "exposure_witnesses": witnesses[root]}
        components.append(component)
        for member in identities:
            assignments[member] = component
    components.sort(key=lambda row: row["component"])
    ledger, retained = [], []
    for row in rows:
        component = assignments[row["id"]]
        locus = row.get("parent_id", row["id"])
        ledger.append({"id": row["id"], "component": component["component"], "locus": locus,
            "component_size": component["size"], "retained_after_prior_exposure": component["retained_after_prior_exposure"],
            "component_exclusion_reasons": component["exclusion_reasons"],
            "row_coordinate_exposure_witnesses": row_hits.get(row["id"], [])})
        if component["retained_after_prior_exposure"]:
            retained.append(dict(row, component=component["component"], locus=locus))
    return {"components": components, "ledger": ledger, "retained_pairs": retained,
            "statistics": {"new_pairs": len(rows), "components": len(components),
                "retained_components": sum(row["retained_after_prior_exposure"] for row in components),
                "retained_pairs": len(retained), "prior_screen_checked_rows_by_stage": checked_rows,
                "homology_candidate_counts": dict(budget.counts), "candidate_budget": budget.limit,
                **new_stats}}


def verify_inventory(directory, expected_manifest_sha256):
    require(HEX.fullmatch(expected_manifest_sha256), "invalid prior manifest SHA256")
    directory = Path(directory)
    require(sha256(directory / "manifest.json") == expected_manifest_sha256, "prior manifest SHA256 mismatch")
    manifest = read_json(directory / "manifest.json")
    require(manifest.get("schema_version") == 1 and isinstance(manifest.get("files"), dict), "invalid prior manifest")
    required = set(EXPLICIT_OLD_FILES) | {FACTORIZED_OLD_FILE, "kircher-loci.json", "hbb-auxiliary-loci.json"}
    require(required <= set(manifest["files"]), "prior inventory lacks required exposures")
    for name, receipt in manifest["files"].items():
        require(Path(name).name == name and not (directory / name).is_symlink(), "unsafe prior inventory path")
        path = directory / name
        require(path.stat().st_size == receipt["bytes"] and sha256(path) == receipt["sha256"],
                "prior inventory member SHA256/size mismatch: " + name)
    return manifest


def preparation(args):
    require(not args.out.exists(), "preparation output already exists")
    require(sha256(args.new_pairs) == args.expected_new_sha256, "new-pair SHA256 mismatch")
    verify_inventory(args.prior_inventory, args.expected_prior_manifest_sha256)
    parameters = {"schema_version": 1, "protocol_id": "variant-transfer-v1",
        "phase": "DNA_only_parameters_before_new_outcome_parsing", "bin_width": BIN_WIDTH,
        "guard_bp_inclusive": GUARD, "sequence_length": LENGTH, "max_hamming_mismatches": MISMATCHES,
        "new_new_segments": 21, "old_explicit_segments": 21, "old_factorized_SNV_segments": 22,
        "candidate_budget": args.candidate_budget,
        "relations": ["both alleles belong to their variant; exact source parent", "same chromosome and floor(GRCh38 start/1000000)",
                      "separate positive half-open interval overlap, including across start-bin boundaries",
                      "exact or reverse-complement DNA; full200nt ungapped >=90% any allele/either strand"],
        "algorithm": "2bit packed Hamming; k+1 disjoint segments guarantee a matching seed for <=k mismatches; all candidates verified exactly",
        "factorized_SNV_rule": "Old reference windows plus observed single-SNV alternatives: <=20 reference mismatches matches reference; at21 require one allowed observed edit repairing a mismatch.22segments guarantee completeness at21.",
        "candidate_budget_rule": "Sum of distinct previous-new-variant candidate pairs outside currently connected buckets, plus oriented-new-allele/old-window candidates across old stages; abort without dropping if exceeded",
        "safe_pruning": "Skip new-new buckets already connected by union-find; after all new relations, coordinate-exposed components or components with one old-homology witness need no more old-screen queries",
        "screening_scope": "Every retained component has every member/allele/orientation screened; witnesses for excluded components are sufficient examples, not an exhaustive hit catalog",
        "new_pair_count_expected": 332719,
        "remaining_gates": ["source technical-control classes", "source/header/endpoint binding", "numeric outcome coverage"]}
    require(type(args.candidate_budget) is int and args.candidate_budget >= 0, "invalid candidate budget")
    inputs = {"new_pairs": {"path": str(args.new_pairs.resolve()), "sha256": args.expected_new_sha256},
              "prior_inventory": {"path": str(args.prior_inventory.resolve()),
                                  "manifest_sha256": args.expected_prior_manifest_sha256}}
    args.out.mkdir(parents=True)
    (args.out / "parameters.json").write_bytes(encode(parameters))
    (args.out / "inputs.json").write_bytes(encode(inputs))
    shutil.copyfile(Path(__file__), args.out / "implementation.py")
    manifest = {"schema_version": 1, "phase": "prepared_not_executed",
                "files": {path.name: sha256(path) for path in sorted(args.out.iterdir())}}
    (args.out / "preparation.json").write_bytes(encode(manifest))
    print(json.dumps({"preparation_sha256": sha256(args.out / "preparation.json")}, sort_keys=True))


def run(args):
    require(not args.out.exists(), "grouping output already exists")
    require(sha256(args.preparation / "preparation.json") == args.expected_preparation_sha256,
            "preparation SHA256 mismatch")
    manifest = read_json(args.preparation / "preparation.json")
    require(manifest.get("phase") == "prepared_not_executed" and set(manifest.get("files", {})) == {
        "parameters.json", "inputs.json", "implementation.py"}, "unexpected preparation schema")
    for name, digest in manifest["files"].items():
        require(sha256(args.preparation / name) == digest, "preparation member SHA256 mismatch")
    require(sha256(__file__) == manifest["files"]["implementation.py"], "run implementation differs from frozen preparation")
    parameters, inputs = read_json(args.preparation / "parameters.json"), read_json(args.preparation / "inputs.json")
    new_path = Path(inputs["new_pairs"]["path"])
    require(sha256(new_path) == inputs["new_pairs"]["sha256"], "new DNA input changed after preparation")
    prior_root = Path(inputs["prior_inventory"]["path"])
    verify_inventory(prior_root, inputs["prior_inventory"]["manifest_sha256"])
    args.out.mkdir(parents=True)
    for name in ["parameters.json", "inputs.json", "implementation.py", "preparation.json"]:
        shutil.copyfile(args.preparation / name, args.out / name)
    budget = Budget(parameters["candidate_budget"])
    try:
        rows = list(read_jsonl(new_path))
        require(len(rows) == parameters["new_pair_count_expected"], "new reference-valid population count changed")
        explicit = [validate_old(row) for name in EXPLICIT_OLD_FILES for row in read_jsonl(prior_root / name)]
        factorized = [validate_old(row) for row in read_jsonl(prior_root / FACTORIZED_OLD_FILE)]
        agarwal = [row for row in explicit if row.get("collection") == "Agarwal2025"]
        guard_sets = {"explored_Agarwal_insert": GuardIndex(agarwal),
                      "explored_Kircher_locus": GuardIndex(read_json(prior_root / "kircher-loci.json")["loci"]),
                      "explored_HBB_locus": GuardIndex(read_json(prior_root / "hbb-auxiliary-loci.json")["loci"])}
        result = group(rows, explicit, factorized, guard_sets, budget,
                       lambda status: print(json.dumps(status, sort_keys=True), flush=True))
        (args.out / "components.json").write_bytes(encode({"schema_version": 1, "components": result["components"]}))
        write_jsonl(args.out / "variant-ledger.jsonl.gz", result["ledger"])
        write_jsonl(args.out / "eligible-pairs.jsonl.gz", result["retained_pairs"])
        summary = {"schema_version": 1, "phase": "DNA_components_and_prior_exposure_complete_no_outcomes",
                   "grouping_and_prior_exposure_gates_passed": True, "source_cohort_gate_passed": False,
                   "remaining_gates": parameters["remaining_gates"], "statistics": result["statistics"],
                   "eligible_ID_order_sha256": hashlib.sha256("\n".join(row["id"] for row in result["retained_pairs"]).encode()).hexdigest(),
                   "preparation_sha256": args.expected_preparation_sha256}
        (args.out / "summary.json").write_bytes(encode(summary))
    except Exception as error:
        (args.out / "failure.json").write_bytes(encode({"schema_version": 1, "error_type": type(error).__name__,
            "error": str(error), "grouping_and_prior_exposure_gates_passed": False,
            "candidate_counts": dict(budget.counts), "preparation_sha256": args.expected_preparation_sha256}))
        finish_manifest(args.out)
        raise
    finish_manifest(args.out)
    print(json.dumps(summary["statistics"], sort_keys=True))


def finish_manifest(directory):
    (directory / "manifest.json").write_bytes(encode({"schema_version": 1,
        "files": {path.name: {"sha256": sha256(path), "bytes": path.stat().st_size}
                  for path in sorted(directory.iterdir()) if path.name != "manifest.json"}}))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--new-pairs", type=Path, required=True)
    prepare.add_argument("--expected-new-sha256", required=True)
    prepare.add_argument("--prior-inventory", type=Path, required=True)
    prepare.add_argument("--expected-prior-manifest-sha256", required=True)
    prepare.add_argument("--candidate-budget", type=int, default=CANDIDATE_BUDGET)
    prepare.add_argument("--out", type=Path, required=True)
    execute = commands.add_parser("run")
    execute.add_argument("--preparation", type=Path, required=True)
    execute.add_argument("--expected-preparation-sha256", required=True)
    execute.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    preparation(args) if args.command == "prepare" else run(args)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, TypeError, OSError, json.JSONDecodeError) as error:
        print(f"variant grouping failed: {error}", file=sys.stderr)
        raise SystemExit(2)
