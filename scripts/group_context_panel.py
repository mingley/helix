#!/usr/bin/env python3
"""Outcome-independent genomic/sequence families and guarded context-study folds."""
from bisect import bisect_right
from collections import defaultdict
import hashlib
import json
from pathlib import Path


class UnionFind:
    def __init__(self, n):
        self.parents = list(range(n))
        self.sizes = [1] * n

    def find(self, i):
        while i != self.parents[i]:
            self.parents[i] = self.parents[self.parents[i]]
            i = self.parents[i]
        return i

    def union(self, i, j):
        i, j = self.find(i), self.find(j)
        if i != j:
            if self.sizes[i] < self.sizes[j]:
                i, j = j, i
            self.parents[j] = i
            self.sizes[i] += self.sizes[j]


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def packed(sequence):
    value = 0
    for base in sequence:
        value = (value << 2) | "ACGT".index(base)
    return value


def packed_hamming(left, right, mask):
    difference = left ^ right
    return ((difference | (difference >> 1)) & mask).bit_count()


def sequence_components(records, bin_width=1_000_000, candidate_budget=100_000_000):
    """Union all declared relations. No numerical outcome field is accessed.

    Equal-length <=floor(L/10)-mismatch neighbors must share one of k+1
    disjoint segments. Both strand orientations are indexed, making the
    candidate search exhaustive for the declared ungapped Hamming relation.
    """
    rows = sorted(records, key=lambda row: row["id"])
    if len({r["id"] for r in rows}) != len(rows) or not rows:
        raise ValueError("empty or duplicate sequence IDs")
    for row in rows:
        seq = row["sequence"]
        if (not isinstance(seq, str) or len(seq) < 3 or set(seq) - set("ACGT")
                or "\n" in row["id"] or row["start"] < 0
                or row["end"] <= row["start"]
                or row["end"] - row["start"] != len(seq)):
            raise ValueError("invalid exact genomic sequence contract: " + row["id"])
    uf = UnionFind(len(rows))
    first_bin, first_parent, first_sequence = {}, {}, {}
    unique = []
    for i, row in enumerate(rows):
        parent = row.get("parent_id") or row["id"]
        if parent in first_parent:
            uf.union(i, first_parent[parent])
        first_parent.setdefault(parent, i)
        for index in range(row["start"] // bin_width, (row["end"] - 1) // bin_width + 1):
            key = row["chromosome"], index
            if key in first_bin:
                uf.union(i, first_bin[key])
            first_bin.setdefault(key, i)
        canonical = min(row["sequence"], reverse_complement(row["sequence"]))
        if canonical in first_sequence:
            uf.union(i, first_sequence[canonical])
        else:
            unique.append(i)
            first_sequence[canonical] = i
    # Positive interval overlaps necessarily share at least one touched bin.
    segment_index = defaultdict(list)
    encodings, masks = {}, {}
    candidates = comparisons = joined = 0
    # Index all source IDs, so the frozen candidate-pair budget has the same
    # meaning even when several IDs have exact/RC-identical sequences.
    for position, i in enumerate(range(len(rows))):
        seq = rows[i]["sequence"]
        length = len(seq)
        tolerance = length // 10
        segments = tolerance + 1
        bounds = [(j * length // segments, (j + 1) * length // segments)
                  for j in range(segments)]
        peers = set()
        for j, (start, end) in enumerate(bounds):
            peers.update(segment_index.get((length, j, seq[start:end]), ()))
        candidates += len(peers)
        if candidates > candidate_budget:
            raise ValueError("homology candidate budget exceeded; abort without dropping candidates")
        forward, reverse = packed(seq), packed(reverse_complement(seq))
        encodings[i] = forward, reverse
        mask = masks.setdefault(length, int("01" * length, 2))
        for peer in sorted(peers):
            if uf.find(i) == uf.find(peer):
                continue
            comparisons += 1
            p, q = encodings[peer]
            if min(packed_hamming(forward, p, mask), packed_hamming(forward, q, mask)) <= tolerance:
                uf.union(i, peer)
                joined += 1
        for oriented in {seq, reverse_complement(seq)}:
            for j, (start, end) in enumerate(bounds):
                segment_index[(length, j, oriented[start:end])].append(i)
        if position and position % 10_000 == 0:
            print(f"Sequence-family indexing: {position}/{len(rows)}; {candidates} candidates", flush=True)
    members = defaultdict(list)
    for i, row in enumerate(rows):
        members[uf.find(i)].append(row["id"])
    components = []
    for ids in members.values():
        ids.sort()
        digest = hashlib.sha256("\n".join(ids).encode()).hexdigest()
        components.append({"group": digest, "ids": ids, "size": len(ids)})
    components.sort(key=lambda item: (-item["size"], item["group"]))
    return rows, components, {
        "records": len(rows), "unique_exact_or_rc_sequences": len(unique),
        "candidate_pairs": candidates, "cross_component_hamming_comparisons": comparisons,
        "homology_union_edges": joined, "components": len(components),
        "bin_width": bin_width, "similarity": "equal-length ungapped >=90% identity in either strand",
        "candidate_budget": candidate_budget,
    }


def guarded_partition(rows, components, folds=5, guard_bp=10_000,
                      min_components=100, max_component_fraction=0.4,
                      min_training_fraction=0.5):
    if (len(components) < min_components
            or components[0]["size"] / len(rows) > max_component_fraction):
        raise ValueError("component feasibility gate failed before fitting")
    assignments = [[] for _ in range(folds)]
    loads = [0] * folds
    groups = {}
    for component in components:
        selected = min(range(folds), key=lambda f: (loads[f], f))
        assignments[selected].append(component["group"])
        loads[selected] += component["size"]
        for identity in component["ids"]:
            groups[identity] = component["group"]
    result = []
    for fold, test_groups in enumerate(assignments):
        test_set = set(test_groups)
        intervals = defaultdict(list)
        for row in rows:
            if groups[row["id"]] in test_set:
                intervals[row["chromosome"]].append((row["start"], row["end"]))
        search = {}
        for chromosome, spans in intervals.items():
            spans.sort()
            starts, ends, maximum = [], [], -1
            for start, end in spans:
                starts.append(start)
                maximum = max(maximum, end)
                ends.append(maximum)
            search[chromosome] = starts, ends
        excluded_groups = set()
        for row in rows:
            group = groups[row["id"]]
            if group in test_set or row["chromosome"] not in search:
                continue
            starts, ends = search[row["chromosome"]]
            j = bisect_right(starts, row["end"] + guard_bp) - 1
            if j >= 0 and ends[j] >= row["start"] - guard_bp:
                excluded_groups.add(group)
        test_ids = [row["id"] for row in rows if groups[row["id"]] in test_set]
        guard_ids = [row["id"] for row in rows if groups[row["id"]] in excluded_groups]
        train_ids = [row["id"] for row in rows
                     if groups[row["id"]] not in test_set | excluded_groups]
        if not test_ids or len(train_ids) / len(rows) < min_training_fraction:
            raise ValueError("guarded training-support gate failed before fitting")
        minimum_distances = {}
        train_set = set(train_ids)
        for row in rows:
            chromosome = row["chromosome"]
            if row["id"] not in train_set or chromosome not in search:
                continue
            starts, ends = search[chromosome]
            j = bisect_right(starts, row["end"]) - 1
            distances = []
            if j >= 0:
                distances.append(max(0, row["start"] - ends[j]))
            if j + 1 < len(starts):
                distances.append(max(0, starts[j + 1] - row["end"]))
            distance = min(distances)
            minimum_distances[chromosome] = min(minimum_distances.get(chromosome, distance), distance)
        result.append({"fold": fold, "test_groups": sorted(test_groups),
                       "guard_groups": sorted(excluded_groups), "test_ids": test_ids,
                       "guard_ids": guard_ids, "train_ids": train_ids,
                       "minimum_train_test_genomic_gap_bp": min(minimum_distances.values(), default=None),
                       "minimum_train_test_genomic_gap_by_chromosome": minimum_distances})
    return {"schema_version": 1, "records": [{"id": r["id"], "group": groups[r["id"]]}
                                            for r in rows],
            "components": components, "folds": result,
            "guard_bp": guard_bp,
            "guard_convention": "Exclude whole training component if any half-open intervals have gap <= guard_bp; adjacent intervals have gap zero."}


def build(records, **options):
    rows, components, stats = sequence_components(records)
    result = guarded_partition(rows, components, **options)
    result["grouping_report"] = stats
    return result


def self_test():
    import random
    rng = random.Random(221)
    seq = "".join(rng.choice("ACGT") for _ in range(60))
    variants = [seq, reverse_complement(seq)]
    for distance in (6, 7):
        changed = list(seq)
        for i in range(distance):
            changed[i] = next(b for b in "ACGT" if b != changed[i])
        variants.append("".join(changed))
    variants += ["".join(rng.choice("ACGT") for _ in range(60)) for _ in range(30)]
    records = [{"id": str(i).zfill(3), "sequence": s, "chromosome": "chr1",
                "start": i * 2_000_000, "end": i * 2_000_000 + len(s)}
               for i, s in enumerate(variants)]
    rows, components, _ = sequence_components(records)
    independent = UnionFind(len(rows))
    for i, a in enumerate(rows):
        for j, b in enumerate(rows[:i]):
            distance = min(sum(x != y for x, y in zip(a["sequence"], orientation))
                           for orientation in (b["sequence"], reverse_complement(b["sequence"])))
            if distance <= len(a["sequence"]) // 10:
                independent.union(i, j)
    observed = {identity: c["group"] for c in components for identity in c["ids"]}
    for i, a in enumerate(rows):
        for j, b in enumerate(rows):
            if (observed[a["id"]] == observed[b["id"]]) != (independent.find(i) == independent.find(j)):
                raise ValueError("indexed homology differs from exhaustive independent enumeration")
    result = guarded_partition(rows, components, min_components=5, max_component_fraction=0.5,
                               min_training_fraction=0.4)
    if sorted(identity for fold in result["folds"] for identity in fold["test_ids"]) != sorted(observed):
        raise ValueError("test partition does not cover every source identity once")
    try:
        sequence_components(records, candidate_budget=0)
    except ValueError:
        pass
    else:
        raise ValueError("candidate budget did not abort")
    # A physical element crossing a bin boundary must unite both neighborhoods.
    crossing = [dict(records[i + 10], start=s, end=s + 60) for i, s in enumerate((999_950, 1_000_001, 999_800))]
    _, physical, _ = sequence_components(crossing)
    if len(physical) != 1:
        raise ValueError("spanning-bin union failed")
    for length in (3, 20, 21, 29, 30, 31, 60, 61, 200):
        original = "".join(rng.choice("ACGT") for _ in range(length))
        threshold = length // 10
        for distance in (threshold, threshold + 1):
            mutated = list(original)
            for j in range(distance):
                mutated[j] = next(base for base in "ACGT" if base != mutated[j])
            pair = [dict(records[0], id="a", sequence=original, start=0, end=length),
                    dict(records[1], id="b", sequence=reverse_complement("".join(mutated)),
                         start=2_000_000, end=2_000_000 + length)]
            _, actual, _ = sequence_components(pair)
            expected_joined = min(sum(a != b for a, b in zip(original, orientation))
                                  for orientation in (pair[1]["sequence"], reverse_complement(pair[1]["sequence"]))) <= threshold
            if (len(actual) == 1) != expected_joined:
                raise ValueError("isolated reverse-strand threshold failed")
    for gap in (9_999, 10_000, 10_001):
        guard_rows = [dict(records[10], id="test", start=999_900, end=999_960),
                      dict(records[11], id="near", start=999_960 + gap, end=1_000_020 + gap),
                      dict(records[12], id="remote", chromosome="chr2", start=8_000_000, end=8_000_060),
                      dict(records[13], id="far", start=5_000_000, end=5_000_060)]
        components = [{"group": "near-parent", "ids": ["near", "remote"], "size": 2},
                      {"group": "test", "ids": ["test"], "size": 1},
                      {"group": "far", "ids": ["far"], "size": 1}]
        split = guarded_partition(guard_rows, components, folds=3, min_components=3,
                                  max_component_fraction=0.5, min_training_fraction=0.1)
        fold = next(item for item in split["folds"] if "test" in item["test_ids"])
        if set(fold["guard_ids"]) != ({"near", "remote"} if gap <= 10_000 else set()):
            raise ValueError("inclusive whole-component guard boundary failed")
    print("Grouping matches brute-force strand-aware Hamming; budget, bin and fold controls passed.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--panel", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--protocol", type=Path, default=Path("experiments/context-specificity-v2.json"))
    args = parser.parse_args()
    if args.self_test:
        self_test()
    elif args.panel and args.out:
        if args.out.exists():
            raise ValueError("preserve existing grouping outputs")
        protocol_sha = hashlib.sha256(args.protocol.read_bytes()).hexdigest()
        if protocol_sha != "8f0576845fecb5cfcf9bfa354fab446d849244448a3489905a6ce33bdeb21cec":
            raise ValueError("frozen context protocol differs")
        try:
            document = json.loads(args.panel.read_bytes())
            result = build(document["records"])
            result["protocol_sha256"] = protocol_sha
            result["panel_sha256"] = hashlib.sha256(args.panel.read_bytes()).hexdigest()
            result["implementation_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            args.out.write_text(json.dumps(result, indent=2) + "\n")
        except Exception as error:
            receipt = args.out.with_suffix(args.out.suffix + ".aborted.json")
            with receipt.open("x") as handle:
                json.dump({"status": "aborted_before_fitting", "error_type": type(error).__name__,
                           "reason": str(error), "protocol_sha256": protocol_sha,
                           "panel_sha256": hashlib.sha256(args.panel.read_bytes()).hexdigest()}, handle, indent=2)
                handle.write("\n")
            raise
    else:
        parser.error("provide --self-test, or --panel and --out")
