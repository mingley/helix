"""Synthetic exhaustive comparisons for DNA-only components and exposure guards."""
import copy
import hashlib
import importlib.util
from pathlib import Path
import random
import unittest

SPEC = importlib.util.spec_from_file_location("variant_group", Path(__file__).resolve().parents[2] / "scripts" / "group_variant_transfer.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def random_dna(seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(200))


def substitute(sequence, positions):
    sequence = list(sequence)
    for position in positions:
        sequence[position] = next(base for base in "ACGT" if base != sequence[position])
    return "".join(sequence)


def pair(identity, sequence, start, chromosome="chr1", alternate=None, parent=None):
    row = {"id": identity, "sequence_ref": sequence,
           "sequence_alt": substitute(sequence, [99]) if alternate is None else alternate,
           "hg38": {"chromosome": chromosome, "start": start, "end": start + 200, "strand": "+"}}
    row["gc_preserving"] = (row["sequence_ref"][99] in "GC") == (row["sequence_alt"][99] in "GC")
    if parent is not None:
        row["parent_id"] = parent
    return row


def old(identity, sequence, chromosome="chr2", start=100, alternatives=None):
    row = {"id": identity, "sequence": sequence,
           "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
           "canonical_rc_sha256": hashlib.sha256(min(sequence, M.rc(sequence)).encode()).hexdigest(),
           "assembly": "GRCh38", "coordinate_convention": "zero_based_half_open",
           "chromosome": chromosome, "start": start, "end": start + 200}
    if alternatives is not None:
        row["observed_snv_alts"] = alternatives
    return row


def brute_components(rows):
    uf = M.UnionFind(len(rows))
    for i, a in enumerate(rows):
        for j, b in enumerate(rows[:i]):
            ai, bi = a["hg38"], b["hg38"]
            related = ai["chromosome"] == bi["chromosome"] and (
                ai["start"] // M.BIN_WIDTH == bi["start"] // M.BIN_WIDTH
                or max(ai["start"], bi["start"]) < min(ai["end"], bi["end"]))
            related |= a.get("parent_id", a["id"]) == b.get("parent_id", b["id"])
            distance = min(sum(x != y for x, y in zip(s, t))
                           for s in [a["sequence_ref"], a["sequence_alt"]]
                           for t in [b["sequence_ref"], b["sequence_alt"], M.rc(b["sequence_ref"]), M.rc(b["sequence_alt"])])
            if related or distance <= 20:
                uf.union(i, j)
    return {frozenset(row["id"] for j, row in enumerate(rows) if uf.find(j) == uf.find(i))
            for i in range(len(rows))}


class VariantGroupingTests(unittest.TestCase):
    def test_packed_hamming_matches_literal_bases_and_reverse_strand(self):
        for distance in [0, 1, 20, 21, 199, 200]:
            a = random_dna(4)
            b = substitute(a, range(distance))
            self.assertEqual(M.hamming(M.packed(a), M.packed(b)), distance)
            self.assertEqual(M.hamming(M.packed(a), M.packed(M.rc(M.rc(b)))), distance)

    def test_indexed_relations_match_exhaustive_both_allele_both_strand(self):
        base = random_dna(10)
        rows = [pair("a", base, 0), pair("b", M.rc(substitute(base, range(20))), 2000000),
                pair("c", substitute(base, range(21)), 4000000)]
        rows += [pair(f"r{j:02}", random_dna(100 + j), 8000000 + j * 2000000) for j in range(15)]
        # Transitive source-parent and actual cross-bin overlaps are also present.
        rows += [pair("parent1", random_dna(50), 60000000, parent="shared"),
                 pair("parent2", random_dna(51), 64000000, parent="shared"),
                 pair("cross", random_dna(52), 999950),
                 pair("nextbin", random_dna(53), 1000100)]
        rows = M.validate_new(rows)
        uf, _ = M.connect_new(rows, M.Budget(100000))
        actual = {frozenset(row["id"] for j, row in enumerate(rows) if uf.find(j) == uf.find(i))
                  for i in range(len(rows))}
        self.assertEqual(actual, brute_components(rows))

    def test_isolated20_vs21_mismatch_threshold_and_RC(self):
        base = random_dna(55)
        for distance in [20, 21]:
            changed = substitute(base, range(distance))
            for orientation in [changed, M.rc(changed)]:
                rows = M.validate_new([pair("a", base, 0), pair("b", orientation, 2000000)])
                uf, _ = M.connect_new(rows, M.Budget(100000))
                expected = len(brute_components(rows)) == 1
                self.assertEqual(uf.find(0) == uf.find(1), expected)

    def test_union_all_before_exposure_removes_whole_transitive_component(self):
        a, b, c = random_dna(1), random_dna(2), random_dna(3)
        rows = [pair("A", a, 999950), pair("B", b, 1000100, parent="parent"),
                pair("C", c, 6000000, parent="parent"), pair("safe", random_dna(4), 10000000)]
        guards = {"explored": M.GuardIndex([old("old", random_dna(9), "chr1", 999000)])}
        result = M.group(rows, [], [], guards, M.Budget(100000))
        self.assertEqual([row["id"] for row in result["retained_pairs"]], ["safe"])
        excluded = next(c for c in result["components"] if "A" in c["variant_ids"])
        self.assertEqual(excluded["variant_ids"], ["A", "B", "C"])
        self.assertFalse(excluded["retained_after_prior_exposure"])

    def test_inclusive10000_guard_and_halfopen_adjacency(self):
        index = M.GuardIndex([old("old", random_dna(1), "chr1", 100000)])
        for gap in [0, 9999, 10000, 10001]:
            witness = index.witness({"chromosome": "chr1", "start": 100200 + gap, "end": 100400 + gap})
            self.assertEqual(witness is not None, gap <= 10000)
            if witness:
                self.assertEqual(witness["gap_bp"], gap)

    def test_alternate_only_prior_homology_excludes_component(self):
        reference = random_dna(8)
        alternate = substitute(reference, [99])
        previous = substitute(alternate, range(20))
        self.assertEqual(sum(a != b for a, b in zip(reference, previous)), 21)
        rows = [pair("new", reference, 0, alternate=alternate)]
        result = M.group(rows, [old("old", previous)], [], {}, M.Budget(100000))
        self.assertEqual(result["retained_pairs"], [])
        witness = result["components"][0]["exposure_witnesses"][0]
        self.assertEqual(witness["new_allele"], "sequence_alt")
        self.assertEqual(witness["mismatches"], 20)

    def test22segments_find_factorized21_mismatch_repair_missed_by21(self):
        reference = random_dna(12)
        positions = [j * 200 // 21 for j in range(21)]
        query = substitute(reference, positions)
        self.assertEqual(len(set(M.seeds(M.packed(reference), M.SEGMENTS21)) &
                             set(M.seeds(M.packed(query), M.SEGMENTS21))), 0)
        previous = old("factorized", reference, alternatives=[[positions[0], query[positions[0]]]])
        result = M.group([pair("new", query, 0)], [], [previous], {}, M.Budget(100000))
        self.assertEqual(result["retained_pairs"], [])
        witness = result["components"][0]["exposure_witnesses"][0]
        self.assertEqual(witness["old_DNA_kind"], "observed_single_SNV_window")
        self.assertEqual(witness["old_window_position_0based"], positions[0])

    def test_factorized_family_matches_explicit_enumeration_both_strands(self):
        reference = random_dna(28)
        alternatives = [[0, "".join(b for b in "ACGT" if b != reference[0])],
                        [199, "".join(b for b in "ACGT" if b != reference[199])]]
        previous = old("factorized", reference, alternatives=alternatives)
        family = [reference] + [reference[:p] + base + reference[p + 1:] for p, alts in alternatives for base in alts]
        for changed_positions in [range(20), range(21), range(22), range(1, 22)]:
            query = substitute(reference, changed_positions)
            expected = min(sum(a != b for a, b in zip(query, old)) for old in family) <= 20
            for oriented in [query, M.rc(query)]:
                result = M.group([pair("new", oriented, 0)], [], [previous], {}, M.Budget(100000))
                # The new alternate also must be checked by the exhaustive oracle.
                p = pair("new", oriented, 0)
                expected_pair = min(sum(a != b for a, b in zip(newseq, oldseq))
                    for newseq in [p["sequence_ref"], p["sequence_alt"], M.rc(p["sequence_ref"]), M.rc(p["sequence_alt"])]
                    for oldseq in family) <= 20
                self.assertEqual(not result["retained_pairs"], expected_pair)

    def test_old_coordinates_do_not_join_unrelated_new1Mb_components(self):
        previous = old("old", random_dna(5), "chr1", 0)
        rows = [pair("new", random_dna(6), 500000)]
        result = M.group(rows, [previous], [], {"explored": M.GuardIndex([previous])}, M.Budget(100000))
        self.assertEqual(len(result["retained_pairs"]), 1)

    def test_start_bin_crossing_does_not_join_nonoverlapping_next_bin(self):
        crossing = pair("crossing", random_dna(70), 999950)
        for start, expected in [(1000100, True), (1000150, False), (1000200, False)]:
            rows = M.validate_new([crossing, pair("neighbor", random_dna(71), start)])
            uf, stats = M.connect_new(rows, M.Budget(100000))
            self.assertEqual(uf.find(0) == uf.find(1), expected)
            self.assertEqual(stats["positive_interval_overlap_union_edges"], int(expected))

    def test_candidate_budget_aborts_and_never_silently_drops(self):
        base = random_dna(45)
        with self.assertRaisesRegex(ValueError, "budget exceeded"):
            M.group([pair("a", base, 0), pair("b", substitute(base, range(20)), 2000000)],
                    [], [], {}, M.Budget(0))

    def test_invalid_DNA_duplicate_ids_unbound_build_and_outcomes_rejected(self):
        good = pair("a", random_dna(5), 0)
        for mutation in ["dna", "position", "outcome", "duplicate"]:
            rows = [copy.deepcopy(good)]
            if mutation == "dna": rows[0]["sequence_alt"] = "N" * 200
            elif mutation == "position": rows[0]["sequence_alt"] = substitute(rows[0]["sequence_ref"], [100])
            elif mutation == "outcome": rows[0]["observed"] = 0.1
            else: rows.append(copy.deepcopy(good))
            with self.subTest(mutation=mutation), self.assertRaises(ValueError): M.validate_new(rows)
        bad = old("bad", random_dna(1))
        bad["assembly"] = "GRCh37"
        with self.assertRaises(ValueError): M.GuardIndex([bad])


if __name__ == "__main__":
    unittest.main()
