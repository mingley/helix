import struct
import tempfile
import unittest
from pathlib import Path

from verify import TwoBit, expected_record, parse_chain, square


def example_square():
    key = ("chr1:100:A:C", "chr1:103:A:G", "var1", "middle")
    reference = "A" * 200
    first = reference[:20] + "C" + reference[21:]
    second = reference[:23] + "G" + reference[24:]
    double = first[:23] + "G" + first[24:]
    return key, [(3, 3, "RR", reference), (4, 3, "AR", first),
                 (5, 3, "RA", second), (6, 3, "AA", double)]


def make_twobit(path, endian):
    # Unaligned extraction, N-block restoration and soft masking are tested.
    sequence = "TCAGACGTTAGA"
    values = {b: i for i, b in enumerate("TCAG")}
    packed = bytes(sum(values[sequence[i + j]] << (6 - 2 * j) for j in range(4))
                   for i in range(0, len(sequence), 4))
    def ints(*values):
        return struct.pack(endian + str(len(values)) + "I", *values)
    name = b"chr1"
    header = ints(0x1A412743, 0, 1, 0)
    index = bytes([len(name)]) + name + ints(len(header) + 1 + len(name) + 4)
    # One N block at index 9; one soft-mask block (ignored by this audit).
    record = ints(len(sequence), 1, 9, 1, 1, 4, 2, 0) + packed
    path.write_bytes(header + index + record)


class DesignControls(unittest.TestCase):
    def test_square_and_wrong_diplotype(self):
        key, rows = example_square()
        record = square(key, rows)
        self.assertEqual((record["position1_zero_based"], record["position2_zero_based"]), (20, 23))
        swapped = [(n, d, {"AR": "RA", "RA": "AR"}.get(label, label), sequence)
                   for n, d, label, sequence in rows]
        self.assertIsNone(square(key, swapped))
        self.assertIsNone(square(key, rows[:-1]))
        self.assertIsNone(square((key[0], "chr2:103:A:G", *key[2:]), rows))

    def test_double_must_be_exact(self):
        key, rows = example_square()
        n, d, label, seq = rows[-1]
        rows[-1] = (n, d, label, seq[:150] + "T" + seq[151:])
        self.assertIsNone(square(key, rows))


class ChainControls(unittest.TestCase):
    def test_plus_minus_and_boundaries(self):
        index = parse_chain("""chain 0 chr1 1000 + 100 120 chr2 2000 + 300 320 1
20
chain 0 chr3 1000 + 100 120 chr4 2000 - 300 320 2
20
""".splitlines())
        self.assertEqual(index.map("chr1", 102, 107),
                         [{"chromosome": "chr2", "start": 302, "end": 307, "strand": "+"}])
        self.assertEqual(index.map("chr3", 102, 107),
                         [{"chromosome": "chr4", "start": 1693, "end": 1698, "strand": "-"}])
        self.assertEqual(index.map("chr1", 100, 120)[0]["start"], 300)
        self.assertEqual(index.map("chr1", 99, 120), [])
        self.assertEqual(index.map("chr1", 100, 121), [])

    def test_gap_overlap_and_duplicate_mappings(self):
        index = parse_chain("""chain 0 chr1 1000 + 100 125 chr2 2000 + 300 326 1
10 5 6
10
chain 0 chr1 1000 + 100 110 chr2 2000 + 300 310 2
10
chain 0 chr1 1000 + 100 110 chr3 2000 + 400 410 3
10
""".splitlines())
        mappings = index.map("chr1", 100, 110)
        self.assertEqual(len(mappings), 2)  # exact duplicate removed; distinct map retained
        self.assertEqual([x["chromosome"] for x in mappings], ["chr2", "chr3"])
        self.assertEqual(index.map("chr1", 109, 116), [])  # must not bridge a gap
        self.assertEqual(index.map("chr1", 115, 125)[0]["start"], 316)

    def test_reject_bad_or_unterminated_chain(self):
        with self.assertRaises(ValueError):
            parse_chain(["chain 0 chr1 1000 + 100 120 chr2 2000 + 300 320 1", "19"])
        with self.assertRaises(ValueError):
            parse_chain(["chain 0 chr1 1000 + 100 120 chr2 2000 + 300 320 1", "10 1 1"])


class ReferenceControls(unittest.TestCase):
    def test_both_byte_orders_n_blocks_and_reverse_strand(self):
        with tempfile.TemporaryDirectory() as directory:
            for endian in (">", "<"):
                path = Path(directory) / ("ref-" + ("big" if endian == ">" else "little"))
                make_twobit(path, endian)
                reader = TwoBit(path)
                try:
                    self.assertEqual(reader.fetch("chr1", 0, 12, "+"), "TCAGACGTTNGA")
                    self.assertEqual(reader.fetch("chr1", 1, 10, "+"), "CAGACGTTN")
                    self.assertEqual(reader.fetch("chr1", 1, 10, "-"), "NAACGTCTG")
                    with self.assertRaises(ValueError):
                        reader.fetch("chr1", 0, 13, "+")
                    with self.assertRaises(ValueError):
                        reader.fetch("chr2", 0, 2, "+")
                finally:
                    reader.close()

    def test_outcome_statuses_and_full_mismatch_positions(self):
        key, rows = example_square()
        item = square(key, rows)
        class Chains:
            def __init__(self, candidates):
                self.candidates = candidates
            def map(self, _chrom, _start, _end):
                return self.candidates
        class Reference:
            def __init__(self, value):
                self.value = value
            def fetch(self, *_args):
                return self.value
        primary = {"chromosome": "chr1", "start": 10, "end": 210, "strand": "+"}
        alt = {"chromosome": "chr1_alt", "start": 10, "end": 210, "strand": "+"}
        ref = Reference("A" * 200)
        self.assertEqual(expected_record(item, Chains([]), ref)["status"], "no_map")
        self.assertEqual(expected_record(item, Chains([primary, alt]), ref)["status"], "ambiguous_map")
        self.assertEqual(expected_record(item, Chains([alt]), ref)["status"], "nonprimary")
        self.assertEqual(expected_record(item, Chains([primary]), ref)["status"], "valid")
        mismatch = expected_record(item, Chains([primary]), Reference("T" + "A" * 199))
        self.assertEqual(mismatch["status"], "full_reference_mismatch")
        self.assertEqual(mismatch["mismatches"], [{"position_zero_based": 0, "source": "A", "hg38": "T"}])
        item["position2_zero_based"] += 1
        self.assertEqual(expected_record(item, Chains([primary]), ref)["status"], "inconsistent_start")


if __name__ == "__main__":
    unittest.main()
