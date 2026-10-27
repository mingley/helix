"""Small, independently packed twoBit format and interval golden fixtures."""
import importlib.util
import pathlib
import struct
import sys
import tempfile
import unittest


PATH = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "twobit_reference.py"
SPEC = importlib.util.spec_from_file_location("independent_twobit_reference", PATH)
REFERENCE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = REFERENCE
SPEC.loader.exec_module(REFERENCE)


def fixture(records, endian=">", index_order=None):
    """Encode explicit fixture fields separately from the production decoder."""
    # Each item is name, packed ACGT sequence, explicit N intervals, mask intervals.
    def number(value):
        return struct.pack(endian + "I", value)

    chunks = []
    offsets = {}
    position = 16 + sum(1 + len(name) + 4 for name, *_ in records)
    for name, sequence, ns, masks in records:
        values = [len(sequence), len(ns)] + [start for start, _ in ns] + [size for _, size in ns]
        values += [len(masks)] + [start for start, _ in masks] + [size for _, size in masks] + [0]
        # Independent bit-string construction checks both byte packing and bit order.
        pairs = {"T": "00", "C": "01", "A": "10", "G": "11"}
        bits = "".join(pairs[base] for base in sequence)
        bits += "0" * ((-len(bits)) % 8)
        packed = bytes(int(bits[i:i + 8], 2) for i in range(0, len(bits), 8))
        chunk = b"".join(number(value) for value in values) + packed
        offsets[name] = position
        position += len(chunk)
        chunks.append(chunk)
    names = index_order or [record[0] for record in records]
    index = b"".join(bytes([len(name)]) + name.encode("ascii") + number(offsets[name]) for name in names)
    header = b"".join(number(value) for value in [0x1A412743, 0, len(records), 0])
    return header + index + b"".join(chunks), offsets


class TwoBitGoldenTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.path = pathlib.Path(self.temporary.name) / "reference.2bit"

    def open(self, raw):
        self.path.write_bytes(raw)
        reader = REFERENCE.TwoBitReader(self.path)
        self.addCleanup(reader.close)
        return reader

    def rejected(self, raw):
        self.path.write_bytes(raw)
        with self.assertRaises(REFERENCE.TwoBitError):
            REFERENCE.TwoBitReader(self.path)

    def test_all_byte_alignments_both_byte_orders_and_partial_last_byte(self):
        sequence = "TCAGACGTG"
        for endian, label in [(">", "big"), ("<", "little")]:
            raw, offsets = fixture([("chrTest", sequence, [], [])], endian)
            self.assertEqual(raw[offsets["chrTest"] + 16], 0x1B)  # TCAG = 00 01 10 11.
            reference = self.open(raw)
            self.assertEqual(reference.byte_order, label)
            for start in range(len(sequence) + 1):
                for end in range(start, len(sequence) + 1):
                    self.assertEqual(reference.fetch("chrTest", start, end), sequence[start:end])

    def test_n_and_softmask_blocks_and_reverse_strand(self):
        raw, _ = fixture([("chrTest", "TCAGACGTG", [(2, 3), (8, 1)], [(1, 4), (7, 2)])])
        reference = self.open(raw)
        self.assertEqual(reference.fetch("chrTest"), "TCNNNCGTN")
        self.assertEqual(reference.fetch("chrTest", masked=True), "TcnnnCGtn")
        self.assertEqual(reference.fetch("chrTest", strand="-", masked=True), "naCGnnngA")
        self.assertEqual(reference.fetch("chrTest", 3, 8), "NNCGT")
        self.assertEqual(reference.fetch("chrTest", 3, 8, masked=True), "nnCGt")

    def test_sequence_index_order_can_differ_from_physical_record_order(self):
        raw, _ = fixture([("a", "AAA", [], []), ("b", "CCCCG", [], [])], index_order=["b", "a"])
        reference = self.open(raw)
        self.assertEqual(list(reference.sequence_sizes), ["b", "a"])
        self.assertEqual(reference.fetch("a"), "AAA")
        self.assertEqual(reference.fetch("b"), "CCCCG")

    def test_zero_length_sequence_and_zero_length_intervals(self):
        raw, _ = fixture([("empty", "", [], []), ("one", "A", [], [])])
        reference = self.open(raw)
        self.assertEqual(reference.fetch("empty"), "")
        self.assertEqual(reference.fetch("one", 1, 1), "")
        self.assertEqual(reference.fetch("one", 0, 1, strand="-"), "T")

    def test_coordinates_fail_without_clamping_or_chromosome_aliases(self):
        raw, _ = fixture([("chr1", "ACGT", [], [])])
        reference = self.open(raw)
        for start, end in [(-1, 2), (0, 5), (3, 2), (True, 3), (0, 1.0)]:
            with self.subTest(start=start, end=end), self.assertRaises(ValueError):
                reference.fetch("chr1", start, end)
        with self.assertRaises(KeyError):
            reference.fetch("1", 0, 1)
        with self.assertRaises(ValueError):
            reference.fetch("chr1", 0, 1, strand="unknown")
        with self.assertRaises(ValueError):
            reference.fetch("chr1", 0, 1, masked="yes")

    def test_closed_reader_rejects_fetch_and_close_is_idempotent(self):
        raw, _ = fixture([("chr1", "ACGT", [], [])])
        reference = self.open(raw)
        reference.close()
        reference.close()
        with self.assertRaises(ValueError):
            reference.fetch("chr1")

    def test_bad_signature_version_reserved_and_short_headers_fail(self):
        raw, _ = fixture([("chr1", "ACGT", [], [])])
        self.rejected(b"not two bit")
        self.rejected(b"XXXX" + raw[4:])
        self.rejected(raw[:4] + struct.pack(">I", 1) + raw[8:])
        self.rejected(raw[:12] + struct.pack(">I", 1) + raw[16:])

    def test_duplicate_names_and_record_offsets_fail(self):
        raw, offsets = fixture([("aa", "ACGT", [], []), ("bb", "TCAG", [], [])])
        second_name = 16 + 1 + 2 + 4 + 1
        self.rejected(raw[:second_name] + b"aa" + raw[second_name + 2:])
        second_offset = second_name + 2
        self.rejected(raw[:second_offset] + struct.pack(">I", offsets["aa"]) + raw[second_offset + 4:])

    def test_offsets_inside_index_and_beyond_file_fail(self):
        raw, _ = fixture([("chr1", "ACGT", [], [])])
        offset_position = 16 + 1 + 4
        for bad in [0, 16, len(raw), 0xFFFFFFFF]:
            self.rejected(raw[:offset_position] + struct.pack(">I", bad) + raw[offset_position + 4:])

    def test_truncated_packed_data_and_overlapping_record_offsets_fail(self):
        raw, _ = fixture([("a", "ACGTACGT", [], [])])
        self.rejected(raw[:-1])
        raw, offsets = fixture([("aa", "ACGTACGT", [], []), ("bb", "ACGT", [], [])])
        second_offset = 16 + 1 + 2 + 4 + 1 + 2
        self.rejected(raw[:second_offset] + struct.pack(">I", offsets["aa"] + 4) + raw[second_offset + 4:])

    def test_n_blocks_must_be_sorted_disjoint_nonempty_and_within_sequence(self):
        for blocks in [[(3, 2)], [(0, 0)], [(2, 1), (1, 1)], [(0, 2), (1, 1)]]:
            with self.subTest(blocks=blocks):
                raw, _ = fixture([("a", "ACGT", blocks, [])])
                self.rejected(raw)

    def test_masks_have_same_bounds_and_order_requirements(self):
        for masks in [[(4, 1)], [(0, 0)], [(2, 1), (1, 1)], [(0, 3), (2, 1)]]:
            with self.subTest(masks=masks):
                raw, _ = fixture([("a", "ACGT", [], masks)])
                self.rejected(raw)

    def test_truncated_block_arrays_and_nonzero_record_reserved_fail(self):
        raw, offsets = fixture([("a", "ACGT", [(1, 1)], [])])
        record = offsets["a"]
        self.rejected(raw[:record + 12])
        # Record has length/count/starts/sizes/maskcount/reserved then packed DNA.
        reserved_offset = record + 20
        self.rejected(raw[:reserved_offset] + struct.pack(">I", 1) + raw[reserved_offset + 4:])

    def test_malicious_counts_fail_before_large_allocation(self):
        raw, offsets = fixture([("a", "ACGT", [], [])])
        self.rejected(raw[:8] + struct.pack(">I", 0xFFFFFFFF) + raw[12:])
        ncount_offset = offsets["a"] + 4
        self.rejected(raw[:ncount_offset] + struct.pack(">I", 0xFFFFFFFF) + raw[ncount_offset + 4:])

    def test_empty_whitespace_and_non_ascii_sequence_names_fail(self):
        raw, _ = fixture([("a", "ACGT", [], [])])
        self.rejected(raw[:16] + b"\0" + raw[17:])
        for bad in [b" ", b"\n", b"\xff"]:
            self.rejected(raw[:17] + bad + raw[18:])


if __name__ == "__main__":
    unittest.main()
