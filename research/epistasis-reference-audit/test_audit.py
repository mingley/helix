"""Synthetic cases for reference-audit decisions; no reporter outcomes."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('reference_audit', Path(__file__).with_name('audit.py'))
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
helpers = audit.mapper()


class FakeReference:
    def __init__(self, sequence):
        self.sequence = sequence
        self.calls = []

    def fetch(self, chrom, start, end, *, strand):
        self.calls.append((chrom, start, end, strand))
        return self.sequence


def example():
    rr = 'A' * 200
    ar = rr[:1] + 'C' + rr[2:]
    ra = rr[:4] + 'G' + rr[5:]
    aa = ar[:4] + 'G' + ar[5:]
    entry = {'variant1': 'chr1:102:A:C', 'variant2': 'chr1:105:A:G',
             'center': 'var1', 'window': 'left', 'distance': 3,
             'position1_zero_based': 1, 'position2_zero_based': 4,
             'rows': {'RR': 3, 'AR': 4, 'RA': 5, 'AA': 6}}
    return entry, {'RR': rr, 'AR': ar, 'RA': ra, 'AA': aa}


def blocks(*rows):
    return {'chr1': helpers.Blocks(rows)}


class ReferenceAuditControls(unittest.TestCase):
    def setUp(self):
        self.entry, self.seqs = example()
        self.ref = FakeReference(self.seqs['RR'])
        self.plus = (0, 400, 'chr2', 10000, '+', 1000, 'a')

    def run_one(self, index):
        return audit.audit_one(self.entry, self.seqs, index, self.ref, helpers)

    def test_valid_forward_and_reverse_strand(self):
        result = self.run_one(blocks(self.plus))
        self.assertEqual(result['status'], 'valid')
        self.assertEqual(result['mismatches'], [])
        self.assertEqual(self.ref.calls[-1], ('chr2', 1100, 1300, '+'))
        reverse = (0, 400, 'chr2', 10000, '-', 1000, 'b')
        result = self.run_one(blocks(reverse))
        self.assertEqual(result['status'], 'valid')
        self.assertEqual(self.ref.calls[-1], ('chr2', 8700, 8900, '-'))

    def test_signed_coordinate_disagreement_rejected(self):
        # The absolute distances still agree, as in the previous audit.
        self.entry['variant2'] = 'chr1:99:A:G'
        result = self.run_one(blocks(self.plus))
        self.assertEqual(result['status'], 'inconsistent_start')
        self.assertIsNone(result['hg19']['end'])
        self.assertIsNone(result['mismatches'])
        self.assertEqual(self.ref.calls, [])

    def test_negative_common_start_rejected(self):
        self.entry['variant1'] = 'chr1:1:A:C'
        self.entry['variant2'] = 'chr1:4:A:G'
        result = self.run_one(blocks(self.plus))
        self.assertEqual(result['status'], 'inconsistent_start')
        self.assertIsNone(result['hg19']['end'])
        self.assertEqual(result['mapping_candidates'], [])

    def test_partial_block_is_not_a_whole_mapping(self):
        partial = (0, 299, 'chr2', 10000, '+', 1000, 'a')
        self.assertEqual(self.run_one(blocks(partial))['status'], 'no_map')
        self.assertEqual(self.run_one({})['status'], 'no_map')
        self.assertEqual(self.ref.calls, [])

    def test_duplicate_mapping_is_deduplicated_but_distinct_mapping_is_ambiguous(self):
        duplicate = (0, 400, 'chr2', 10000, '+', 1000, 'another-chain')
        self.assertEqual(self.run_one(blocks(self.plus, duplicate))['status'], 'valid')
        other = (0, 400, 'chr3', 10000, '+', 2000, 'c')
        result = self.run_one(blocks(self.plus, other))
        self.assertEqual(result['status'], 'ambiguous_map')
        self.assertEqual(len(result['mapping_candidates']), 2)
        self.assertIsNone(result['mismatches'])

    def test_nonprimary_and_full_insert_mismatch_are_preserved(self):
        nonprimary = (0, 400, 'chr2_alt', 10000, '+', 1000, 'a')
        self.assertEqual(self.run_one(blocks(nonprimary))['status'], 'nonprimary')
        self.ref.sequence = self.seqs['RR'][:10] + 'T' + self.seqs['RR'][11:]
        result = self.run_one(blocks(self.plus))
        self.assertEqual(result['status'], 'full_reference_mismatch')
        self.assertEqual(result['mismatches'], [{'position_zero_based': 10, 'source': 'A', 'hg38': 'T'}])

    def test_malformed_square_is_a_binding_error(self):
        self.seqs['AA'] = self.seqs['RR']
        with self.assertRaisesRegex(ValueError, 'four-construct square'):
            self.run_one(blocks(self.plus))


if __name__ == '__main__':
    unittest.main()
