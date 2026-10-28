import unittest
from verify import eligible


def example():
    key = ('chr1:100:A:C', 'chr1:103:A:G', 'var1', 'middle')
    ref = 'A' * 200
    one = ref[:20] + 'C' + ref[21:]
    two = ref[:23] + 'G' + ref[24:]
    both = one[:23] + 'G' + one[24:]
    items = [(3, 3, 'RR', ref), (4, 3, 'AR', one),
             (5, 3, 'RA', two), (6, 3, 'AA', both)]
    return key, items


class SyntheticSquares(unittest.TestCase):
    def test_valid_and_label_permutation(self):
        key, items = example()
        result = eligible(key, list(reversed(items)))
        self.assertEqual((20, 23), (result['position1_zero_based'], result['position2_zero_based']))
        self.assertEqual({'RR': 3, 'AR': 4, 'RA': 5, 'AA': 6}, result['rows'])

    def test_missing_or_repeated_labels(self):
        key, items = example()
        self.assertIsNone(eligible(key, items[:3]))
        items[3] = (6, 3, 'RA', items[3][3])
        self.assertIsNone(eligible(key, items))

    def test_double_must_equal_combined_edits(self):
        key, items = example()
        number, distance, label, seq = items[3]
        items[3] = (number, distance, label, seq[:100] + 'T' + seq[101:])
        self.assertIsNone(eligible(key, items))
        items[3] = (number, distance, label, items[1][3])
        self.assertIsNone(eligible(key, items))

    def test_source_distance_and_literal_alleles_must_agree(self):
        key, items = example()
        self.assertIsNone(eligible(('chr1:100:A:G', *key[1:]), items))
        self.assertIsNone(eligible((key[0], 'chr1:102:A:G', *key[2:]), items))
        self.assertIsNone(eligible((key[0], 'chr2:103:A:G', *key[2:]), items))
        self.assertIsNone(eligible((key[0], 'chr1:103:AA:G', *key[2:]), items))


if __name__ == '__main__':
    unittest.main()
