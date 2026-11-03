from fractions import Fraction
from itertools import product
import unittest

import analyze


class GCTrimerControls(unittest.TestCase):
    def test_all_short_sequences_and_reverse_complements(self):
        for length in (4, 5, 6):
            for bases in product("ACGT", repeat=length):
                sequence = "".join(bases)
                counts, g, t, boundary = analyze.features(sequence)
                self.assertEqual(3 * g, t + boundary)
                reverse = sequence.translate(analyze.COMPLEMENT)[::-1]
                self.assertEqual(analyze.features(reverse), (counts, g, t, boundary))

    def test_gc_is_redundant_for_interior_edits_but_not_at_edges(self):
        baseline = "A" * 200
        original = analyze.features(baseline)
        for index in (2, 100, 197):
            changed = baseline[:index] + "C" + baseline[index + 1:]
            _, g, t, boundary = analyze.features(changed)
            self.assertEqual(boundary, original[3])
            self.assertEqual(3 * (g - original[1]), t - original[2])
        for index in (0, 1, 198, 199):
            changed = baseline[:index] + "C" + baseline[index + 1:]
            _, g, t, boundary = analyze.features(changed)
            self.assertNotEqual(boundary, original[3])
            self.assertEqual(3 * (g - original[1]) - (t - original[2]), boundary - original[3])

    def test_same_trimers_can_have_distinct_gc_if_boundaries_change(self):
        first, second = "AAC" * 66 + "AA", "ACA" * 66 + "AC"
        c1, g1, t1, b1 = analyze.features(first)
        c2, g2, t2, b2 = analyze.features(second)
        self.assertEqual(c1, c2)
        self.assertEqual(t1, t2)
        self.assertEqual((g1, g2, b1, b2), (66, 67, 0, 3))

    def test_linear_scores_reparameterize_exactly(self):
        weights = {word: Fraction(i - 10, 17) for i, word in enumerate(analyze.CLASSES)}
        gc_weight = Fraction(-9, 11)
        for sequence in ("ACGT", "ACGTTAA", "AAC" * 66 + "AA", "ACA" * 66 + "AC"):
            left, right = analyze.exact_scores(sequence, weights, gc_weight)
            self.assertEqual(left, right)

    def test_rejects_invalid_sequences(self):
        for sequence in ("", "ACG", "ACGN", "acgt"):
            with self.assertRaises(ValueError):
                analyze.features(sequence)


if __name__ == "__main__":
    unittest.main()
