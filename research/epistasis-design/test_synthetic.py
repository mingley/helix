"""Synthetic checks for the source-agnostic interaction design."""

from itertools import product
import unittest


def contrast(y):
    return y[3] - y[1] - y[2] + y[0]


def rc(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical_counts(sequence):
    counts = {}
    for i in range(len(sequence) - 2):
        kmer = sequence[i:i + 3]
        key = min(kmer, rc(kmer))
        counts[key] = counts.get(key, 0) + 1
    counts["GC"] = sum(base in "GC" for base in sequence)
    return counts


def replace(sequence, i, base):
    return sequence[:i] + base + sequence[i + 1:]


class SyntheticInteractionTests(unittest.TestCase):
    def test_common_offset_and_signed_cancellation(self):
        square = (0.2, 0.7, -0.4, 1.0)
        self.assertAlmostEqual(contrast(square), contrast(tuple(x + 8 for x in square)))
        effects = (contrast(square), -contrast(square))
        self.assertAlmostEqual(sum(effects) / len(effects), 0)
        self.assertGreater(sum(x * x for x in effects) / len(effects), 0)

    def test_independent_errors_and_shared_error_bias(self):
        # Exactly average all equiprobable ±1 error pairs; no simulation error.
        for theta in (0, 0.5, -2):
            products = [(theta + a) * (theta + b) for a, b in product((-1, 1), repeat=2)]
            self.assertEqual(sum(products) / len(products), theta * theta)
            shared = [(theta + a) ** 2 for a in (-1, 1)]
            self.assertEqual(sum(shared) / len(shared), theta * theta + 1)
        # A sample of the valid product estimator is allowed to be negative.
        self.assertLess((0 + 1) * (0 - 1), 0)

    def test_separated_changes_have_zero_mixed_feature_difference(self):
        # All length-four sequences, all alternative end substitutions (distance 3).
        for bases in product("ACGT", repeat=4):
            reference = "".join(bases)
            for a in set("ACGT") - {reference[0]}:
                for b in set("ACGT") - {reference[3]}:
                    x10 = replace(reference, 0, a)
                    x01 = replace(reference, 3, b)
                    x11 = replace(x10, 3, b)
                    features = [canonical_counts(s) for s in (reference, x10, x01, x11)]
                    for key in set().union(*features):
                        self.assertEqual(contrast([f.get(key, 0) for f in features]), 0)

    def test_distance_two_counterexample(self):
        features = [canonical_counts(s).get("CAC", 0) for s in ("AAA", "CAA", "AAC", "CAC")]
        self.assertEqual(contrast(features), 1)


if __name__ == "__main__":
    unittest.main()
