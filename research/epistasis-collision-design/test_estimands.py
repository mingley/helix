"""Exact finite controls for the collision methods note; no assay data."""
from fractions import Fraction as Q
from itertools import combinations, product
import unittest


def cross_product(differences):
    n = len(differences)
    if n < 2:
        raise ValueError("at least two matched replicates are required")
    return Q(sum(differences) ** 2 - sum(x * x for x in differences), n * (n - 1))


def class_variance(values, weights):
    if len(values) != len(weights) or sum(weights) != 1:
        raise ValueError("one normalized weight is required for each value")
    mean = sum(w * x for w, x in zip(weights, values))
    return sum(w * (x - mean) ** 2 for w, x in zip(weights, values))


class CollisionEstimandControls(unittest.TestCase):
    def test_weighted_class_identity_and_overlap(self):
        values = [Q(1), Q(-1), Q(3)]
        weights = [Q(1, 2), Q(1, 3), Q(1, 6)]
        floor = class_variance(values, weights)
        pair_formula = sum(weights[j] * weights[k] * (values[j] - values[k]) ** 2
                           for j, k in combinations(range(len(values)), 2))
        prediction = Q(5, 3)
        mean = sum(w * x for w, x in zip(weights, values))
        risk = sum(w * (prediction - x) ** 2 for w, x in zip(weights, values))
        self.assertEqual(floor, pair_formula)
        self.assertEqual(risk, (prediction - mean) ** 2 + floor)

        four = [Q(0), Q(0), Q(1), Q(1)]
        four_floor = class_variance(four, [Q(1, 4)] * 4)
        sum_of_pair_floors = sum((four[j] - four[k]) ** 2 / 4
                                 for j, k in combinations(range(4), 2))
        self.assertEqual(four_floor, Q(1, 4))
        self.assertEqual(sum_of_pair_floors, 4 * four_floor)

    def test_independent_heterogeneous_replicate_errors(self):
        delta = 2
        for sizes in ((1, 2), (1, 2, 3), (1, 2, 3, 4)):
            statistics = [cross_product([delta + magnitude * sign
                                         for magnitude, sign in zip(sizes, signs)])
                          for signs in product((-1, 1), repeat=len(sizes))]
            self.assertEqual(sum(statistics) / len(statistics), delta ** 2)
        with self.assertRaisesRegex(ValueError, "at least two"):
            cross_product([1])

    def test_shared_construct_specific_input_bias(self):
        delta = 2
        statistics = []
        for shared, e1, e2, e3 in product((-1, 1), repeat=4):
            statistics.append(cross_product([delta + shared + e1,
                                             delta + shared + e2,
                                             delta + shared + e3]))
        self.assertEqual(sum(statistics) / len(statistics), delta ** 2 + 1)
        # Conditional on a specific input, the cross-product estimates its
        # conditional squared difference, not the marginal-mean difference.
        conditional = [cross_product([delta + 1 + e1, delta + 1 + e2])
                       for e1, e2 in product((-1, 1), repeat=2)]
        self.assertEqual(sum(conditional) / len(conditional), (delta + 1) ** 2)

    def test_true_common_offsets_cancel_within_replicates(self):
        delta = 3
        statistics = []
        for e11, e21, e12, e22, common_input in product((-1, 1), repeat=5):
            y11 = delta + common_input + 7 + e11
            y21 = common_input + 7 + e21
            y12 = delta + common_input - 4 + e12
            y22 = common_input - 4 + e22
            self.assertEqual(y11 - y21, delta + e11 - e21)
            self.assertEqual(y12 - y22, delta + e12 - e22)
            statistics.append(cross_product([y11 - y21, y12 - y22]))
        self.assertEqual(sum(statistics) / len(statistics), delta ** 2)

    def test_negative_covariance_and_negative_sample_estimates(self):
        delta = 2
        statistics = [cross_product([delta + z, delta - z]) for z in (-1, 1)]
        self.assertEqual(sum(statistics) / len(statistics), delta ** 2 - 1)
        self.assertEqual(cross_product([1, -1]), -1)
        # A realized positive estimate can exceed a true zero floor.
        self.assertEqual(cross_product([1, 1]), 1)

    def test_outcome_dependent_selection_can_bias_cross_products(self):
        samples = list(product((-1, 1), repeat=2))
        unselected = [cross_product(pair) for pair in samples]
        selected = [cross_product(pair) for pair in samples if pair[0] > 0 and pair[1] > 0]
        self.assertEqual(sum(unselected) / len(unselected), 0)
        self.assertEqual(sum(selected) / len(selected), 1)


if __name__ == "__main__":
    unittest.main()
