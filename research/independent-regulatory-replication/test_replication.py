"""Synthetic golden tests for the independent implementation, without assay fitting."""
import importlib.util
import pathlib
import sys
import unittest
import itertools
import csv
import io

import numpy as np

PATH = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "replicate_regulatory_study.py"
SPEC = importlib.util.spec_from_file_location("independent_regulatory_replication", PATH)
REPLICATION = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = REPLICATION
SPEC.loader.exec_module(REPLICATION)


class IndependentGoldenTests(unittest.TestCase):
    def test_splitmix64_matches_published_integer_stream_and_frozen_seeds(self):
        rng = REPLICATION.SplitMix64(0)
        self.assertEqual([rng.next() for _ in range(3)], [
            0xE220A8397B1DCDAF, 0x6E789E6AA1B965F4, 0x06C45D188009454F])
        rng = REPLICATION.SplitMix64(2026102702)
        self.assertEqual([rng.index(9) for _ in range(6)], [7, 0, 6, 4, 7, 7])

    def test_sparse_positions_are_retained_and_quality_rules_are_input_only(self):
        stream = io.StringIO()
        writer = csv.writer(stream)
        writer.writerow(["accession", "hgvs_nt", "score", "unique_tags", "dna_sequences", "rna_sequences"])
        writer.writerows([
            ["a", "n.1A>C", 1.5, 10, 1, 1],
            ["b", "n.1A>G", 2.5, 10, 1, 1],
            ["c", "n.1A>T", "", 10, 1, 1],
            ["d", "n.2A>C", 3.5, 9, 1, 1],
            ["e", "n.=", 0, 10, 1, 1],
        ])
        rows, accounting = REPLICATION.parse_rows(
            stream.getvalue().encode(), "AA", "assay", expected_rows=5)
        self.assertEqual([row.accession for row in rows], ["a", "b"])
        self.assertEqual(accounting["position_allele_counts"], {1: 2})
        self.assertEqual(accounting["incomplete_positions"], [1])
        self.assertEqual(len(accounting["excluded_rows"]), 3)

    def test_grouped_metrics_weight_positions_not_rows(self):
        rows = [REPLICATION.Observation("a", str(i), p, "A", a, 0, f"n.{p}A>{a}")
                for i, (p, a) in enumerate([(1, "C"), (1, "G"), (1, "T"), (2, "C")])]
        result = REPLICATION.grouped_metrics(rows, np.zeros(4), np.array([1, 1, 1, 3]))
        self.assertEqual(result["mse"], 5.0)
        self.assertEqual(result["mae"], 2.0)

    def test_sparse_group_sizes_do_not_create_ranks_for_constant_predictions(self):
        rows = [REPLICATION.Observation("a", str(i), p, "A", a, 0, f"n.{p}A>{a}")
                for i, (p, a) in enumerate([(1, "C"), (1, "G"), (1, "T"),
                                           (2, "C"), (2, "G"), (3, "C")])]
        result = REPLICATION.grouped_metrics(
            rows, np.arange(6, dtype=float), np.full(6, 0.1))
        self.assertIsNone(result["spearman_rows"])
        self.assertIsNone(result["spearman_positions"])
        self.assertIsNone(REPLICATION.metrics(np.full(6, 0.1), np.arange(6))["r_squared"])

    def test_observed_position_ranks_preserve_exact_source_decimal_ties(self):
        rows = [
            REPLICATION.Observation("a", "1", 1, "A", "C", 0.1, "n.1A>C", score_raw="0.1"),
            REPLICATION.Observation("a", "2", 1, "A", "G", 0.2, "n.1A>G", score_raw="0.2"),
            REPLICATION.Observation("a", "3", 2, "A", "C", 0.15, "n.2A>C", score_raw="0.15"),
        ]
        result = REPLICATION.grouped_metrics(rows, np.array([0.1, 0.2, 0.15]), np.array([0, 0, 1]))
        self.assertIsNone(result["spearman_positions"])

    def test_identical_feature_inputs_have_one_exact_prediction_per_fit(self):
        train_x = np.arange(48, dtype=float).reshape(4, 12) / 13.0
        test_x = np.tile(train_x, (257, 1))
        fit = REPLICATION.CachedFit.build(
            np.vstack([train_x, test_x]), np.arange(4), np.arange(4, 4 + len(test_x)),
            list(range(12)), alpha=10)
        values = np.arange(4 + len(test_x), dtype=float)
        prediction, _, _ = fit.predict(values)
        for row in range(4):
            self.assertEqual(len(set(prediction[row::4])), 1)

    def test_rotations_are_chemistry_stratified_and_singletons_consume_no_draw(self):
        specs = [("a", 1, "A", "C"), ("a", 2, "A", "C"), ("a", 3, "A", "C"),
                 ("a", 1, "A", "G"), ("a", 2, "A", "G"), ("a", 1, "A", "T"),
                 ("b", 1, "C", "A"), ("b", 2, "C", "A"), ("b", 3, "C", "A")]
        rows = [REPLICATION.Observation(a, str(i), p, r, t, i, f"n.{p}{r}>{t}")
                for i, (a, p, r, t) in enumerate(specs)]
        protocol = {"spatial_null": {"seed": 0, "replicates": 1},
                    "features": {"mutation_identity_order": [
                        r + ">" + a for r in "ACGT" for a in "ACGT" if r != a]}}
        rotated, _ = REPLICATION.rotation_targets(rows, np.arange(9), protocol)
        np.testing.assert_array_equal(rotated[:, 0], [2, 0, 1, 4, 3, 5, 8, 6, 7])

    def test_feature_vocabulary_and_signed_three_mer_changes(self):
        protocol = {
            "features": {
                "mutation_identity_order": [r + ">" + a for r in "ACGT" for a in "ACGT" if r != a],
                "delta3mer_order": ["".join(s) for s in itertools.product("ACGT", repeat=3)]},
            "models": [{"id": "full", "feature_sets": [
                "mutation_identity", "position", "delta3mer", "motif_deltas"]}],
        }
        rows = [REPLICATION.Observation("a", "id", 2, "A", "C", 0, "n.2A>C")]
        profiles = [{"matrix_id": str(i), "pfm": {base: [1] for base in "ACGT"}}
                    for i in range(4)]
        matrix, models, names = REPLICATION.feature_matrix("AAA", rows, profiles, protocol)
        self.assertEqual(matrix.shape, (1, 82))
        self.assertEqual(len(models["full"]), 82)
        self.assertEqual(matrix[0, names.index("mutation:A>C")], 1)
        self.assertEqual(matrix[0, names.index("position:x")], 0)
        self.assertEqual(matrix[0, names.index("delta3mer:AAA")], -1)
        self.assertEqual(matrix[0, names.index("delta3mer:ACA")], 1)
        self.assertEqual(matrix[0, 14:78].sum(), 0)

    def test_cached_batch_refits_recompute_each_training_mean(self):
        x = np.array([[0.0], [1.0], [2.0], [3.0]])
        fit = REPLICATION.CachedFit.build(x, [0, 1, 2], [3], [0], alpha=3)
        y = np.column_stack([np.arange(4), np.arange(4) + 5])
        prediction, means, _ = fit.predict(y)
        np.testing.assert_allclose(prediction, [[2.0, 7.0]], atol=1e-12)
        np.testing.assert_array_equal(means, [1.0, 6.0])

    def test_numerically_constant_decimal_feature_uses_unit_scale(self):
        prediction, fit = REPLICATION.ridge_predict(
            np.full((1000, 1), 0.1), np.full(1000, 7.0), np.array([[100.0]]), 10)
        self.assertEqual(fit["training_feature_scale"], [1.0])
        np.testing.assert_array_equal(prediction, [7.0])

    def test_folds_group_substitutions_and_exclude_complete_guard_zone(self):
        positions = np.repeat(np.arange(1, 13), 3)
        folds = REPLICATION.blocked_folds(positions, [(1, 4), (5, 8), (9, 12)], guard=1)
        self.assertEqual(set(positions[folds[1]["test_indices"]]), {5, 6, 7, 8})
        self.assertEqual(set(positions[folds[1]["guard_indices"]]), {4, 9})
        self.assertEqual(set(positions[folds[1]["train_indices"]]), {1, 2, 3, 10, 11, 12})
        self.assertEqual(len(folds[1]["test_indices"]), 12)
        with self.assertRaises(ValueError):
            REPLICATION.blocked_folds(positions, [(1, 6), (6, 12)], guard=0)

    def test_motif_delta_has_hand_computed_one_column_value(self):
        counts = {"A": [3.0], "C": [0.0], "G": [0.0], "T": [0.0]}
        result = REPLICATION.motif_delta("A", 0, "C", counts)
        # Both orientations are scored. A/T max has p=.8125; C/G both p=.0625.
        self.assertAlmostEqual(result["reference_max"], np.log2(3.25))
        self.assertAlmostEqual(result["alternate_max"], -2.0)
        self.assertAlmostEqual(result["delta"], -2.0 - np.log2(3.25))
        self.assertEqual(result["coordinate_orientation_candidates"], 2)

    def test_dual_strand_motif_delta_is_reverse_complement_invariant(self):
        sequence = "AACGTTAG"
        counts = {"A": [9.0, 1.0, 0.0], "C": [0.0, 9.0, 1.0],
                  "G": [1.0, 0.0, 9.0], "T": [0.0, 0.0, 0.0]}
        forward = REPLICATION.motif_delta(sequence, 3, "A", counts)
        reverse = REPLICATION.motif_delta(
            sequence.translate(REPLICATION.COMPLEMENT)[::-1], len(sequence) - 1 - 3,
            "A".translate(REPLICATION.COMPLEMENT), counts,
        )
        self.assertAlmostEqual(forward["delta"], reverse["delta"])
        self.assertEqual(forward["coordinate_orientation_candidates"],
                         reverse["coordinate_orientation_candidates"])

    def test_ridge_penalty_and_training_only_scaling_have_closed_form(self):
        prediction, fit = REPLICATION.ridge_predict(
            np.array([[0.0], [1.0], [2.0]]),
            np.array([0.0, 1.0, 2.0]),
            np.array([[3.0], [1000.0]]),
            alpha=3.0,
        )
        np.testing.assert_allclose(prediction, [2.0, 500.5], atol=1e-12)
        np.testing.assert_allclose(fit["training_feature_mean"], [1.0], atol=1e-12)
        np.testing.assert_allclose(fit["training_feature_scale"], [np.sqrt(2.0 / 3.0)], atol=1e-12)

    def test_intercept_is_not_penalized_and_constant_features_do_not_explode(self):
        prediction, fit = REPLICATION.ridge_predict(
            np.array([[5.0], [5.0], [5.0]]), np.array([12.0, 12.0, 12.0]),
            np.array([[-999.0], [999.0]]), alpha=1e9,
        )
        np.testing.assert_array_equal(prediction, [12.0, 12.0])
        self.assertEqual(fit["training_feature_scale"], [1.0])

    def test_empty_feature_model_is_training_mean(self):
        prediction, _ = REPLICATION.ridge_predict(
            np.empty((3, 0)), np.array([2.0, 4.0, 9.0]), np.empty((2, 0)), alpha=1.0,
        )
        np.testing.assert_array_equal(prediction, [5.0, 5.0])

    def test_metrics_match_exact_errors_and_average_rank_ties(self):
        result = REPLICATION.metrics(np.array([1.0, 2.0, 2.0, 4.0]), np.array([4.0, 3.0, 2.0, 1.0]))
        self.assertAlmostEqual(result["mse"], 19.0 / 4.0)
        self.assertAlmostEqual(result["mae"], 7.0 / 4.0)
        self.assertAlmostEqual(result["spearman"], -0.9486832980505138)
        self.assertIsNone(REPLICATION.metrics(np.ones(3), np.zeros(3))["pearson"])

    def test_fit_rejects_bad_data_and_unfrozen_outcomes_before_file_reads(self):
        with self.assertRaises(ValueError):
            REPLICATION.ridge_predict(np.zeros((2, 1)), np.array([1.0, np.nan]), np.zeros((1, 1)), 1.0)
        with self.assertRaises(ValueError):
            REPLICATION.checked_assay(pathlib.Path("/not/an/input"), protocol_frozen=False)


if __name__ == "__main__":
    unittest.main()
