"""Independent synthetic golden cases; no source assay data is loaded."""
import importlib.util
import pathlib
import hashlib
import gzip
import json
import random
import sys
import tempfile
import unittest
import zipfile

import numpy as np


PATH = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "replicate_context_study.py"
SPEC = importlib.util.spec_from_file_location("independent_context_replication", PATH)
REPLICATION = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = REPLICATION
SPEC.loader.exec_module(REPLICATION)


class ContextGoldenTests(unittest.TestCase):
    def test_three_mer_vocabulary_has_32_distinct_reverse_complement_pairs(self):
        words = REPLICATION.KMER_ORDER
        self.assertEqual(len(words), 32)
        self.assertEqual(tuple(sorted(words)), words)
        self.assertIn("AAA", words)
        self.assertNotIn("TTT", words)
        self.assertIn("AAC", words)
        self.assertNotIn("GTT", words)
        self.assertEqual(len(REPLICATION.FEATURE_NAMES), 33)

    def test_three_mer_and_gc_features_have_hand_counted_values(self):
        values = REPLICATION.sequence_features("AAAC")
        expected = np.zeros(33)
        expected[REPLICATION.KMER_ORDER.index("AAA")] = 0.5
        expected[REPLICATION.KMER_ORDER.index("AAC")] = 0.5
        expected[-1] = 0.25
        np.testing.assert_array_equal(values, expected)

    def test_features_are_exact_reverse_complement_invariant(self):
        for sequence in ["AAA", "ACG", "AAACGTTGGCAT", "CGCGCGCG", "A" * 200]:
            np.testing.assert_array_equal(
                REPLICATION.sequence_features(sequence),
                REPLICATION.sequence_features(REPLICATION.reverse_complement(sequence)))

    def test_ambiguous_lowercase_and_short_sequences_are_rejected(self):
        for sequence in ["", "AC", "ACN", "acg", "AC-G"]:
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                REPLICATION.sequence_features(sequence)

    def test_standardized_ridge_has_closed_form_and_unpenalized_intercept(self):
        model = REPLICATION.fit_context_model(
            [[-1], [0], [1]], [-2, 0, 2], ["a"] * 3, "shared_dna", alpha=3)
        prediction = model.predict([[2]], ["a"])
        # X variance 2/3; z'z=3; z'y=4*sqrt(3/2); beta=z'y/6.
        # Ridge shrinks the unpenalized slope 2 to slope 1, so x=2 predicts 2.
        np.testing.assert_allclose(prediction, [2], rtol=0, atol=2e-15)
        shifted = REPLICATION.fit_context_model(
            [[-1], [0], [1]], [98, 100, 102], ["a"] * 3, "shared_dna", alpha=3)
        np.testing.assert_allclose(shifted.predict([[2]], ["a"]), [102], rtol=0, atol=2e-14)

    def test_scaling_uses_training_rows_only_and_zero_variance_is_safe(self):
        model = REPLICATION.fit_context_model(
            [[0, 7], [2, 7]], [0, 4], ["a"] * 2, "shared_dna", alpha=2)
        np.testing.assert_array_equal(model.shared.feature_means, [1, 7])
        np.testing.assert_array_equal(model.shared.feature_scales, [1, 1])
        before = model.diagnostics()
        model.predict([[1e12, -1e12]], ["a"])
        self.assertEqual(before, model.diagnostics())

    def test_context_mean_uses_exact_training_context_only(self):
        model = REPLICATION.fit_context_model(
            [[0], [1], [100], [200]], [8, 12, -8, -6], ["a", "a", "b", "b"],
            "context_mean", alpha=10)
        np.testing.assert_array_equal(model.predict([[1e8], [-1e8]], ["a", "b"]), [10, -7])
        with self.assertRaises(ValueError):
            model.predict([[0]], ["unseen"])

    def test_additive_model_matches_augmented_unpenalized_dummy_normal_equations(self):
        # Deliberately unbalanced X by context: y-centering alone is incorrect.
        x = np.array([[0., 2.], [1., 1.], [4., 0.], [10., 2.], [13., 3.]])
        y = np.array([1., 4., -2., 9., 5.])
        labels = ["a", "a", "a", "b", "b"]
        model = REPLICATION.fit_context_model(x, y, labels, "additive_context_dna", alpha=7)
        mean, scale = x.mean(axis=0), x.std(axis=0)
        design = np.column_stack([
            np.array([label == "a" for label in labels], dtype=float),
            np.array([label == "b" for label in labels], dtype=float), (x - mean) / scale])
        penalty = np.diag([0, 0, 7, 7])
        augmented_coef = np.linalg.solve(design.T @ design + penalty, design.T @ y)
        expected = design @ augmented_coef
        np.testing.assert_allclose(model.predict(x, labels), expected, rtol=0, atol=2e-14)
        np.testing.assert_allclose(model.shared.standardized_coefficients, augmented_coef[2:], atol=2e-14)

    def test_additive_context_offsets_never_receive_ridge_penalty(self):
        x = np.array([[0.], [1.], [2.], [4.]])
        labels = ["a", "a", "b", "b"]
        y = np.array([1., 5., -3., 2.])
        offsets = np.array([100., 100., -70., -70.])
        base = REPLICATION.fit_context_model(x, y, labels, "additive_context_dna", alpha=100)
        shift = REPLICATION.fit_context_model(x, y + offsets, labels, "additive_context_dna", alpha=100)
        np.testing.assert_allclose(shift.predict(x, labels) - base.predict(x, labels), offsets, atol=2e-14)
        np.testing.assert_allclose(shift.shared.standardized_coefficients,
                                   base.shared.standardized_coefficients, atol=2e-14)

    def test_context_specific_fits_capture_opposite_synthetic_slopes(self):
        x = np.array([[-1.], [0.], [1.]])
        y = np.array([[8., -5.], [10., -7.], [12., -9.]])
        models = REPLICATION.fit_balanced_panel(x, y, ["a", "b"], alpha=3)
        np.testing.assert_allclose(models["independent_context_dna"].predict([[2], [2]], ["a", "b"]),
                                   [12, -9], rtol=0, atol=2e-14)
        np.testing.assert_array_equal(models["additive_context_dna"].predict([[2], [2]], ["a", "b"]),
                                      [10, -7])

    def test_balanced_shared_and_additive_reuse_exact_slopes(self):
        x = np.array([[0., 1.], [1., 0.], [3., 4.], [7., 6.]])
        y = np.array([[1., 3., 5.], [4., 3., 2.], [5., 0., 7.], [8., 2., 3.]])
        models = REPLICATION.fit_balanced_panel(x, y, ["a", "b", "c"], alpha=10)
        self.assertIs(models["shared_dna"].shared.standardized_coefficients,
                      models["additive_context_dna"].shared.standardized_coefficients)
        np.testing.assert_array_equal(models["shared_dna"].shared.feature_scales,
                                      models["independent_context_dna"].by_context["a"].feature_scales)

    def test_aligned_shared_penalty_matches_pooled_context_count_times_alpha(self):
        x = np.array([[-2.], [-1.], [1.], [2.]])
        y = np.array([[1., 4.], [3., -1.], [6., 7.], [5., 8.]])
        models = REPLICATION.fit_balanced_panel(x, y, ["a", "b"], alpha=10)
        pooled = REPLICATION.fit_context_model(np.repeat(x, 2, axis=0), y.ravel(),
                                              ["a", "b"] * 4, "shared_dna", alpha=20)
        np.testing.assert_allclose(models["shared_dna"].shared.standardized_coefficients,
                                   pooled.shared.standardized_coefficients, rtol=0, atol=2e-15)
        pooled_wrong = REPLICATION.fit_context_model(np.repeat(x, 2, axis=0), y.ravel(),
                                                    ["a", "b"] * 4, "shared_dna", alpha=10)
        self.assertGreater(abs(pooled_wrong.shared.standardized_coefficients[0]
                               - models["shared_dna"].shared.standardized_coefficients[0]), 0.1)

    def test_average_context_specific_slope_and_prediction_equal_shared(self):
        x = np.array([[0., 1.], [1., 0.], [3., 4.], [7., 6.]])
        y = np.array([[1., 3., 5.], [4., 3., 2.], [5., 0., 7.], [8., 2., 3.]])
        contexts = ["a", "b", "c"]
        models = REPLICATION.fit_balanced_panel(x, y, contexts, alpha=10)
        separate_beta = np.array([fit.standardized_coefficients
                                  for fit in models["independent_context_dna"].by_context.values()])
        np.testing.assert_allclose(separate_beta.mean(axis=0),
                                   models["shared_dna"].shared.standardized_coefficients,
                                   rtol=0, atol=2e-15)
        query = np.array([[2., 5.], [9., -1.]])
        separate_prediction = np.array([
            models["independent_context_dna"].predict(query, [context] * len(query))
            for context in contexts])
        shared_prediction = models["shared_dna"].predict(query, [contexts[0]] * len(query))
        np.testing.assert_allclose(separate_prediction.mean(axis=0), shared_prediction,
                                   rtol=0, atol=2e-14)

    def test_interaction_mse_improvement_matches_pairwise_contrast_identity(self):
        x = np.array([[-2.], [-1.], [1.], [2.]])
        y = np.array([[1., 4., 7.], [3., -1., 4.], [6., 7., 1.], [5., 8., -2.]])
        contexts = ["a", "b", "c"]
        models = REPLICATION.fit_balanced_panel(x, y, contexts, alpha=10)
        query = np.array([[-1.5], [1.5]])
        observed = np.array([[2., 3., 6.], [7., 5., -1.]])
        predictions = {}
        for model in ["additive_context_dna", "independent_context_dna"]:
            predictions[model] = np.column_stack([
                models[model].predict(query, [context] * len(query)) for context in contexts])
        loss_improvement = np.mean((observed - predictions["additive_context_dna"]) ** 2) - np.mean(
            (observed - predictions["independent_context_dna"]) ** 2)
        pair_improvements = []
        for left in range(3):
            for right in range(left + 1, 3):
                true_delta = observed[:, left] - observed[:, right]
                additive_delta = predictions["additive_context_dna"][:, left] - predictions["additive_context_dna"][:, right]
                specific_delta = predictions["independent_context_dna"][:, left] - predictions["independent_context_dna"][:, right]
                pair_improvements.append(np.mean((true_delta - additive_delta) ** 2)
                                         - np.mean((true_delta - specific_delta) ** 2))
        self.assertAlmostEqual(loss_improvement, (2 / 6) * np.mean(pair_improvements), places=13)

    def test_complete_panel_preserves_missingness_and_rejects_duplicates(self):
        result = REPLICATION.paired_panel_indices(["x", "x", "y"], ["b", "a", "a"], ["a", "b"])
        self.assertEqual(result, {"x": (1, 0)})
        with self.assertRaises(ValueError):
            REPLICATION.paired_panel_indices(["x", "x"], ["a", "a"], ["a", "b"])

    def test_identity_and_genomic_overlap_form_transitive_indivisible_families(self):
        element = REPLICATION.SequenceElement
        elements = [
            element("a", "AAAC", "hg38", "chr1", 0, 4),
            element("b", "CCCC", "hg38", "chr1", 3, 7),
            element("c", "GGGG", "hg38", "chr2", 100, 104),  # RC identity links b.
            element("d", "ACAC", "hg38", "chr1", 7, 11),  # Touching alone is not overlap.
            element("e", "ACTG", "hg19", "chr1", 0, 4),   # Assembly is distinct.
        ]
        families = REPLICATION.identity_overlap_families(elements)
        self.assertEqual(families, {"a": "a", "b": "a", "c": "a", "d": "d", "e": "e"})
        self.assertEqual(families, REPLICATION.identity_overlap_families(list(reversed(elements))))

    def test_family_holdout_keeps_all_element_context_rows_together(self):
        families = ["a", "a", "a", "b", "b", "b", "c", "c", "c"]
        folds = REPLICATION.grouped_folds(families, {"a": 0, "b": 1, "c": 2}, 3)
        self.assertEqual(folds[0]["test_indices"], [0, 1, 2])
        self.assertEqual(folds[0]["train_indices"], [3, 4, 5, 6, 7, 8])
        for fold in folds:
            self.assertFalse(set(fold["train_families"]) & set(fold["test_families"]))
        with self.assertRaises(ValueError):
            REPLICATION.grouped_folds(families, {"a": 0, "b": 1}, 3)

    def test_identical_raw_inputs_predict_exactly_identically(self):
        x = np.array([[0.13, 0.7], [0.42, 0.1], [0.72, 0.9]])
        model = REPLICATION.fit_context_model(x, [1, 3, 2], ["a"] * 3, "shared_dna")
        values = model.predict(np.tile(x, (257, 1)), ["a"] * (257 * 3))
        for row in range(3):
            self.assertEqual(len(set(values[row::3])), 1)

    def test_primary_metrics_and_invalid_numeric_inputs(self):
        self.assertEqual(REPLICATION.errors([1, 3], [0, 1]), {"count": 2, "mse": 2.5, "mae": 1.5})
        with self.assertRaises(ValueError):
            REPLICATION.errors([1, np.nan], [0, 1])
        with self.assertRaises(ValueError):
            REPLICATION.fit_context_model([[0], [1]], [0, 1], ["a", "a"], "shared_dna", alpha=0)

    def test_new_hex_seed_is_distinct_and_matches_frozen_uint64_stream(self):
        rng = REPLICATION.SplitMix64(0x2026102702)
        self.assertEqual([rng.next_u64() for _ in range(3)], [
            0x311730414543F7D5, 0xB0EE1155B9089752, 0xFEA23A01DC5344F1])
        self.assertNotEqual(0x2026102702, 2026102702)

    def test_component_bootstrap_pools_elements_instead_of_group_means(self):
        loss = np.array([[1., 0, 0, 0, 0], [3., 0, 0, 0, 0], [9., 0, 0, 0, 0]])
        result = REPLICATION.paired_component_bootstrap(
            ["a", "a", "b"], loss, loss, [("global_mean", "context_mean")], 10, 0)
        # Equal-element mean=(1+3+9)/3; equal-group mean=((1+3)/2+9)/2.
        self.assertAlmostEqual(result["equal_sequences"]["comparisons"][0]["observed_mse_improvement"], 13 / 3)
        self.assertEqual(result["equal_groups"]["comparisons"][0]["observed_mse_improvement"], 5.5)
        for draw, value in zip(result["resampling_indices"], result["equal_sequences"]["sampled_mse_improvements"]):
            numerator = sum([4, 9][index] for index in draw)
            denominator = sum([2, 1][index] for index in draw)
            self.assertAlmostEqual(value[0], numerator / denominator)

    def test_common_bootstrap_draws_reused_for_all_models_and_contrasts(self):
        loss = np.array([[2., 1., 3., 4., 2.], [4., 2., 6., 8., 4.]])
        result = REPLICATION.paired_component_bootstrap(
            ["a", "b"], loss, loss / 2, [("global_mean", "context_mean"),
                                         ("additive_context_dna", "independent_context_dna")], 7, 5)
        for mse, mae in zip(result["equal_sequences"]["sampled_mse_improvements"],
                            result["equal_sequences"]["sampled_mae_improvements"]):
            self.assertAlmostEqual(mse[1], 2 * mse[0])
            np.testing.assert_allclose(mae, np.array(mse) / 2, rtol=0, atol=0)

    def test_complete_panel_oof_evaluation_keeps_test_targets_out_of_fit(self):
        rows = [
            {"id": str(i), "group": str(i), "sequence": seq,
             "observations": {"a": float(i), "b": -float(i), "c": float(i) / 2}}
            for i, seq in enumerate(["AAAA", "AAAC", "AACC", "ACCC", "CCCC", "CCCG"])]
        first = {"train_ids": ["3", "4", "5"], "test_ids": ["0", "1", "2"], "guard_ids": []}
        second = {"train_ids": ["0", "1", "2"], "test_ids": ["3", "4", "5"], "guard_ids": []}
        baseline = REPLICATION.evaluate_complete_panel(rows, ["a", "b", "c"], [first, second])
        changed = [dict(row, observations={key: value + 1000 if row["id"] in first["test_ids"] else value
                                          for key, value in row["observations"].items()}) for row in rows]
        shifted = REPLICATION.evaluate_complete_panel(changed, ["a", "b", "c"], [first, second])
        self.assertEqual(baseline["folds"][0], shifted["folds"][0])
        self.assertEqual(len(baseline["predictions"]), 18)
        self.assertLess(baseline["balanced_design_checks"]["max_abs_per_sequence_mse_contrast_identity_residual"], 1e-12)

    def test_oof_evaluation_rejects_split_family_and_incomplete_test_coverage(self):
        rows = [{"id": str(i), "group": "shared" if i < 2 else str(i), "sequence": "AAAC",
                 "observations": {"a": i, "b": -i}} for i in range(4)]
        folds = [{"train_ids": ["0", "2"], "test_ids": ["1", "3"], "guard_ids": []},
                 {"train_ids": ["1", "3"], "test_ids": ["0", "2"], "guard_ids": []}]
        with self.assertRaisesRegex(ValueError, "family leaks"):
            REPLICATION.evaluate_complete_panel(rows, ["a", "b"], folds)
        for row in rows:
            row["group"] = row["id"]
        with self.assertRaisesRegex(ValueError, "exactly once"):
            REPLICATION.evaluate_complete_panel(rows, ["a", "b"], [folds[0], folds[0]])

    def workbook_fixture(self, extra_cell=""):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        path = pathlib.Path(temporary.name) / "source.xlsx"
        main = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
        rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("xl/workbook.xml", f'<workbook xmlns="{main}" xmlns:r="{rel}"><sheets><sheet name="source" r:id="r1"/></sheets></workbook>')
            archive.writestr("xl/_rels/workbook.xml.rels", '<Relationships><Relationship Id="r1" Target="worksheets/sheet1.xml"/></Relationships>')
            archive.writestr("xl/sharedStrings.xml", f'<sst xmlns="{main}"><si><r><t>AA</t></r><r><t>CG</t></r></si></sst>')
            archive.writestr("xl/worksheets/sheet1.xml", f'<worksheet xmlns="{main}"><sheetData><row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1"><v>0.10000000000000001</v></c><c r="C1" t="e"><v>#N/A</v></c><c r="D1" t="b"><v>1</v></c><c r="E1" t="inlineStr"><is><t> Cell </t></is></c>{extra_cell}</row></sheetData></worksheet>')
        return path

    def test_independent_xlsx_preserves_raw_numeric_strings_types_and_rich_text(self):
        with REPLICATION.IndependentWorkbook(self.workbook_fixture()) as book:
            row, cells = next(book.rows("source"))
        self.assertEqual(row, 1)
        self.assertEqual(cells[0], "AACG")
        self.assertEqual(cells[1].literal, "0.10000000000000001")
        self.assertEqual(cells[1].finite_float(), 0.1)
        self.assertEqual(cells[2], REPLICATION.ExcelError("#N/A"))
        self.assertIs(cells[3], True)
        self.assertEqual(cells[4], " Cell ")

    def test_xlsx_formula_is_not_silently_treated_as_an_observed_number(self):
        path = self.workbook_fixture('<c r="F1"><f>B1+1</f><v>1.1</v></c>')
        with REPLICATION.IndependentWorkbook(path) as book, self.assertRaisesRegex(ValueError, "formula"):
            list(book.rows("source"))

    def test_full_grouping_matches_independent_exhaustive_small_graph_oracle(self):
        rng = random.Random(731)
        sequences = ["".join(rng.choice("ACGT") for _ in range(60)) for _ in range(10)]
        near = list(sequences[0])
        for index in range(6):
            near[index] = next(base for base in "ACGT" if base != near[index])
        sequences += ["".join(near), REPLICATION.reverse_complement(sequences[1])]
        rows = [{"id": str(i), "sequence": sequence, "chromosome": "chr" + str(i),
                 "start": 0, "end": 60, "parent_id": str(i)} for i, sequence in enumerate(sequences)]
        rows[3]["parent_id"] = rows[2]["parent_id"]
        rows[4].update(chromosome="bin", start=999980, end=1000040)
        rows[5].update(chromosome="bin", start=1000200, end=1000260)
        edges = {i: set() for i in range(len(rows))}
        for i, left in enumerate(rows):
            for j, right in enumerate(rows[:i]):
                same_bin = (left["chromosome"] == right["chromosome"] and
                            bool(set(range(left["start"] // 1000000, (left["end"] - 1) // 1000000 + 1)) &
                                 set(range(right["start"] // 1000000, (right["end"] - 1) // 1000000 + 1))))
                direct = sum(a != b for a, b in zip(left["sequence"], right["sequence"]))
                reverse = sum(a != b for a, b in zip(left["sequence"], REPLICATION.reverse_complement(right["sequence"])))
                if same_bin or left["parent_id"] == right["parent_id"] or min(direct, reverse) <= 6:
                    edges[i].add(j)
                    edges[j].add(i)
        unseen, expected = set(edges), []
        while unseen:
            pending = [min(unseen)]
            component = set()
            while pending:
                index = pending.pop()
                if index not in component:
                    component.add(index)
                    pending.extend(edges[index] - component)
            unseen -= component
            members = sorted(rows[i]["id"] for i in component)
            expected.append({"ids": members, "size": len(members),
                             "group": hashlib.sha256("\n".join(members).encode()).hexdigest()})
        expected.sort(key=lambda item: (-item["size"], item["group"], item["ids"]))
        actual = REPLICATION.reconstruct_grouping(rows, minimum_components=1, largest_component_fraction=1,
                                                  minimum_training_fraction=0, fold_count=2, guard_bp=0)
        self.assertEqual(actual["components"], expected)
        self.assertEqual(actual, REPLICATION.reconstruct_grouping(list(reversed(rows)), minimum_components=1,
                         largest_component_fraction=1, minimum_training_fraction=0, fold_count=2, guard_bp=0))

    def test_guard_boundary_is_inclusive_and_uses_whole_components(self):
        rng = random.Random(11)
        rows = [{"id": identifier, "sequence": "".join(rng.choice("ACGT") for _ in range(200)),
                 "chromosome": "chr1", "start": start, "end": start + 200, "parent_id": identifier}
                for identifier, start in [("a", 999800), ("b", 1010000), ("c", 20000000)]]
        options = {"minimum_components": 1, "largest_component_fraction": 1,
                   "minimum_training_fraction": 0, "fold_count": 2, "guard_bp": 10000}
        actual = REPLICATION.reconstruct_grouping(rows, **options)
        fold = next(fold for fold in actual["folds"] if "a" in fold["test_ids"])
        self.assertIn("b", fold["guard_ids"])
        rows[1]["start"] += 1
        rows[1]["end"] += 1
        actual = REPLICATION.reconstruct_grouping(rows, **options)
        fold = next(fold for fold in actual["folds"] if "a" in fold["test_ids"])
        self.assertIn("b", fold["train_ids"])

    def test_exact_duplicate_source_ids_still_count_against_candidate_budget(self):
        rows = [{"id": str(i), "sequence": "A" * 60, "chromosome": "chr" + str(i),
                 "start": 0, "end": 60, "parent_id": str(i)} for i in range(2)]
        with self.assertRaisesRegex(ValueError, "candidate budget"):
            REPLICATION.reconstruct_grouping(rows, candidate_budget=0, minimum_components=1,
                                             largest_component_fraction=1, minimum_training_fraction=0, fold_count=2)

    def test_unfrozen_or_changed_protocol_fails_before_source_loading(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "protocol.json"
            path.write_text('{"status":"frozen"}')
            with self.assertRaisesRegex(ValueError, "separately frozen"):
                REPLICATION.checked_protocol(path)

    def test_exact_source_sum_audit_distinguishes_decimal_identity_from_float_roundoff(self):
        rows = [{"id": "zero", "score_source_row": 2,
                 "source_literals": {"a": "0.1", "b": "0.2", "c": "-0.3"},
                 "observations": {"a": 0.1, "b": 0.2, "c": -0.3}},
                {"id": "nonzero", "score_source_row": 3,
                 "source_literals": {"a": "1", "b": "2", "c": "-2.99999"},
                 "observations": {"a": 1, "b": 2, "c": -2.99999}}]
        before = json.dumps(rows, sort_keys=True)
        audit = REPLICATION.published_ess_semantic_audit(rows, ["a", "b", "c"])
        self.assertEqual(audit["count_abs_sum_gt1e-12"], 1)
        self.assertEqual(audit["count_abs_sum_gt1e-8"], 1)
        self.assertEqual(audit["count_abs_sum_gt1e-6"], 1)
        self.assertEqual(audit["maximum_source_row"]["sum_exact_numerator"], "1")
        self.assertEqual(audit["maximum_source_row"]["sum_exact_denominator"], "100000")
        self.assertEqual(before, json.dumps(rows, sort_keys=True))

    def test_comparison_checks_shapes_nonfinite_and_exact_bounds_without_rank_hacks(self):
        comparison = REPLICATION.NumericalComparison()
        comparison.numeric("within", [0], [1e-7], 1e-7)
        comparison.numeric("beyond", [0], [1.01e-7], 1e-7)
        comparison.numeric("shape", [0], [[0]], 1e-7)
        comparison.numeric("nonfinite", [np.nan], [np.nan], 1e-7)
        report = comparison.report()
        self.assertFalse(report["matched"])
        self.assertEqual(report["failure_count"], 3)
        self.assertTrue(report["checks"][0]["matched"])

    def test_gzip_reader_uses_magic_and_commits_raw_and_decoded_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "no-gzip-extension.data"
            decoded = b'{"value":1}\n'
            raw = gzip.compress(decoded, mtime=0)
            path.write_bytes(raw)
            value, commitment = REPLICATION.read_json_payload(path)
            self.assertEqual(value, {"value": 1})
            self.assertEqual(commitment["raw_sha256"], hashlib.sha256(raw).hexdigest())
            self.assertEqual(commitment["decoded_sha256"], hashlib.sha256(decoded).hexdigest())

    def test_output_manifest_checksum_and_escape_checks_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            (root / "data").write_bytes(b"good")
            manifest = {"schema_version": 1, "files": {"data": hashlib.sha256(b"good").hexdigest()}}
            (root / "manifest.json").write_text(json.dumps(manifest))
            REPLICATION.verified_output_manifest(root, "manifest.json", 1)
            (root / "data").write_bytes(b"bad")
            with self.assertRaisesRegex(ValueError, "checksum"):
                REPLICATION.verified_output_manifest(root, "manifest.json", 1)
            manifest["files"] = {"../outside": "0" * 64}
            (root / "manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "escapes"):
                REPLICATION.verified_output_manifest(root, "manifest.json", 1)


if __name__ == "__main__":
    unittest.main()
