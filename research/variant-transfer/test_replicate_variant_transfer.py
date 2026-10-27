"""Independent transfer controls use synthetic numbers and generated workbooks."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import openpyxl

SPEC = importlib.util.spec_from_file_location("independent_variant", Path(__file__).resolve().parents[2]
                                           / "scripts/replicate_variant_transfer.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def rows():
    result = []
    for i, (component, observed, predicted) in enumerate([
        ("a", [1.0, 3.0], [0.5, 2.0]),
        ("b", [2.0, 4.0], [1.0, 2.0]),
        ("b", [0.0, -2.0], [0.0, -1.0]),
    ]):
        result.append({"id": str(i), "component": component, "locus": str(i),
            "gc_preserving": i == 0, "observed": dict(zip(M.CELLS, observed)),
            "predictions": {"I": dict(zip(M.CELLS, predicted)), "S": dict.fromkeys(M.CELLS, 0.0)}})
    return result


def workbook(path, kind, values=None, duplicate=False, header_error=False):
    book = openpyxl.Workbook()
    sheet = book.active
    identifier = "chr1:100:A:C"
    if kind == "dna":
        sheet.title = "S2.Variant Oligo Seqs"
        headers = {1: "Variant", 2: "Allele 1 Oligo", 3: "Allele 2 Oligo"}
        data = {1: identifier, 2: "A" * 200, 3: "A" * 99 + "C" + "A" * 100}
    else:
        sheet.title = "S4.MPRA results"
        headers = {1: "Variant", 40: "log2Skew in HEPG2", 41: "log2Skew in K562"}
        data = {1: identifier, 40: 1.0, 41: -2.0}
    if header_error:
        headers[1] = "Unexpected ID"
    for col, value in headers.items():
        sheet.cell(2, col, value)
    data.update(values or {})
    for col, value in data.items():
        sheet.cell(3, col, value)
        if duplicate:
            sheet.cell(4, col, value)
    book.save(path)
    book.close()
    return {"sha256": M.sha256(path), "bytes": Path(path).stat().st_size}


class IndependentTransferTests(unittest.TestCase):
    def test_exact_feature_counts_rc_and_delta(self):
        reference = "ACGT" * 50
        self.assertEqual(M.exact_counts(reference).sum(), 298)
        np.testing.assert_array_equal(M.exact_counts(reference), M.exact_counts(M.rc(reference)))
        alternate = reference[:99] + "A" + reference[100:]
        delta = M.exact_counts(alternate) - M.exact_counts(reference)
        self.assertEqual(delta[:-1].sum(), 0)
        self.assertTrue(np.any(delta))
        for bad in ["A" * 199, "a" * 200, "N" * 200]:
            with self.assertRaises(ValueError):
                M.exact_counts(bad)

    def test_zero_feature_change_is_exact(self):
        reference = "A" * 97 + "AAATT" + "C" * 98
        alternate = reference[:99] + "T" + reference[100:]
        np.testing.assert_array_equal(M.exact_counts(reference), M.exact_counts(alternate))

    def test_numpy_five_fold_prediction_and_contexts(self):
        reference = "A" * 200
        alternate = reference[:99] + "C" + reference[100:]
        data = {"x": {"delta_counts": (M.exact_counts(alternate) - M.exact_counts(reference)).tolist()}}
        slopes = {model: {cell: np.zeros((5, 33)) for cell in M.CELLS} for model in ("I", "S")}
        for model in slopes:
            for cell in M.CELLS:
                slopes[model][cell][:, -1] = np.arange(5)
        names, predictions, folds = M.predict_fixed(data, slopes)
        self.assertEqual(names, ["x"])
        for model in slopes:
            for cell in M.CELLS:
                self.assertAlmostEqual(predictions[model][cell][0], 0.01)
                np.testing.assert_array_equal(folds[model][cell][0], np.arange(5) * 0.005)

    def test_macro_component_weighting_has_hand_computed_value(self):
        result = M.evaluate(rows(), replicates=8, seed=0)
        macro = result["observed"]["equal_component"]
        self.assertEqual(macro["losses"]["paired"]["zero"]["mse"], 5.5)
        self.assertEqual(macro["losses"]["paired"]["I"]["mse"], 1.0625)
        self.assertEqual(macro["gains_vs_zero"]["paired"]["I"]["mse"], 4.4375)
        self.assertAlmostEqual(result["observed"]["equal_variant"]["losses"]["paired"]["I"]["mse"],
                               (0.625 + 2.5 + 0.5) / 3)
        self.assertFalse(result["primary"]["minimum_support_passed"])

    def test_macro_is_invariant_to_whole_component_replication(self):
        original = rows()
        replicated = original + [copy.deepcopy(row) | {"id": "z" + row["id"]} for row in original[1:]]
        a = M.evaluate(original, replicates=3)
        b = M.evaluate(replicated, replicates=3)
        self.assertEqual(a["observed"]["equal_component"], b["observed"]["equal_component"])
        self.assertNotEqual(a["observed"]["equal_variant"], b["observed"]["equal_variant"])

    def test_rng_known_splitmix64_values_and_draw_resource_bound(self):
        indices = M.splitmix_indices(13, replicates=1, seed=0)
        expected = [0xE220A8397B1DCDAF, 0x6E789E6AA1B965F4, 0x06C45D188009454F]
        self.assertEqual(indices[0, :3].tolist(), [value % 13 for value in expected])
        with self.assertRaises(ValueError):
            M.splitmix_indices(25001, replicates=2000)

    def test_type7_interpolation(self):
        self.assertEqual(M.type7([0, 10, 20, 30], [0.25, 0.5, 0.75]), [7.5, 15.0, 22.5])

    def test_all_draws_and_stats_are_retained(self):
        result = M.evaluate(rows(), replicates=12, seed=5)
        bootstrap = result["resampling"]
        self.assertEqual(len(bootstrap["resampling_indices"]), 12)
        self.assertEqual(len(bootstrap["sampled_statistics"]), 12)
        for draw, record in zip(bootstrap["resampling_indices"], bootstrap["sampled_statistics"]):
            losses = [result["per_component"][i]["losses"]["paired"]["I"]["mse"] for i in draw]
            self.assertEqual(np.mean(losses), record["equal_component"]["losses"]["paired"]["I"]["mse"])

    def test_association_retains_zeros_and_constant_vector_is_nonestimable(self):
        same = M.association([-1.0, 0.0, 1.0], [-2.0, 0.0, 2.0])
        self.assertEqual(same["sign_agreement_count"], 3)
        self.assertEqual(same["exact_zero_predictions"], 1)
        constant = M.association([0.1] * 7, [0.0] * 7)
        self.assertFalse(constant["pearson"]["estimable"])
        self.assertEqual(constant["sign_agreement_count"], 0)

    def test_gc_subgroups_preserve_and_list_original_components(self):
        result = M.evaluate(rows(), replicates=3)
        a = result["gc_diagnostics_descriptive_only"]["GC_preserving"]
        self.assertEqual(a["represented_component_ids"], ["a"])
        self.assertEqual(a["empty_original_component_ids"], ["b"])
        self.assertEqual(a["variants"], 1)

    def test_source_dna_original_xlsx_and_snv_validation(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "S2.xlsx"
            commitment = workbook(path, "dna")
            found, report = M.read_s2_dna(path, commitment, {"chr1:100:A:C"})
            self.assertEqual(report["required_ids"], 1)
            self.assertFalse(found["chr1:100:A:C"]["gc_preserving"])
            self.assertEqual(found["chr1:100:A:C"]["delta_counts"][-1], 1)
            commitment = workbook(path, "dna", values={3: "A" * 98 + "C" + "A" * 101})
            with self.assertRaises(ValueError):
                M.read_s2_dna(path, commitment, {"chr1:100:A:C"})

    def test_source_header_hash_and_duplicate_ids_are_strict(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "S2.xlsx"
            commitment = workbook(path, "dna")
            with self.assertRaises(ValueError):
                M.read_s2_dna(path, commitment | {"sha256": "0" * 64}, {"chr1:100:A:C"})
            commitment = workbook(path, "dna", header_error=True)
            with self.assertRaises(ValueError):
                M.read_s2_dna(path, commitment, {"chr1:100:A:C"})
            commitment = workbook(path, "dna", duplicate=True)
            with self.assertRaises(ValueError):
                M.read_s2_dna(path, commitment, {"chr1:100:A:C"})

    def test_source_effects_exact_bits_and_bound_missing_policy(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "S4.xlsx"
            commitment = workbook(path, "effects")
            found, report = M.read_s4_effects(path, commitment, {"chr1:100:A:C"}, [])
            self.assertEqual(found["chr1:100:A:C"]["observed"]["K562"], -2.0)
            self.assertEqual(found["chr1:100:A:C"]["observed_float64_little_endian"]["HepG2"],
                             M.float_bits(1.0))
            commitment = workbook(path, "effects", values={40: "NA"})
            with self.assertRaises(ValueError):
                M.read_s4_effects(path, commitment, {"chr1:100:A:C"}, [])
            found, report = M.read_s4_effects(path, commitment, {"chr1:100:A:C"}, ["NA"])
            self.assertEqual(found, {})
            self.assertEqual(report["incomplete"]["chr1:100:A:C"]["missing_contexts"], ["HepG2"])

    def test_source_formula_error_bool_and_numeric_strings_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "S4.xlsx"
            for invalid in ["=1+1", "#DIV/0!", True, "1.5", None]:
                commitment = workbook(path, "effects", values={40: invalid})
                with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                    M.read_s4_effects(path, commitment, {"chr1:100:A:C"}, [])

    def test_missing_source_id_requires_separately_bound_absence_rule(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "S4.xlsx"
            commitment = workbook(path, "effects")
            with self.assertRaises(ValueError):
                M.read_s4_effects(path, commitment, {"absent"}, [])
            found, report = M.read_s4_effects(path, commitment, {"absent"}, [], absence_allowed=True)
            self.assertEqual(report["absent_ids"], ["absent"])

    def test_comparison_tolerance_and_exact_metadata_failures(self):
        original = M.evaluate(rows(), replicates=3)
        self.assertTrue(M.compare_trees(original, copy.deepcopy(original))["matched"])
        changed = copy.deepcopy(original)
        changed["primary"]["observed_improvement"] += 1e-11
        self.assertTrue(M.compare_trees(original, changed)["matched"])
        changed["primary"]["observed_improvement"] += 1e-8
        self.assertFalse(M.compare_trees(original, changed)["matched"])
        changed = copy.deepcopy(original)
        changed["per_variant"][0]["observed_float64_little_endian"]["K562"] = M.float_bits(1.00000000001)
        self.assertFalse(M.compare_trees(original, changed)["observed_float_bits_and_integer_indices_exact"])
        changed = copy.deepcopy(original)
        changed["resampling"]["resampling_indices"][0][0] = 1 - changed["resampling"]["resampling_indices"][0][0]
        self.assertFalse(M.compare_trees(original, changed)["matched"])
        changed = copy.deepcopy(original)
        changed["variant_count"] = float(changed["variant_count"])
        self.assertFalse(M.compare_trees(original, changed)["matched"])


if __name__ == "__main__":
    unittest.main()
