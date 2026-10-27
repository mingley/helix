"""Controls for literal measured-activity binding, independent of fitted models."""
import copy
import itertools
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from build_context_activity import activity_records, finite_source, validate_recovery_variants


class BindingControls(unittest.TestCase):
    def fixture(self):
        return {"tables": {"all_cell_types_summary": {
            "id_header": "name",
            "headers": {"A": "name", "B": "HepG2 [log2(rna/dna)]",
                        "C": "K562 [log2(rna/dna)]", "D": "WTC11 [log2(rna/dna)]"},
            "rows": [{"id": name, "source_row": i + 2, "values": {
                "HepG2 [log2(rna/dna)]": "2.125e0", "K562 [log2(rna/dna)]": "-4.5",
                "WTC11 [log2(rna/dna)]": "0.000125"}}
                for i, name in enumerate(["b", "a", "outside"]) ]}}}

    def test_literal_values_and_fixed_id_order(self):
        records, coverage = activity_records(self.fixture(), {"a", "b"})
        self.assertEqual([r["id"] for r in records], ["a", "b"])
        self.assertEqual(records[0]["observations"], {
            "HepG2": 2.125, "K562": -4.5, "WTC11": 0.000125})
        self.assertEqual(coverage["source_ids_outside_fixed_cohort"], 1)

    def test_missing_required_id_blocks_whole_cohort(self):
        with self.assertRaisesRegex(ValueError, "incomplete"):
            activity_records(self.fixture(), {"a", "b", "absent"})

    def test_duplicate_source_id_even_outside_cohort_rejected(self):
        fixture = self.fixture()
        fixture["tables"]["all_cell_types_summary"]["rows"].append(
            copy.deepcopy(fixture["tables"]["all_cell_types_summary"]["rows"][-1]))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            activity_records(fixture, {"a", "b"})

    def test_specificity_cannot_replace_activity_header(self):
        fixture = self.fixture()
        fixture["tables"]["all_cell_types_summary"]["headers"]["B"] = "Element Specificity Score"
        with self.assertRaisesRegex(ValueError, "headers"):
            activity_records(fixture, {"a", "b"})

    def test_id_alias_is_not_guessed(self):
        fixture = self.fixture()
        fixture["tables"]["all_cell_types_summary"]["id_header"] = "Sequence ID"
        with self.assertRaisesRegex(ValueError, "ID header"):
            activity_records(fixture, {"a", "b"})

    def test_missing_context_is_not_imputed(self):
        fixture = self.fixture()
        del fixture["tables"]["all_cell_types_summary"]["rows"][0]["values"]["K562 [log2(rna/dna)]"]
        with self.assertRaisesRegex(ValueError, "lacks"):
            activity_records(fixture, {"a", "b"})

    def test_nonfinite_and_missing_tokens_rejected(self):
        for value in ["NaN", "Infinity", "NA", "", "1e9999", None, True, 2.0]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                finite_source(value)

    def test_numeric_zero_and_negative_are_not_filtered(self):
        fixture = self.fixture()
        for row in fixture["tables"]["all_cell_types_summary"]["rows"]:
            row["values"]["HepG2 [log2(rna/dna)]"] = "0"
            row["values"]["K562 [log2(rna/dna)]"] = "-100"
        records, _ = activity_records(fixture, {"a", "b"})
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["observations"]["K562"], -100)

    def grid_fixture(self):
        return {"variants": [{"id": str(i), "activity_base": base,
                              "pre_zscore_population": population,
                              "pre_zscore_ddof": ddof, "post_ESS_step": post}
                             for i, (base, population, ddof, post) in enumerate(itertools.product(
            ("published_activity_aggregate", "mean_three_published_replicates",
             "mean_three_after_explicit_replicate_median_subtraction"),
            ("complete_three_context_source", "all_finite_per_context"), (0, 1),
            ("none", "column_zscore_ddof0", "column_zscore_ddof1")))]}

    def test_complete_grid_and_duplicate_or_dropped_route(self):
        validate_recovery_variants(self.grid_fixture())
        fixture = self.grid_fixture()
        fixture["variants"][-1] = copy.deepcopy(fixture["variants"][0])
        with self.assertRaisesRegex(ValueError, "36-route"):
            validate_recovery_variants(fixture)
        fixture = self.grid_fixture()
        fixture["variants"].pop()
        with self.assertRaisesRegex(ValueError, "36-route"):
            validate_recovery_variants(fixture)

    def test_false_metadata_cannot_substitute_a_route(self):
        fixture = self.grid_fixture()
        fixture["variants"][0]["pre_zscore_ddof"] = False
        with self.assertRaisesRegex(ValueError, "36-route"):
            validate_recovery_variants(fixture)


if __name__ == "__main__":
    unittest.main()
