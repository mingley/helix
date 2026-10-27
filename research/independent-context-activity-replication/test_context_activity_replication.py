"""New v3 controls: literal activity binding and common-component prediction.

All workbooks, labels and model targets in this file are synthetic.
"""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
import openpyxl

PATH = Path(__file__).resolve().parents[2] / "scripts/replicate_context_activity.py"
SPEC = importlib.util.spec_from_file_location("independent_context_activity_v3", PATH)
V3 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = V3
SPEC.loader.exec_module(V3)

COLUMNS = {c: c + " [log2(rna/dna)]" for c in V3.CONTEXTS}


class ActivitySourceTests(unittest.TestCase):
    def workbook(self, directory, rows, headers=None):
        path = Path(directory) / "source.xlsx"
        workbook = openpyxl.Workbook()
        sheet = workbook.active
        sheet.title = "all_cell_types_summary"
        sheet.append(headers or ["name", *COLUMNS.values(), "ESS", "prediction"])
        for row in rows:
            sheet.append(row)
        workbook.save(path)
        workbook.close()
        return path

    def read(self, path, required=None, **kwargs):
        return V3.read_activity_cells(
            path, sheet_name="all_cell_types_summary", id_header="name",
            columns=COLUMNS, required_ids=required or {"one"}, **kwargs)

    def test_activity_literals_remain_uncentered_and_predictions_are_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = self.workbook(tmp, [["one", 10.25, -2.5, 99, -9999, 8888]])
            values, report = self.read(source)
            self.assertEqual(values["one"]["observations"],
                             {"HepG2": 10.25, "K562": -2.5, "WTC11": 99.0})
            self.assertEqual(report["activity_transform"], "none; source float values preserved unchanged")
            self.assertNotEqual(sum(values["one"]["observations"].values()), 0)

    def test_missing_required_element_blocks_instead_of_selecting_subset(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = self.workbook(tmp, [["one", 1, 2, 3]])
            with self.assertRaisesRegex(ValueError, "coverage"):
                self.read(source, {"one", "two"})

    def test_duplicate_required_id_blocks_even_when_values_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = self.workbook(tmp, [["one", 1, 2, 3], ["one", 1, 2, 3]])
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                self.read(source)

    def test_no_case_or_suffix_alias_for_source_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = self.workbook(tmp, [["ONE", 1, 2, 3], ["one_F", 1, 2, 3]])
            with self.assertRaisesRegex(ValueError, "coverage"):
                self.read(source)

    def test_headers_are_exact_not_substring_or_favorable_aliases(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = self.workbook(tmp, [["one", 1, 2, 3]],
                                   ["name", "prediction " + COLUMNS["HepG2"],
                                    COLUMNS["K562"], COLUMNS["WTC11"]])
            with self.assertRaisesRegex(ValueError, "exact header"):
                self.read(source)

    def test_duplicated_header_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = self.workbook(tmp, [["one", 1, 2, 3, 4]],
                                   ["name", *COLUMNS.values(), COLUMNS["HepG2"]])
            with self.assertRaisesRegex(ValueError, "duplicated"):
                self.read(source)

    def test_formula_boolean_error_and_numeric_text_are_not_observations(self):
        for invalid in ["=1+2", True, "#N/A", "3.0", None]:
            with self.subTest(invalid=invalid), tempfile.TemporaryDirectory() as tmp:
                source = self.workbook(tmp, [["one", invalid, 2, 3]])
                with self.assertRaisesRegex(ValueError, "coverage"):
                    self.read(source)

    def test_nonfinite_literals_are_rejected_without_xlsx_coercion(self):
        for invalid in [float("nan"), float("inf"), float("-inf")]:
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, "finite"):
                V3.literal_numeric(SimpleNamespace(data_type="n", value=invalid))

    def test_extraneous_source_ids_are_audited_without_changing_fixed_panel(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = self.workbook(tmp, [["not_selected", 1000, 2000, 3000], ["one", 1, 2, 3]])
            values, report = self.read(source)
            self.assertEqual(set(values), {"one"})
            self.assertEqual(report["extra_source_IDs"], ["not_selected"])
            self.assertEqual(report["source_nonempty_ID_rows"], 2)

    def test_every_central_gate_is_required_before_source_file_access(self):
        names = ["archive_integrity_verified", "empirical_activity_schema_bound", "coverage_complete",
                 "design_identity_verified", "replicate_aggregation_status_recorded"]
        for failed in names:
            with self.subTest(failed=failed), tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                manifest = {"schema_version": 1, "protocol_sha256": V3.FROZEN_V3_SHA256,
                            "endpoint": V3.ENDPOINT, "gates": {key: key != failed for key in names}}
                (directory / "manifest.json").write_text(json.dumps(manifest))
                # There is deliberately no workbook: the gate must fail before opening it.
                with self.assertRaisesRegex(ValueError, "central activity source gate"):
                    V3.verified_dataset_gate(directory, {}, V3.FROZEN_V3_SHA256)

    def test_aggregation_uncertainty_must_be_recorded_even_when_gate_flags_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest = {"schema_version": 1, "protocol_sha256": V3.FROZEN_V3_SHA256,
                        "endpoint": V3.ENDPOINT, "gates": {key: True for key in [
                            "archive_integrity_verified", "empirical_activity_schema_bound", "coverage_complete",
                            "design_identity_verified", "replicate_aggregation_status_recorded"]}}
            (directory / "manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "aggregation uncertainty"):
                V3.verified_dataset_gate(directory, {}, V3.FROZEN_V3_SHA256)

    def test_old_specificity_endpoint_cannot_pass_activity_manifest_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "manifest.json").write_text(json.dumps({
                "schema_version": 1, "protocol_sha256": V3.FROZEN_V3_SHA256,
                "endpoint": "published_element_specificity_score"}))
            with self.assertRaisesRegex(ValueError, "protocol/endpoint/schema"):
                V3.verified_dataset_gate(directory, {}, V3.FROZEN_V3_SHA256)

    def test_additional_provenance_receipts_are_verified_without_rejecting_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            names = ["activity.json.gz", "source-binding.json", "source-export.json.gz",
                     "source-ledger.jsonl.gz", "recovery-summary.json", "archive-validation.json",
                     "extraction-receipt.json", "builder.py"]
            for name in names:
                (directory / name).write_bytes(b"synthetic receipt")
            source = directory / "source.xlsx"
            source.write_bytes(b"synthetic byte commitment, not an outcome workbook")
            source_hash = V3.BASE.file_hashes(source)
            protocol = {"source": {"jointdata_expected_sha256": source_hash["sha256"],
                                   "expected_bytes": source_hash["bytes"]}}
            manifest = {"schema_version": 1, "protocol_sha256": V3.FROZEN_V3_SHA256,
                        "endpoint": V3.ENDPOINT, "gates": {key: True for key in [
                            "archive_integrity_verified", "empirical_activity_schema_bound", "coverage_complete",
                            "design_identity_verified", "replicate_aggregation_status_recorded"]},
                        "aggregation_status": {"scope": "synthetic"}, "schema_binding": {"scope": "synthetic"},
                        "files": {name: V3.BASE.file_hashes(directory / name)["sha256"] for name in names},
                        "source": {"raw_path": "source.xlsx", "sha256": source_hash["sha256"], "bytes": source_hash["bytes"]}}
            (directory / "manifest.json").write_text(json.dumps(manifest))
            V3.verified_dataset_gate(directory, protocol, V3.FROZEN_V3_SHA256)
            (directory / "builder.py").write_bytes(b"changed extra provenance")
            with self.assertRaisesRegex(ValueError, "commitment mismatch"):
                V3.verified_dataset_gate(directory, protocol, V3.FROZEN_V3_SHA256)


class CommonActivityControls(unittest.TestCase):
    def evaluation(self, target, models, query, contexts):
        predictions = []
        for index, row in enumerate(target):
            for j, context in enumerate(contexts):
                predictions.append({
                    "sequence_id": str(index), "group": "g" + str(index),
                    "fold_index_0based": 0, "context": context, "observed": float(row[j]),
                    "predictions": {name: float(model.predict([query[index]], [context])[0])
                                    for name, model in models.items()},
                })
        return {"contexts": contexts, "predictions": predictions}

    def test_common_activity_gain_has_hand_computed_ridge_value(self):
        # Training x variance 2.5, standardized z'z=4. Alpha10 shrinks slope5
        # to 5*4/(4+10)=10/7, with unpenalized mean2.
        contexts = ["a", "b", "c"]
        x = np.array([[-2.], [-1.], [1.], [2.]])
        offsets = np.array([-3., 1., 2.])
        y = 2 + 5 * x + offsets[None, :]
        models = V3.BASE.fit_balanced_panel(x, y, contexts, alpha=10)
        query = np.array([[-1.5], [1.5]])
        target = 2 + 5 * query + offsets[None, :]
        result = V3.common_activity_diagnostic(self.evaluation(target, models, query, contexts))
        expected_global = 7.5 ** 2
        expected_shared = (75 / 14) ** 2
        self.assertAlmostEqual(result["global_mean_mse"], expected_global, places=12)
        self.assertAlmostEqual(result["shared_dna_mse"], expected_shared, places=12)
        self.assertAlmostEqual(result["mean_activity_mse_improvement"],
                               expected_global - expected_shared, places=12)
        self.assertAlmostEqual(result["primary_cell_mse_improvement"],
                               result["mean_activity_mse_improvement"], places=12)
        self.assertLess(result["maximum_absolute_per_element_identity_residual"], 1e-12)

    def test_row_centering_erases_shared_activity_that_v3_targets(self):
        contexts = ["a", "b", "c"]
        x = np.array([[-2.], [-1.], [1.], [2.]])
        offsets = np.array([-3., 1., 2.])
        raw = 2 + 5 * x + offsets[None, :]
        centered = raw - raw.mean(axis=1, keepdims=True)
        models = V3.BASE.fit_balanced_panel(x, centered, contexts, alpha=10)
        for context in contexts:
            np.testing.assert_array_equal(
                models["context_mean"].predict(x, [context] * len(x)),
                models["additive_context_dna"].predict(x, [context] * len(x)))
        np.testing.assert_array_equal(models["shared_dna"].shared.standardized_coefficients, [0.])

    def test_missing_context_prevents_mean_activity_identity_audit(self):
        evaluation = {"contexts": ["a", "b", "c"], "predictions": []}
        with self.assertRaisesRegex(ValueError, "empty"):
            V3.common_activity_diagnostic(evaluation)

    def test_v3_primary_bootstrap_is_C_minus_A_in_retained_order(self):
        losses = [[9, 8, 7, 6, 5]]
        bootstrap = V3.BASE.paired_component_bootstrap(
            ["only_group"], losses, losses, V3.COMPARISONS, replicates=2, seed=1)
        primary = bootstrap["equal_sequences"]["comparisons"][1]
        self.assertEqual(primary["comparison"], {"baseline": "context_mean", "candidate": "additive_context_dna"})
        self.assertEqual(primary["observed_mse_improvement"], 2)
        self.assertEqual(primary["mse_interval"], [2, 2])


if __name__ == "__main__":
    unittest.main()
