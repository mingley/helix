#!/usr/bin/env python3
"""Hand-computed synthetic cases; no recovered source values or downloads."""
import importlib.util
import itertools
import json
import math
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("activity_audit", REPOSITORY / "scripts/audit_recovered_context.py")
AUDIT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = AUDIT
SPEC.loader.exec_module(AUDIT)


def endpoint(name="fixture", values=None, ids=None, raw=None):
    values = np.asarray(values if values is not None else [[-1, 0, 1], [0, 1, -1], [1, -1, 0]], dtype=float)
    return AUDIT.Endpoint(name, ids or list("abc")[:len(values)], values, "fixture", raw=raw)


def table(headers, rows, id_header="name"):
    return {"id_header": id_header, "headers": dict(zip("ABCDEFGHIJKLMNOPQRSTUVWXYZ", [id_header] + headers)),
            "selected_value_headers": headers,
            "rows": [{"id": identifier, "source_row": index + 2, "values": dict(zip(headers, values))}
                     for index, (identifier, values) in enumerate(rows)]}


class ArithmeticTests(unittest.TestCase):
    def test_unbound_schema_blocks_before_source_values_are_opened(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binding = root / "binding.json"
            binding.write_text(json.dumps({"protocol_sha256": AUDIT.PROTOCOL_SHA256,
                "status": "not_bound", "semantic_binding_complete": False}))
            arguments = SimpleNamespace(out=root/"new-output", protocol=REPOSITORY/"experiments/context-endpoint-recovery-v1.json",
                                        binding=binding, source=root/"nonexistent-values.gz")
            with self.assertRaisesRegex(ValueError, "schema gate"):
                AUDIT.run(arguments)
            self.assertFalse(arguments.out.exists())

    def test_existing_attempt_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            marker = root / "previous-result"
            marker.write_text("retain me")
            with self.assertRaises(FileExistsError):
                AUDIT.run(SimpleNamespace(out=root))
            self.assertEqual(marker.read_text(), "retain me")

    def test_hand_computed_zscores_and_row_centering(self):
        original = endpoint()
        # Every column has mean0, population variance2/3, sample variance1;
        # every row already sums0. Thus all derived values are known by hand.
        for ddof, factor in ((0, math.sqrt(1.5)), (1, 1.0)):
            record, result = AUDIT.pipeline(original, AUDIT.POPULATIONS[0], ddof, "none")
            self.assertTrue(record["estimable"])
            np.testing.assert_allclose(result.values, original.values * factor, atol=1e-15)
            np.testing.assert_allclose(result.values.sum(axis=1), 0, atol=1e-15)

    def test_sample_population_equivalence_after_post_standardization(self):
        original = endpoint()
        for ddof, post in itertools.product((0, 1), AUDIT.POST[1:]):
            _, result = AUDIT.pipeline(original, AUDIT.POPULATIONS[0], ddof, post)
            expected = original.values * (math.sqrt(1.5) if post.endswith("ddof0") else 1)
            np.testing.assert_allclose(result.values, expected, atol=1e-15)

    def test_missing_context_changes_only_prespecified_pre_population(self):
        original = endpoint(values=[[-1, 0, 1], [0, 1, -1], [1, -1, 0], [10, math.nan, math.nan]], ids=list("abcd"))
        complete, result1 = AUDIT.pipeline(original, AUDIT.POPULATIONS[0], 0, "none")
        percell, result2 = AUDIT.pipeline(original, AUDIT.POPULATIONS[1], 0, "none")
        self.assertEqual(result1.ids, list("abc"))
        self.assertEqual(result2.ids, list("abc"))
        self.assertEqual(complete["pre_parameters"][0]["count"], 3)
        self.assertEqual(percell["pre_parameters"][0]["count"], 4)
        self.assertEqual(percell["pre_parameters"][0]["mean"], 2.5)
        self.assertFalse(np.allclose(result1.values, result2.values))

    def test_all_36_unique_tuples_even_when_nonestimable(self):
        bases = [endpoint(name=base) for base in AUDIT.BASES]
        records = [AUDIT.pipeline(*definition)[0]
                   for definition in itertools.product(bases, AUDIT.POPULATIONS, (0, 1), AUDIT.POST)]
        self.assertEqual(len(records), 36)
        self.assertEqual(len({r["id"] for r in records}), 36)
        _, bad = AUDIT.pipeline(endpoint(values=[[1, 1, 1]], ids=["a"]), AUDIT.POPULATIONS[0], 1, "none")
        self.assertIsNone(bad)
        bad, _ = AUDIT.pipeline(endpoint(values=[[1, 1, 1], [1, 1, 1]], ids=["a", "b"]), AUDIT.POPULATIONS[0], 0, "none")
        self.assertIn("nonpositive", bad["nonestimable_reason"])

    def test_exact_decimal_constant_is_not_a_tiny_positive_sd(self):
        # A binary-float mean of seven copies of0.1 can differ by one ULP;
        # np.std then gives a tiny positive value for an exactly constant input.
        source = endpoint(values=[[0.1, i, 6-i] for i in range(7)], ids=[str(i) for i in range(7)])
        for ddof in (0, 1):
            report, derived = AUDIT.pipeline(source, AUDIT.POPULATIONS[0], ddof, "none")
            self.assertFalse(report["estimable"])
            self.assertIsNone(derived)
            self.assertIn("exactly constant", report["nonestimable_reason"])
        report, derived = AUDIT.code_hypothesis(source, source, source.ids, source.ids)
        self.assertFalse(report["estimable"])
        self.assertIsNone(derived)
        self.assertIn("exactly constant", report["nonestimable_reason"])

    def test_exact_decimal_is_only_direct_source_comparison(self):
        raw1 = {"a": ["1.000", "2", "3"]}
        raw2 = {"a": ["1", "2.0", "3.000"]}
        direct1 = endpoint(values=[[1, 2, 3]], ids=["a"], raw=raw1)
        direct2 = endpoint(values=[[1, 2, 3]], ids=["a"], raw=raw2)
        report = AUDIT.compare(direct1, direct2, ["a"], ["a"])["cohorts"]["full_source_population"]
        self.assertEqual(report["complete_three_context_pooled"]["exact_decimal"]["equal_count"], 3)
        derived = endpoint(values=[[1, 2, 3]], ids=["a"])
        report = AUDIT.compare(derived, direct2, ["a"], ["a"])["cohorts"]["full_source_population"]
        self.assertFalse(report["complete_three_context_pooled"]["exact_decimal"]["applicable"])
        self.assertIsNone(report["complete_three_context_pooled"]["exact_decimal"]["equal_count"])

    def test_fixed_tolerance_uses_unchanged_target_and_type7_quantile(self):
        result = AUDIT.metric(np.array([0.0, 1.0, 2.0, 3.0]), np.zeros(4))
        self.assertEqual(result["within_fixed_tolerance_count"], 1)
        self.assertEqual(result["absolute_error_type7_quantiles"]["0.5"], 1.5)
        self.assertAlmostEqual(result["absolute_error_type7_quantiles"]["0.95"], 2.85)
        self.assertEqual(AUDIT.metric(np.array([1e8 + 0.5]), np.array([1e8]))["within_fixed_tolerance_count"], 1)

    def test_admitted_intersection_does_not_refit_normalization(self):
        original = endpoint()
        _, result = AUDIT.pipeline(original, AUDIT.POPULATIONS[0], 0, "none")
        target = endpoint(values=result.values.copy())
        comparison = AUDIT.compare(result, target, list("abc"), ["a"])
        selected = comparison["cohorts"]["fixed_old_admitted_intersection"]
        self.assertEqual(selected["complete_three_context_paired_ids"], 1)
        self.assertEqual(selected["complete_three_context_pooled"]["maximum_absolute_error"], 0)
        self.assertAlmostEqual(selected["candidate_moments_on_complete_comparison"]["HepG2"]["mean"], -math.sqrt(1.5))

    def test_missing_and_exact_ids_have_explicit_ledgers(self):
        headers = list(AUDIT.CONTEXTS)
        tables = {"scores": table(headers, [("a", ["1.000", "NA", None]), (None, ["2", "3", "4"]), ("a_F", ["2", "3", "4"])])}
        role = {"sheet": "scores", "id_header": "name", "columns": dict(zip(AUDIT.CONTEXTS, headers))}
        result = AUDIT.direct_endpoint("fixture", tables, role, "source")
        self.assertEqual(result.ids, ["a", "a_F"])
        self.assertEqual(len(result.ledger["rows"]), 3)
        self.assertEqual(result.ledger["rows"][1]["identifier_status"], "missing_or_invalid")
        self.assertTrue(math.isnan(result.values[0, 1]))
        self.assertEqual(result.raw["a"], ["1.000", "NA", None])

    def test_duplicate_source_ids_block_instead_of_selecting(self):
        headers = list(AUDIT.CONTEXTS)
        tables = {"scores": table(headers, [("a", ["1", "2", "3"]), ("a", ["4", "5", "6"])])}
        role = {"sheet": "scores", "id_header": "name", "columns": dict(zip(AUDIT.CONTEXTS, headers))}
        source = AUDIT.direct_endpoint("fixture", tables, role, "source")
        self.assertIsNotNone(source.blocked)
        self.assertFalse(AUDIT.compare(source, source, ["a"], ["a"])["estimable"])

    def test_long_replicates_medians_and_missing_key(self):
        tables, roles = {}, {}
        for c, context in enumerate(AUDIT.CONTEXTS):
            rows = []
            for identifier, shift in (("a", 0), ("b", 4)):
                for replicate, offset in ((1, 10), (2, 20), (3, 30)):
                    rows.append((identifier, [str(replicate), str(offset + shift + c)]))
            tables[context] = table(["replicate", "activity"], rows)
            roles[context] = {"sheet": context, "id_header": "name", "replicate_header": "replicate", "activity_header": "activity"}
        plain, centered = AUDIT.replicate_bases(tables, roles)
        self.assertIsNone(plain.blocked)
        np.testing.assert_allclose(plain.values, [[20, 21, 22], [24, 25, 26]])
        np.testing.assert_allclose(centered.values, [[-2, -2, -2], [2, 2, 2]])
        # Repeated IDs across distinct replicates are valid. Repeated keys block.
        tables["HepG2"]["rows"].append(tables["HepG2"]["rows"][0].copy())
        plain, _ = AUDIT.replicate_bases(tables, roles)
        self.assertIn("duplicate", plain.blocked)
        tables["HepG2"]["rows"].pop()
        tables["HepG2"]["rows"].pop(0)
        plain, _ = AUDIT.replicate_bases(tables, roles)
        self.assertTrue(math.isnan(plain.values[0, 0]))
        self.assertIn({"id": "a", "context": "HepG2", "replicate": 1, "reason": "unobserved_empirical_key"}, plain.ledger["missing_required_empirical_keys"])

    def test_code_hypothesis_uses_entire_exact_s12_subset(self):
        source = endpoint(values=[[-2, 0, 2], [0, 3, -3], [2, -3, 1], [100, 100, 100]], ids=list("abcd"))
        selected = source.values[:3]
        target = endpoint(values=(selected-selected.mean(axis=0))/selected.std(axis=0, ddof=1))
        report, result = AUDIT.code_hypothesis(source, target, list("abc"), ["a"])
        self.assertTrue(report["estimable"])
        self.assertTrue(report["fixed_1e-12_test_by_cohort"]["full_s12_source_id_set"]["matched"])
        self.assertEqual(report["normalization_population"]["count"], 3)
        self.assertEqual(report["normalization_parameters"]["HepG2"]["mean"], 0)
        self.assertFalse(report["historical_pipeline_certified"])
        incomplete = endpoint(values=[[0, 0, 0], [1, 1, 1]], ids=["a", "missing"])
        report, result = AUDIT.code_hypothesis(source, incomplete, list("abc"), ["a"])
        self.assertFalse(report["estimable"])
        self.assertIsNone(result)
        self.assertEqual(report["missing_or_incomplete_S11_ids"], ["missing"])


if __name__ == "__main__":
    unittest.main()
