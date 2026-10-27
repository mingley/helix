"""Synthetic-only controls for sealed variant-transfer evaluation."""
import copy
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import types
import unittest

SPEC = importlib.util.spec_from_file_location("variant_evaluate", Path(__file__).resolve().parents[2] / "scripts" / "evaluate_variant_transfer.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def row(name, group, observed, prediction, preserving=True):
    return {"id": name, "component": group, "locus": "locus-" + group,
            "gc_preserving": preserving,
            "observed": dict(zip(M.CONTEXTS, observed)),
            "predictions": {"I": dict(zip(M.CONTEXTS, prediction)),
                            "S": dict.fromkeys(M.CONTEXTS, 0.0)}}


def fixture():
    return [row("a1", "A", [2.0, 0.0], [1.0, 0.0]),
            row("a2", "A", [0.0, 2.0], [0.0, 1.0]),
            row("b1", "B", [1.0, 1.0], [3.0, 3.0], False)]


def payload(variants=None):
    return {"schema_version": 1, "protocol_sha256": "1" * 64,
            "binding_sha256": "2" * 64, "prediction_I_sha256": "3" * 64,
            "prediction_S_sha256": "4" * 64, "variants": fixture() if variants is None else variants}


def protocol():
    return {"protocol_id": "variant-transfer-v1", "status": "frozen",
            "contexts": {"primary": M.CONTEXTS}, "analysis_contract": copy.deepcopy(M.ANALYSIS_CONTRACT),
            "pretrained_model": {"source_artifacts_sha256": "a" * 64,
                                 "source_protocol_sha256": "b" * 64}}


class VariantEvaluationTests(unittest.TestCase):
    def test_hand_computed_macro_differs_from_pooled_with_unequal_components(self):
        result = M.evaluate(fixture(), replicates=8)
        macro = result["observed"]["equal_component"]
        pooled = result["observed"]["equal_variant"]
        self.assertEqual(macro["losses"]["paired"]["zero"]["mse"], 1.5)
        self.assertEqual(macro["losses"]["paired"]["I"]["mse"], 2.25)
        self.assertEqual(macro["gains_vs_zero"]["paired"]["I"]["mse"], -0.75)
        self.assertEqual(macro["gains_vs_zero"]["paired"]["I"]["mae"], -0.25)
        self.assertAlmostEqual(pooled["gains_vs_zero"]["paired"]["I"]["mse"], 0.0)
        self.assertEqual(result["variant_order"], ["a1", "a2", "b1"])
        self.assertEqual(result["resampling"]["component_order"], ["A", "B"])

    def test_every_draw_carries_whole_components_and_both_cells(self):
        result = M.evaluate(fixture(), replicates=32)
        resampling = result["resampling"]
        self.assertEqual(len(resampling["resampling_indices"]), 32)
        for indices, statistics in zip(resampling["resampling_indices"], resampling["sampled_statistics"]):
            self.assertEqual(len(indices), 2)
            expected = sum([1.5, -3.0][i] for i in indices) / 2
            self.assertAlmostEqual(statistics["equal_component"]["gains_vs_zero"]["paired"]["I"]["mse"], expected)
            for cell in M.CONTEXTS:
                self.assertAlmostEqual(statistics["equal_component"]["gains_vs_zero"][cell]["I"]["mse"], expected)
            self.assertIn("equal_variant", statistics)

    def test_splitmix_known_uint64_vector_and_exact_index_rule(self):
        rng = M.SplitMix64(0)
        self.assertEqual([rng.next_u64() for _ in range(3)],
                         [0xe220a8397b1dcdaf, 0x6e789e6aa1b965f4, 0x06c45d188009454f])
        result = M.evaluate(fixture(), replicates=4)
        independent = M.SplitMix64(M.SEED)
        expected = [[independent.next_u64() % 2 for _ in range(2)] for _ in range(4)]
        self.assertEqual(result["resampling"]["resampling_indices"], expected)

    def test_bootstrap_is_deterministic_and_type7_is_interpolated(self):
        a, b = M.evaluate(fixture(), replicates=16), M.evaluate(fixture(), replicates=16)
        self.assertEqual(a, b)
        self.assertEqual(M.quantile([0, 10, 20, 30], 0.25), 7.5)
        self.assertAlmostEqual(M.quantile([0, 10, 20, 30], 0.975), 29.25)

    def test_zero_signs_are_retained_in_denominator_and_confusion(self):
        summary = M.association([0.0, -2.0, 1.0, 0.0], [0.0, 1.0, 1.0, -1.0])
        self.assertEqual(summary["exact_sign_agreement_including_zero"], 0.5)
        self.assertEqual(summary["exact_zero_observations"], 2)
        self.assertEqual(summary["exact_zero_predictions"], 1)
        self.assertEqual(summary["sign_confusion_observed_then_predicted"]["zero"]["negative"], 1)
        constant = M.association([0.0, 1.0], [0.0, 0.0])
        self.assertFalse(constant["pearson"]["estimable"])
        self.assertIsNone(constant["pearson"]["value"])
        self.assertAlmostEqual(M.association([1, 2, 3], [3, 2, 1])["pearson"]["value"], -1)

    def test_gc_diagnostic_preserves_original_components_and_empty_list(self):
        result = M.evaluate(fixture(), replicates=4)
        preserved = result["gc_diagnostics_descriptive_only"]["GC_preserving"]
        self.assertEqual(preserved["variant_ids"], ["a1", "a2"])
        self.assertEqual(preserved["represented_component_ids"], ["A"])
        self.assertEqual(preserved["empty_original_component_ids"], ["B"])
        self.assertEqual(result["loci"][1]["variant_ids"], ["b1"])

    def test_nonzero_exact_constants_remain_nonestimable_after_rounding(self):
        for value, count in [(0.1, 11), (0.3, 37), (0.7, 35)]:
            constant = [value] * count
            varying = list(range(count))
            for observed, predicted in [(constant, varying), (varying, constant)]:
                with self.subTest(value=value, count=count):
                    result = M.association(observed, predicted)
                    self.assertFalse(result["pearson"]["estimable"])
                    self.assertIsNone(result["pearson"]["value"])

    def test_support_gate_is_separate_from_positive_performance(self):
        variants = [row(f"v{j:03}", f"g{j // 10}", [1.0, 2.0], [1.0, 2.0]) for j in range(100)]
        result = M.evaluate(variants, replicates=8)
        self.assertTrue(result["primary"]["useful_transfer_rule_passed"])
        self.assertTrue(result["primary"]["minimum_support_passed"])
        small = M.evaluate(variants[:99], replicates=8)
        self.assertTrue(small["primary"]["positive_lower_endpoint"])
        self.assertFalse(small["primary"]["useful_transfer_rule_passed"])
        self.assertFalse(M.evaluate(fixture(), replicates=8)["primary"]["positive_observed_gain"])

    def test_reject_nonfinite_bools_numeric_strings_missing_cells_duplicates(self):
        for invalid in [float("nan"), float("inf"), True, "1.0", None]:
            p = payload()
            p["variants"][0]["observed"]["K562"] = invalid
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                M.validate_variants(p)
        p = payload()
        del p["variants"][0]["observed"]["HepG2"]
        with self.assertRaises(ValueError):
            M.validate_variants(p)
        with self.assertRaises(ValueError):
            M.validate_variants(payload([fixture()[0], fixture()[0]]))
        with self.assertRaises(ValueError):
            M.validate_variants(payload([]))

    def test_loss_overflow_and_resource_overrun_fail_without_dropping(self):
        with self.assertRaises(ValueError):
            M.losses(row("huge", "A", [1e308, 1.0], [0.0, 0.0]))
        components = M.evaluate(fixture(), replicates=2)["per_component"]
        old_budget = M.MAX_DRAW_INDICES
        try:
            M.MAX_DRAW_INDICES = 3
            with self.assertRaises(ValueError):
                M.bootstrap_components(components, replicates=2)
        finally:
            M.MAX_DRAW_INDICES = old_budget

    def test_protocol_rejects_wrong_seed_weighting_status_context(self):
        M.check_protocol(protocol())
        for mutation in ["seed", "weighting", "status", "contexts"]:
            p = protocol()
            if mutation == "seed":
                p["analysis_contract"]["bootstrap"]["seed"] += 1
            elif mutation == "weighting":
                p["analysis_contract"]["primary_weighting"] = "equal_variant"
            elif mutation == "status":
                p["status"] = "draft"
            else:
                p["contexts"]["primary"] = ["K562"]
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                M.check_protocol(p)

    def make_handoff(self, root):
        proto = protocol()
        (root / "protocol.json").write_bytes(M.encoded(proto))
        binding = {"protocol_id": "variant-transfer-v1", "protocol_sha256": M.sha256(root / "protocol.json"),
                   "source_cohort_gate": {"passed": True, "description": "synthetic-only fixture"}}
        (root / "binding.json").write_bytes(M.encoded(binding))
        pins = {"schema_version": 1, "protocol_sha256": M.sha256(root / "protocol.json"),
                "binding_sha256": M.sha256(root / "binding.json")}
        variants = fixture()
        for model in ["I", "S"]:
            directory = root / model
            directory.mkdir()
            training = {"artifacts_sha256": "a" * 64, "protocol_sha256": "b" * 64}
            prediction = {"schema_version": 1, "model": model,
                "endpoint_direction": "alternate minus reference",
                "commitments": training | {"pairs_sha256": "c" * 64},
                "pairs": [{"id": row["id"] + "::" + context, "context": context,
                           "predicted_alternate_minus_reference": row["predictions"][model][context]}
                          for row in variants for context in M.CONTEXTS]}
            (directory / "predictions.json").write_bytes(M.encoded(prediction))
            manifest = {"schema_version": 1, "training": training, "pairs_sha256": "c" * 64,
                        "files": {"predictions.json": M.sha256(directory / "predictions.json")}}
            (directory / "manifest.json").write_bytes(M.encoded(manifest))
            pins[f"prediction_{model}_sha256"] = M.sha256(directory / "predictions.json")
            pins[f"prediction_{model}_manifest_sha256"] = M.sha256(directory / "manifest.json")
        normalized = payload(variants) | {key: pins[key] for key in ["protocol_sha256", "binding_sha256",
            "prediction_I_sha256", "prediction_S_sha256"]}
        (root / "input.json").write_bytes(M.encoded(normalized))
        pins["input_sha256"] = M.sha256(root / "input.json")
        (root / "pins.json").write_bytes(M.encoded(pins))
        return types.SimpleNamespace(input=root / "input.json", protocol=root / "protocol.json",
            binding=root / "binding.json", pins=root / "pins.json", prediction_I=root / "I", prediction_S=root / "S")

    def test_handoff_verifies_every_pin_and_exact_prediction_bits(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            args = self.make_handoff(root)
            variants, pins = M.verify_handoff(args)
            self.assertEqual(len(variants), 3)
            changed = M.load_json(args.input)
            # Numerically equal signed zero must still preserve exact producer bits.
            changed["variants"][0]["predictions"]["S"]["K562"] = -0.0
            args.input.write_bytes(M.encoded(changed))
            pins["input_sha256"] = M.sha256(args.input)
            args.pins.write_bytes(M.encoded(pins))
            with self.assertRaisesRegex(ValueError, "normalized/predictor"):
                M.verify_handoff(args)

    def test_gate_failure_and_prediction_tamper_are_not_performance_failures(self):
        with tempfile.TemporaryDirectory() as name:
            args = self.make_handoff(Path(name))
            (args.prediction_I / "predictions.json").write_bytes(b"tampered")
            with self.assertRaisesRegex(ValueError, "prediction file pin mismatch"):
                M.verify_handoff(args)
        with tempfile.TemporaryDirectory() as name:
            args = self.make_handoff(Path(name))
            binding = M.load_json(args.binding)
            binding["source_cohort_gate"]["passed"] = False
            args.binding.write_bytes(M.encoded(binding))
            pins = M.load_json(args.pins)
            pins["binding_sha256"] = M.sha256(args.binding)
            args.pins.write_bytes(M.encoded(pins))
            with self.assertRaisesRegex(ValueError, "handoff gate not passed"):
                M.verify_handoff(args)

    def test_shared_locus_cannot_be_split_across_components(self):
        p = payload()
        p["variants"][2]["locus"] = p["variants"][0]["locus"]
        with self.assertRaisesRegex(ValueError, "split across components"):
            M.validate_variants(p)

    def test_cli_retains_all_inputs_observations_and_2000_shared_draws(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            args = self.make_handoff(root)
            output = root / "evaluation"
            M.main(["--input", str(args.input), "--protocol", str(args.protocol),
                    "--binding", str(args.binding), "--pins", str(args.pins),
                    "--prediction-I", str(args.prediction_I), "--prediction-S", str(args.prediction_S),
                    "--out", str(output)])
            manifest = M.load_json(output / "manifest.json")
            for filename, digest in manifest["files"].items():
                self.assertEqual(M.sha256(output / filename), digest)
            self.assertEqual((output / "normalized-input.json").read_bytes(), args.input.read_bytes())
            with gzip.open(output / "evaluation.json.gz", "rt") as handle:
                result = json.load(handle)
            self.assertEqual(len(result["per_variant"]), 3)
            self.assertEqual(len(result["resampling"]["resampling_indices"]), 2000)
            self.assertEqual(len(result["resampling"]["sampled_statistics"]), 2000)
            self.assertEqual(result["primary"]["conditional_95_percent_range"], [-3.0, 1.5])


if __name__ == "__main__":
    unittest.main()
