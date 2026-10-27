"""Synthetic controls for transfer arithmetic and strict outcome-free inputs."""
import copy
import importlib.util
import itertools
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("variant_predict", Path(__file__).resolve().parents[2] / "scripts" / "predict_variant_effects.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def synthetic_fits():
    ridge = {"training_rows": 10, "feature_count": 33, "alpha": 10.0,
             "intercept": 7.0, "training_feature_means": [0.03] * 33,
             "training_feature_scales": [0.1 + j / 100 for j in range(33)],
             "standardized_coefficients": [(j - 16) / 10 for j in range(33)]}
    return [{"shared_dna_fit": copy.deepcopy(ridge),
             "training_context_means": {c: 10.0 * j for j, c in enumerate(M.CONTEXTS)},
             "independent_context_fits": {c: copy.deepcopy(ridge) for c in M.CONTEXTS}}
            for _ in range(5)]


def pair(reference=None, position=100, alt="C", context="K562"):
    reference = reference or "A" * 200
    return {"id": "synthetic", "context": context, "position_0based": position,
            "reference": reference,
            "alternate": reference[:position] + alt + reference[position + 1:]}


class VariantTransferTests(unittest.TestCase):
    def score(self, p, f=None, model="I"):
        M.validate_pairs({"schema_version": 1, "pairs": [p]})
        return M.score_pair(p, f or synthetic_fits(), model)

    def test_feature_exact_counts_and_reverse_complement(self):
        s = "AAACCCGT" * 25
        x = M.features(s)
        self.assertEqual(len(x), 33)
        self.assertAlmostEqual(sum(x[:-1]), 1)
        self.assertEqual(x, M.features(M.reverse_complement(s)))
        self.assertEqual(x[-1], 0.5)

    def test_delta_matches_two_direct_predictions(self):
        scored = self.score(pair())
        self.assertLess(scored["max_abs_direct_subtraction_check_error"], 1e-12)
        self.assertNotEqual(scored["predicted_alternate_minus_reference"], 0)

    def test_saved_exact_numerators_bridge_all_prediction_features(self):
        scored = self.score(pair())
        ref = scored["reference_feature_numerators"]
        alt = scored["alternate_feature_numerators"]
        delta = scored["feature_numerator_delta"]
        self.assertEqual(sum(ref[:-1]), 198)
        self.assertEqual(ref[M.TRIMERS.index("AAA")], 198)
        self.assertEqual(ref[-1], 0)
        self.assertEqual(sum(alt[:-1]), 198)
        self.assertEqual(delta[M.TRIMERS.index("AAA")], -3)
        self.assertEqual(delta[-1], 1)
        self.assertEqual(delta, [b - a for a, b in zip(ref, alt)])
        self.assertTrue(all(type(value) is int for value in ref + alt + delta))

    def test_intercept_and_centering_cancel(self):
        f = synthetic_fits()
        old = self.score(pair(), f)["predicted_alternate_minus_reference"]
        for fit in f:
            for ridge in [fit["shared_dna_fit"], *fit["independent_context_fits"].values()]:
                ridge["intercept"] += 1000
                ridge["training_feature_means"] = [1.0] * 33
        self.assertEqual(old, self.score(pair(), f)["predicted_alternate_minus_reference"])

    def test_local_five_mer_edit_independent_of_remote_context(self):
        a = "A" * 98 + "TACGT" + "A" * 97
        b = "G" * 98 + "TACGT" + "T" * 97
        self.assertEqual(self.score(pair(a, 100, "A"))["predicted_alternate_minus_reference"],
                         self.score(pair(b, 100, "A"))["predicted_alternate_minus_reference"])

    def test_allele_reversal_negates_effect(self):
        p = pair()
        q = p | {"reference": p["alternate"], "alternate": p["reference"]}
        self.assertEqual(self.score(p)["predicted_alternate_minus_reference"],
                         -self.score(q)["predicted_alternate_minus_reference"])

    def test_reverse_complement_pair_same_prediction(self):
        p = pair("ACGT" * 50, 100, "G")
        q = p | {"position_0based": 99, "reference": M.reverse_complement(p["reference"]),
                 "alternate": M.reverse_complement(p["alternate"])}
        self.assertEqual(self.score(p)["predicted_alternate_minus_reference"],
                         self.score(q)["predicted_alternate_minus_reference"])

    def test_boundary_edits_use_actual_complete_windows(self):
        for position in [0, 1, 198, 199]:
            scored = self.score(pair(position=position))
            self.assertLess(scored["max_abs_direct_subtraction_check_error"], 1e-12)
        self.assertNotEqual(self.score(pair(position=0))["predicted_alternate_minus_reference"],
                            self.score(pair(position=100))["predicted_alternate_minus_reference"])

    def test_feature_invisible_local_edit_is_zero(self):
        p = pair("G" * 98 + "AAATT" + "C" * 97, 100, "T")
        self.assertEqual(M.features(p["reference"]), M.features(p["alternate"]))
        self.assertEqual(self.score(p)["predicted_alternate_minus_reference"], 0)

    def test_all_64_feature_invisible_oriented_motifs(self):
        zero = 0
        for word in map("".join, itertools.product("ACGT", repeat=5)):
            for alt in "ACGT":
                if alt == word[2]:
                    continue
                p = pair("A" * 98 + word + "A" * 97, 100, alt)
                zero += M.features(p["reference"]) == M.features(p["alternate"])
        self.assertEqual(zero, 64)

    def test_equal_weight_ensemble_and_gc_ablation_label(self):
        f = synthetic_fits()
        for j, fit in enumerate(f):
            fit["independent_context_fits"]["K562"]["standardized_coefficients"] = [j] * 33
        p = pair()
        scored = self.score(p, f)
        self.assertEqual(scored["predicted_alternate_minus_reference"], sum(scored["fold_predictions"]) / 5)
        expected = (1 / 200) * sum(j / f[j]["independent_context_fits"]["K562"]["training_feature_scales"][-1] for j in range(5)) / 5
        self.assertAlmostEqual(scored["learned_gc_term_only_ablation_not_refitted_gc_model"], expected)

    def test_shared_additive_deltas_identical(self):
        self.assertEqual(self.score(pair(), model="S")["predicted_alternate_minus_reference"],
                         self.score(pair(), model="A")["predicted_alternate_minus_reference"])

    def test_strict_input_rejects_outcomes_duplicate_ids_and_indels(self):
        invalid = [pair() | {"observed": 0}, pair() | {"alternate": "A" * 199},
                   pair() | {"position_0based": True}, pair() | {"position_0based": 101},
                   pair() | {"alternate": "A" * 200}, pair() | {"reference": "a" * 200}]
        for p in invalid:
            with self.subTest(p=p), self.assertRaises(ValueError):
                M.validate_pairs({"schema_version": 1, "pairs": [p]})
        with self.assertRaises(ValueError):
            M.validate_pairs({"schema_version": 1, "pairs": [pair(), pair()]})

    def test_unknown_context_model_and_invalid_scales_fail(self):
        with self.assertRaises(ValueError):
            self.score(pair(context="unknown"))
        with self.assertRaises(ValueError):
            self.score(pair(), model="bad")
        bad = synthetic_fits()[0]["shared_dna_fit"]
        bad["training_feature_scales"][0] = 0
        with self.assertRaises(ValueError):
            M.validate_ridge(bad)

    def test_manifest_pin_and_path_escape_rejected_before_fit_parsing(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            (root / "artifacts.json").write_text('{"schema_version":5,"files":{}}')
            with self.assertRaisesRegex(ValueError, "manifest SHA256 mismatch"):
                M.load_frozen_fits(root, "0" * 64)
            with self.assertRaises(ValueError):
                M.safe_member(root, "../artifacts.json")
            (root / "link").symlink_to(root / "artifacts.json")
            with self.assertRaises(ValueError):
                M.safe_member(root, "link")

    def test_manifest_member_tamper_fails_before_invalid_model_is_parsed(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            files = {}
            for member in ["study.json.gz", "protocol.json", "source-manifest.json", "inputs.json"]:
                (root / member).write_bytes(b"not JSON or gzip")
                files[member] = M.sha256(root / member)
            (root / "artifacts.json").write_bytes(M.canonical_json({"schema_version": 5, "files": files}))
            pin = M.sha256(root / "artifacts.json")
            (root / "study.json.gz").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "artifact SHA256 mismatch: study.json.gz"):
                M.load_frozen_fits(root, pin)

    def test_duplicate_json_keys_and_nonfinite_literals_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            path = Path(name) / "input.json"
            for invalid in ['{"pairs":[],"pairs":[]}', '{"observed":NaN}']:
                path.write_text(invalid)
                with self.assertRaises(ValueError):
                    M.load_json(path)


if __name__ == "__main__":
    unittest.main()
