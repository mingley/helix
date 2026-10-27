import copy
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location("variant_assemble", Path(__file__).resolve().parents[2] /
                                             "scripts/assemble_variant_transfer.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def fixture():
    rows = []
    for identity, complete in (("a", True), ("b", True), ("missing", False)):
        observed = {"K562": -0.0, "HepG2": 0.125}
        rows.append({"id": identity, "component": "one-original-component", "locus": identity,
                     "gc_preserving": True, "complete_case": complete, "observed": observed,
                     "source_float64_le_hex": {cell: M.bits(value) for cell, value in observed.items()}})
    predictions = {model: {identity + "::" + cell: {"context": cell,
                   "predicted_alternate_minus_reference": 0.003}
                   for identity in ("a", "b") for cell in M.CELLS} for model in ("I", "S")}
    return {"variants": rows}, predictions


class AssembleTests(unittest.TestCase):
    def test_exact_floats_component_labels_and_complete_cases_preserved(self):
        observations, predictions = fixture()
        result = M.join(observations, predictions)
        self.assertEqual([r["id"] for r in result], ["a", "b"])
        self.assertEqual({r["component"] for r in result}, {"one-original-component"})
        self.assertEqual(M.bits(result[0]["observed"]["K562"]), M.bits(-0.0))
        self.assertEqual(result[0]["predictions"]["I"]["HepG2"], 0.003)

    def test_changed_source_float_bits_and_duplicate_observation_rejected(self):
        observations, predictions = fixture()
        observations["variants"][0]["observed"]["K562"] = 0.0
        with self.assertRaisesRegex(ValueError, "source float bits"):
            M.join(observations, predictions)
        observations, predictions = fixture()
        observations["variants"].append(copy.deepcopy(observations["variants"][0]))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            M.join(observations, predictions)

    def test_extra_wrong_context_and_nonfinite_predictions_rejected(self):
        for mutation in ("extra", "context", "nan"):
            observations, predictions = fixture()
            if mutation == "extra":
                predictions["I"]["unexplained::K562"] = {"context": "K562",
                    "predicted_alternate_minus_reference": 0.0}
            elif mutation == "context":
                predictions["S"]["a::K562"]["context"] = "HepG2"
            else:
                predictions["I"]["a::K562"]["predicted_alternate_minus_reference"] = float("nan")
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                M.join(observations, predictions)


if __name__ == "__main__":
    unittest.main()
