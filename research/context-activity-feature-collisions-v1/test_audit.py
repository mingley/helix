from fractions import Fraction
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import audit


def snapshot(identifier, sequence):
    counts = audit.features(sequence)
    return {"id": identifier, "length": 200, "canonical_trimer_counts": list(counts[:32]),
            "gc_count": counts[32], "vector": [c / 198 for c in counts[:32]] + [counts[32] / 200]}


class FeatureAndFloorControls(unittest.TestCase):
    def test_exact_counts_and_reverse_complement(self):
        seq = ("ACGTTACGAGTC" * 17)[:200]
        self.assertEqual(audit.features(seq), audit.features(audit.reverse_complement(seq)))
        self.assertEqual(sum(audit.features(seq)[:32]), 198)
        self.assertEqual(audit.features(seq)[32], seq.count("C") + seq.count("G"))
        with self.assertRaisesRegex(ValueError, "literal 200-base"):
            audit.features("A" * 199 + "N")

    def test_exact_binary64_labels_not_decimal_spellings(self):
        first, second = Fraction.from_float(1.838), Fraction.from_float(1.29)
        difference = first - second
        self.assertNotEqual(difference, Fraction(137, 250))  # decimal 0.548 is not exact here
        sse = audit.squared_error_sum([first, second])
        self.assertEqual(sse / 2, difference ** 2 / 4)
        self.assertAlmostEqual(float(sse / 2), 0.075076)
        encoded = audit.rational(sse / 2)
        self.assertTrue(all(isinstance(value, str) for value in encoded.values()))
        self.assertEqual(Fraction(int(encoded["numerator"]), int(encoded["denominator"])), sse / 2)

    def test_class_weighting_and_singletons(self):
        values = [Fraction(0), Fraction(0), Fraction(1), Fraction(1)]
        sse = audit.squared_error_sum(values)
        self.assertEqual(sse / 4, Fraction(1, 4))
        self.assertEqual(audit.squared_error_sum([Fraction(123, 10)]), 0)
        with self.assertRaisesRegex(ValueError, "empty"):
            audit.squared_error_sum([])

    def test_fingerprint_collision_is_detected_before_assignment(self):
        panel = {"a": {"sequence": "A" * 200}, "b": {"sequence": "A" * 199 + "C"}}
        data = {"feature_names": audit.FEATURE_NAMES,
                "records": [snapshot(i, r["sequence"]) for i, r in panel.items()]}
        with patch.object(audit, "read_json", return_value=data), \
             patch.object(audit, "feature_fingerprint", return_value="same-digest"):
            with self.assertRaisesRegex(ValueError, "different exact feature tuples"):
                audit.compute_design(panel)

    def test_group_crossing_test_folds_is_an_error(self):
        data = {"records": [{"id": "a", "group": "g"}, {"id": "b", "group": "g"}],
                "protocol_sha256": audit.EXPECTED["grouping_protocol"],
                "panel_sha256": audit.EXPECTED["panel"],
                "folds": [{"fold": 0, "test_ids": ["a"]}, {"fold": 1, "test_ids": ["b"]}]}
        with patch.object(audit, "read_json", return_value=data):
            with self.assertRaisesRegex(ValueError, "crosses test folds"):
                audit.groups_and_folds({"a", "b"})

    def test_panel_metadata_reader_does_not_interpret_old_observations(self):
        sequence = "A" * 200
        source = {"s7_row_1based": 2, "s9_row_1based": 3, "reference_match_strand": "+",
                  "source_declared_strand": "+", "source_design_metadata_record": "example",
                  "core_sha256": hashlib.sha256(sequence.encode()).hexdigest()}
        record = {"id": "example", "sequence": sequence, "chromosome": "chr1", "start": 10,
                  "end": 210, "parent_id": "example", "observations": {"not_used": "sentinel"},
                  "source_metadata": source}
        header = ('{"schema_version":1,"endpoint":"published_element_specificity_score",'
                  '"contexts":["HepG2","K562","WTC11"],"records":[\n')
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "panel.json"
            path.write_text(header + json.dumps(record, separators=(",", ":")) + "\n]}\n")
            with patch.object(audit, "PANEL", path):
                data = audit.read_panel()
        self.assertEqual(set(data), {"example"})
        self.assertNotIn("observations", data["example"])
        self.assertEqual(data["example"]["sequence"], sequence)


if __name__ == "__main__":
    unittest.main()
