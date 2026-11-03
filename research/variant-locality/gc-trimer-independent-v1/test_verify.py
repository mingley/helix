import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction

from verify import canonical_counts, compare_summary, inspect_sequence, reparameterized_scores, verify_coefficients


class IndependentIdentityTests(unittest.TestCase):
    def test_raw_windows_and_boundary_multiplicity(self):
        self.assertEqual(inspect_sequence("AAAA"), (0, 0, 0))
        self.assertEqual(inspect_sequence("CCCC"), (4, 6, 6))
        self.assertEqual(inspect_sequence("ACGT"), (2, 4, 2))
        self.assertEqual(inspect_sequence("ACGTA"), (2, 5, 1))
        for invalid in ("ACG", "ACGN", "acgt", None):
            with self.subTest(sequence=invalid), self.assertRaises(ValueError):
                inspect_sequence(invalid)

    def test_exact_reparameterization_on_exhaustive_small_gc_patterns(self):
        for length in range(4, 9):
            verify_coefficients(length)
            for letters in itertools.product("AC", repeat=length):
                sequence = "".join(letters)
                gc, windows, edge = inspect_sequence(sequence)
                self.assertEqual(3 * gc, windows + edge)
                old, new = reparameterized_scores(sequence, Fraction(-3, 8))
                self.assertEqual(old, new)

    def test_equal_trimermultisets_can_have_distinct_boundary_and_gc(self):
        first = "".join("AAC"[i % 3] for i in range(200))
        second = "".join("AAC"[(i + 1) % 3] for i in range(200))
        self.assertEqual(canonical_counts(first), canonical_counts(second))
        self.assertEqual(inspect_sequence(first), (66, 198, 0))
        self.assertEqual(inspect_sequence(second), (67, 198, 3))
        old_first, new_first = reparameterized_scores(first)
        old_second, new_second = reparameterized_scores(second)
        self.assertEqual(old_first, new_first)
        self.assertEqual(old_second, new_second)
        self.assertEqual(old_second - old_first, Fraction(5, 7 * 200))

    def test_summary_tamper_is_rejected(self):
        expected = {"boundary_count_histogram": {"0": 2, "1": 1},
                    "correlation_squared_exact": {"numerator": "3", "denominator": "4"},
                    "same_trimers_different_GC_example": {"second_GC_count": 67}}
        compare_summary(json.loads(json.dumps(expected)), expected)
        for key, value in (("boundary_count_histogram", {"0": 3, "1": 1}),
                           ("correlation_squared_exact", {"numerator": "4", "denominator": "4"}),
                           ("correlation_squared_exact", {"numerator": 3, "denominator": 4}),
                           ("same_trimers_different_GC_example", {"second_GC_count": 66})):
            tampered = dict(expected, **{key: value})
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "primary summary differs"):
                compare_summary(tampered, expected)

    def test_tampered_primary_summary_cannot_create_receipt(self):
        base = Path(__file__).resolve().parent
        summary_path = base.parent / "gc-trimer-identity-v1/summary.json"
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            summary = json.loads(summary_path.read_text())
            summary["correlation_sufficient_statistics"]["sum_G"] += 1
            altered = directory / "altered-summary.json"
            altered.write_text(json.dumps(summary))
            receipt = directory / "receipt.json"
            result = subprocess.run([sys.executable, "-O", str(base / "verify.py"),
                                     "--summary", str(altered), "--receipt", str(receipt)],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn("primary summary differs", result.stderr)
            self.assertFalse(receipt.exists())

    def test_receipt_does_not_overwrite_input_under_optimization(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "panel.json"
            path.write_text("original source\n")
            result = subprocess.run([sys.executable, "-O", str(Path(__file__).with_name("verify.py")),
                                     "--panel", str(path), "--receipt", str(path)],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("receipt path already exists", result.stderr)
            self.assertEqual(path.read_text(), "original source\n")


if __name__ == "__main__":
    unittest.main()
