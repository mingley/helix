import gzip
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import verify


def group(center="var1", window="left"):
    return {"variant1": "chr1:100:A:C", "variant2": "chr2:220:G:T",
            "center": center, "window": window}


class NameConstructionControls(unittest.TestCase):
    def test_centered_variant_uses_correct_haplotype_index(self):
        self.assertEqual(verify.syntactic_names(group()), (
            "1:100:A:C:R:wL:Alt_2-220-R",
            "1:100:A:C:A:wL:Alt_2-220-R",
            "1:100:A:C:R:wL:Alt_2-220-A",
            "1:100:A:C:A:wL:Alt_2-220-A",
        ))
        self.assertEqual(verify.syntactic_names(group("var2", "right")), (
            "2:220:G:T:R:wR:Alt_1-100-R",
            "2:220:G:T:R:wR:Alt_1-100-A",
            "2:220:G:T:A:wR:Alt_1-100-R",
            "2:220:G:T:A:wR:Alt_1-100-A",
        ))

    def test_invalid_source_names_are_rejected(self):
        bad = group()
        bad["variant1"] = "chr1:100:A"
        with self.assertRaisesRegex(ValueError, "variant identifier"):
            verify.syntactic_names(bad)
        bad = group("var3")
        with self.assertRaisesRegex(ValueError, "center or window"):
            verify.syntactic_names(bad)


class ParsingAndOverlapControls(unittest.TestCase):
    def test_single_name_flag_zero_requires_both_conditions(self):
        # Flag 0 with multiple names and flag 1 with one name demonstrate that
        # neither predicate alone confers the strict single-name condition.
        names = {"one", "two", "three", "four", "five"}
        rows = [
            ["bc", "one", "x", "x", "0", "x", "x", "x", "x", "x"],
            ["bc", "two", "x", "x", "1", "x", "x", "x", "x", "x"],
            ["bc", "three,four", "x", "x", "0", "x", "x", "x", "x", "x"],
            ["bc", "five", "x", "x", "2", "x", "x", "x", "x", "x"],
        ]
        raw = "".join("\t".join(row) + "\r\n" for row in rows).encode()
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "lookup.gz"
            source.write_bytes(gzip.compress(raw))
            seen, strict, totals, md5, size = verify.scan_lookup(source, names)
        self.assertEqual(seen, names)
        self.assertEqual(strict, {"one"})
        self.assertEqual(totals["rows"], 4)
        self.assertEqual(totals["flag_0_rows"], 2)
        self.assertEqual(totals["multiple_name_rows"], 1)
        self.assertEqual((md5, size), (hashlib.md5(raw).hexdigest(), len(raw)))

    def test_malformed_rows_and_flags_fail(self):
        with self.assertRaisesRegex(ValueError, "ten-column"):
            verify.fields(["x"] * 9)
        with self.assertRaisesRegex(ValueError, "mapping flag"):
            verify.fields(["x", "candidate", "x", "x", "4", "x", "x", "x", "x", "x"])
        with self.assertRaisesRegex(ValueError, "empty candidate"):
            verify.fields(["x", "a,", "x", "x", "0", "x", "x", "x", "x", "x"])

    def test_full_and_partial_windows_and_collision_counts(self):
        first, second = group(), group("var2", "middle")
        classes = {verify.key(first): "first_blind_only", verify.key(second): "no_collision"}
        seen = set(verify.syntactic_names(first)) | set(verify.syntactic_names(second)[:2])
        result = verify.counts([first, second], classes, seen)
        self.assertEqual(result["windows_by_number_of_four_candidate_strings_present"]["4"], 1)
        self.assertEqual(result["windows_by_number_of_four_candidate_strings_present"]["2"], 1)
        self.assertEqual(result["pair_ids_with_at_least_one_four_string_window"], 1)
        self.assertEqual(result["collision_windows_by_number_of_four_candidate_strings_present"]["4"], 1)
        self.assertEqual(result["collision_pair_ids_with_four_string_window"], 1)
        self.assertEqual(result["collision_classifications_for_four_string_windows"], {"first_blind_only": 1})

    def test_numeric_inventory_covers_nested_fields(self):
        self.assertEqual(verify.numeric_paths({"a": 1, "b": {"c": 2}, "d": "2", "e": True}),
                         {("a",): 1, ("b", "c"): 2})


class ReceiptIntegrityControls(unittest.TestCase):
    def test_existing_receipt_is_never_overwritten_under_optimization(self):
        script = Path(__file__).with_name("verify.py")
        with tempfile.TemporaryDirectory() as tmp:
            receipt = Path(tmp) / "receipt.json"
            receipt.write_text("existing\n")
            result = subprocess.run([sys.executable, "-O", str(script), "--lookup", "missing",
                                     "--receipt", str(receipt)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("receipt already exists", result.stderr)
            self.assertEqual(receipt.read_text(), "existing\n")

    def test_exclusive_create_catches_race_after_precheck(self):
        with tempfile.TemporaryDirectory() as tmp:
            receipt = Path(tmp) / "receipt.json"
            def concurrent_creation(args):
                args.receipt.write_text("concurrent receipt\n")
                return {"all_checks_passed": True}
            with patch.object(sys, "argv", ["verify.py", "--lookup", "missing", "--receipt", str(receipt)]), \
                    patch.object(verify, "run", side_effect=concurrent_creation):
                with self.assertRaises(FileExistsError):
                    verify.main()
            self.assertEqual(receipt.read_text(), "concurrent receipt\n")


if __name__ == "__main__":
    unittest.main()
