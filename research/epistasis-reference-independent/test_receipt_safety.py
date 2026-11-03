"""Receipt writes must not overwrite inputs or an existing file."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from verify import require_new_receipt, write_receipt


class ReceiptSafety(unittest.TestCase):
    def test_existing_input_and_aliases_are_rejected_under_optimization(self):
        script_dir = Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            source = directory / "source.xlsx"
            source.write_bytes(b"original source contents\n")
            symbolic = directory / "symbolic.xlsx"
            symbolic.symlink_to(source)
            hardlink = directory / "hardlink.xlsx"
            os.link(source, hardlink)
            reference = directory / "reference.2bit"
            reference.write_bytes(b"original reference contents\n")
            for script in ("verify.py", "verify_hg19.py"):
                for destination in (source, symbolic, hardlink, reference):
                    with self.subTest(script=script, destination=destination.name):
                        run = subprocess.run(
                            [sys.executable, "-O", str(script_dir / script),
                             "--source", str(source), "--reference", str(reference),
                             "--receipt", str(destination)],
                            capture_output=True, text=True,
                        )
                        self.assertNotEqual(run.returncode, 0, run.stdout)
                        self.assertIn("receipt path already exists", run.stderr)
                        self.assertEqual(source.read_bytes(), b"original source contents\n")
                        self.assertEqual(reference.read_bytes(), b"original reference contents\n")

    def test_exclusive_write_protects_against_a_new_file_after_preflight(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "receipt.json"
            require_new_receipt(output)
            output.write_text("created after preflight\n")
            with self.assertRaises(FileExistsError):
                write_receipt(output, {"all_checks_passed": True})
            self.assertEqual(output.read_text(), "created after preflight\n")

    def test_new_receipt_is_created_once(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "receipt.json"
            receipt = {"all_checks_passed": True}
            require_new_receipt(output)
            write_receipt(output, receipt)
            self.assertEqual(json.loads(output.read_text()), receipt)
            with self.assertRaises(FileExistsError):
                require_new_receipt(output)
            self.assertEqual(json.loads(output.read_text()), receipt)


if __name__ == "__main__":
    unittest.main()
