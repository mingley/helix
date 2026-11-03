from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

from verify import (canonical, check_equal, fingerprint, floor_for_equal_predictions, fraction_object,
                    reverse_complement, source_activity, vector)


def synthetic_activity_workbook(path, formula=False):
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    header = ("name", "HepG2 [log2(rna/dna)]", "K562 [log2(rna/dna)]",
              "WTC11 [log2(rna/dna)]", "example")
    strings = '<sst xmlns="' + ns + '">' + ''.join('<si><t>' + s + '</t></si>' for s in header) + '</sst>'
    workbook = ('<workbook xmlns="' + ns + '" xmlns:r="http://schemas.openxmlformats.org/'
                'officeDocument/2006/relationships"><sheets><sheet name="all_cell_types_summary" '
                'sheetId="1" r:id="rId1"/></sheets></workbook>')
    rels = ('<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
    first = ''.join(f'<c r="{col}1" t="s"><v>{i}</v></c>' for i, col in enumerate("ABCD"))
    b_value = '<f>1+1</f><v>2</v>' if formula else '<v>1.8380000000000001</v>'
    rows = (f'<row r="1">{first}</row><row r="2"><c r="A2" t="s"><v>4</v></c>'
            f'<c r="B2">{b_value}</c><c r="C2"><v>1.8759999999999999</v></c>'
            '<c r="D2"><v>3.1419999999999999</v></c></row>')
    sheet = '<worksheet xmlns="' + ns + '"><sheetData>' + rows + '</sheetData></worksheet>'
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", rels)
        archive.writestr("xl/sharedStrings.xml", strings)
        archive.writestr("xl/worksheets/sheet1.xml", sheet)


class IndependentCollisionTests(unittest.TestCase):
    def test_raw_features_reverse_complements_and_non_rc_collisions(self):
        sequence = "AACGTACGTA"
        self.assertEqual(vector(sequence), vector(reverse_complement(sequence)))
        self.assertEqual(vector("ACAC"), vector("CACA"))
        self.assertNotEqual("ACAC", reverse_complement("CACA"))
        self.assertNotEqual(vector("AAAA"), vector("AAAC"))
        self.assertEqual(canonical("ACG"), canonical("CGT"))
        with self.assertRaises(ValueError):
            vector("ACGN")

    def test_exact_binary64_floor_and_invalid_inputs(self):
        self.assertEqual(floor_for_equal_predictions([1.0, 3.0]), Fraction(2))
        first, second = 1.838, 1.29
        delta = Fraction.from_float(first) - Fraction.from_float(second)
        self.assertEqual(floor_for_equal_predictions([first, second]), delta * delta / 2)
        self.assertNotEqual(delta, Fraction(548, 1000))
        encoded = fraction_object(delta * delta / 2)
        self.assertIsInstance(encoded["numerator"], str)
        self.assertIsInstance(encoded["denominator"], str)
        self.assertEqual(Fraction(int(encoded["numerator"]), int(encoded["denominator"])), delta * delta / 2)
        with self.assertRaises(ValueError):
            floor_for_equal_predictions([1.0])

    def test_original_xml_lexemes_decode_to_expected_binary64(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.xlsx"
            synthetic_activity_workbook(path)
            rows = source_activity(path, ["example"])
            self.assertEqual(rows["example"]["row"], 2)
            raw = rows["example"]["lexemes"]
            self.assertEqual(raw["HepG2"], "1.8380000000000001")
            self.assertEqual(float(raw["HepG2"]).hex(), float("1.838").hex())
            self.assertEqual(float(raw["K562"]).hex(), float("1.876").hex())
            synthetic_activity_workbook(path, formula=True)
            with self.assertRaisesRegex(ValueError, "formula"):
                source_activity(path, ["example"])

    def test_altered_assignment_or_summary_is_rejected(self):
        assignment = {"id": "x", "feature_fingerprint": fingerprint(vector("ACAC")),
                      "class_size": 2, "group": "g", "test_fold": 1}
        check_equal(dict(assignment), assignment, "assignment")
        for key, replacement in (("class_size", 1), ("test_fold", 2),
                                 ("feature_fingerprint", fingerprint(vector("AAAA")))):
            altered = dict(assignment, **{key: replacement})
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "differs"):
                check_equal(altered, assignment, "assignment")
        with self.assertRaises(ValueError):
            check_equal({"floor": {"numerator": 1, "denominator": 2}},
                        {"floor": {"numerator": 1, "denominator": 3}}, "summary")

    def test_existing_receipt_cannot_be_overwritten_under_optimization(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "protected.json"
            path.write_text("do not overwrite\n")
            script = Path(__file__).with_name("verify.py")
            result = subprocess.run([sys.executable, "-O", str(script), "--receipt", str(path)],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("receipt path already exists", result.stderr)
            self.assertEqual(path.read_text(), "do not overwrite\n")


if __name__ == "__main__":
    unittest.main()
