"""Tiny synthetic XLSX controls; no empirical S4 payload is opened."""
import gzip
import importlib.util
import json
from pathlib import Path
import struct
import tempfile
import types
import unittest
from unittest import mock
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape
import zipfile

SPEC = importlib.util.spec_from_file_location("variant_observations", Path(__file__).resolve().parents[2] / "scripts" / "build_variant_transfer_observations.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def text(column, row, value):
    return f'<c r="{column}{row}" t="inlineStr"><is><t>{escape(value)}</t></is></c>'


def number(column, row, value, kind=None, formula=None):
    typed = f' t="{kind}"' if kind is not None else ""
    f = f'<f>{escape(formula)}</f>' if formula is not None else ""
    v = f'<v>{escape(value)}</v>' if value is not None else ""
    return f'<c r="{column}{row}"{typed}>{f}{v}</c>'


def workbook(path, rows=None, headers=None, name=None, extra_entries=None):
    headers = M.HEADERS if headers is None else headers
    rows = rows if rows is not None else [(3, text("A", 3, "v1") + number("AN", 3, "-0") + number("AO", 3, "1.25e-1"))]
    header = '<row r="2">' + "".join(text(c, 2, value) for c, value in headers.items()) + '</row>'
    data = "".join(f'<row r="{row}">{cells}</row>' for row, cells in rows)
    sheet = f'<worksheet xmlns="{M.MAIN}"><dimension ref="A1:XFD999999"/><sheetData>{header}{data}</sheetData></worksheet>'
    name = M.SHEET if name is None else name
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("xl/workbook.xml", f'<workbook xmlns="{M.MAIN}" xmlns:r="{M.REL}"><sheets><sheet name="{name}" sheetId="1" r:id="rId1"/></sheets></workbook>')
        archive.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        archive.writestr("xl/worksheets/sheet1.xml", sheet)
        for member, value in (extra_entries or {}).items():
            archive.writestr(member, value)


def metadata(identity="v1", component="group1"):
    return {"id": identity, "component": component, "locus": "locus" + component,
            "sequence_ref": "A" * 200, "sequence_alt": "A" * 99 + "C" + "A" * 100,
            "gc_preserving": False, "position_0based": 99}


def binding(source, cohort, count=1, tokens=None, protocol_id="variant-transfer-v1"):
    result = {"schema_version": 1, "protocol_id": protocol_id,
            "protocol_sha256": M.PROTOCOLS[protocol_id], "model_artifacts_sha256": M.MODEL_SHA256,
            "source_cohort_gate": {"passed": True},
            "source": {"sha256": M.sha256(source), "bytes": source.stat().st_size},
            "sheet": {"name": M.SHEET, "header_row_1based": 2, "headers": M.HEADERS},
            "effects": M.EFFECTS,
            "missing_policy": {"allowed_tokens": [] if tokens is None else tokens,
                               "evidence": "synthetic source fixture statement"},
            "cohort": {"sha256": M.sha256(cohort), "format": "json", "variant_count": count}}
    if protocol_id == "variant-transfer-v2":
        result["operational_missingness"] = dict(M.OPERATIONAL_MISSINGNESS)
        result["exposure"] = {"v1_failure_observed": True}
    return result


class ObservationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.source = self.root / "source.xlsx"
        self.cohort = self.root / "cohort.json"
        workbook(self.source)
        self.cohort.write_bytes(M.encoded({"schema_version": 1, "pairs": [metadata()]}))

    def tearDown(self):
        self.tmp.cleanup()

    def parse(self, cohorts=None, tokens=None):
        ids, cells = [], []
        selected, integrity, inventory = M.parse_source(self.source, cohorts or {"v1": metadata()}, tokens or [], ids, cells)
        return selected, integrity, inventory, ids, cells

    def test_numeric_lexemes_and_negative_zero_bits_are_exact(self):
        selected, integrity, inventory, ids, cells = self.parse()
        self.assertEqual(selected["v1"]["HepG2"]["raw_lexeme"], "-0")
        self.assertEqual(selected["v1"]["HepG2"]["float64_le_hex"], struct.pack("<d", -0.0).hex())
        self.assertEqual(selected["v1"]["K562"]["value"], .125)
        self.assertTrue(integrity["all_member_crc_checked"])
        self.assertEqual(inventory["source_unique_IDs"], 1)
        self.assertEqual(len(cells), 2)

    def test_unselected_response_cells_are_not_decoded_and_ID_is_preserved(self):
        workbook(self.source, [(3, text("A", 3, "unselected") + number("AN", 3, "nonsense", "b", "BAD()")),
                               (8, text("A", 8, "v1") + number("AN", 8, "1") + number("AO", 8, "2"))])
        selected, _, _, ids, cells = self.parse()
        self.assertEqual(len(ids), 2)
        self.assertFalse(ids[0]["selected_by_frozen_metadata"])
        self.assertEqual({row["id"] for row in cells}, {"v1"})
        self.assertEqual(selected["v1"]["HepG2"]["value"], 1)

    def test_formulas_rejected_even_with_finite_cached_value_or_empty_formula(self):
        for formula in ["SUM(1)", ""]:
            with self.subTest(formula=formula):
                workbook(self.source, [(3, text("A", 3, "v1") + number("AN", 3, "1", formula=formula) + number("AO", 3, "2"))])
                with self.assertRaisesRegex(M.ObservationCellError, "formula"):
                    self.parse()

    def test_boolean_error_date_and_numeric_string_rejected(self):
        for kind, raw in [("b", "1"), ("e", "#N/A"), ("d", "2026-01-01"), ("str", "0.25")]:
            with self.subTest(kind=kind):
                cell = ET.fromstring(f'<c xmlns="{M.MAIN}" r="AN3" t="{kind}"><v>{escape(raw)}</v></c>')
                with self.assertRaises(M.ObservationCellError):
                    M.effect_cell(cell, [], ["NA"])

    def test_shared_numeric_string_is_not_coerced(self):
        cell = ET.fromstring(f'<c xmlns="{M.MAIN}" r="AN3" t="s"><v>0</v></c>')
        with self.assertRaises(M.ObservationCellError):
            M.effect_cell(cell, ["1.25"], ["NA"])

    def test_only_explicitly_bound_exact_NA_is_missing(self):
        cell = ET.fromstring(f'<c xmlns="{M.MAIN}" r="AN3" t="inlineStr"><is><t>NA</t></is></c>')
        with self.assertRaises(M.ObservationCellError):
            M.effect_cell(cell, [], [])
        result = M.effect_cell(cell, [], ["NA"])
        self.assertIsNone(result["value"])
        self.assertEqual(result["text"], "NA")
        for value in ["na", " NA", "N/A", "", "1"]:
            cell = ET.fromstring(f'<c xmlns="{M.MAIN}" r="AN3" t="inlineStr"><is><t>{value}</t></is></c>')
            with self.assertRaises(M.ObservationCellError):
                M.effect_cell(cell, [], ["NA"])

    def test_absent_blank_nonfinite_or_malformed_number_rejected(self):
        with self.assertRaises(M.ObservationCellError):
            M.effect_cell(None, [], ["NA"])
        for raw in [None, "NaN", "inf", "1e999", " 1", "1_000", "1,000"]:
            cell = ET.fromstring(f'<c xmlns="{M.MAIN}" r="AN3">' + (f'<v>{raw}</v>' if raw is not None else '') + '</c>')
            with self.assertRaises(M.ObservationCellError):
                M.effect_cell(cell, [], ["NA"])

    def test_duplicate_numeric_payload_is_not_silently_truncated(self):
        cell = ET.fromstring(f'<c xmlns="{M.MAIN}" r="AN3"><v>1</v><v>2</v></c>')
        with self.assertRaisesRegex(M.ObservationCellError, "duplicate effect"):
            M.effect_cell(cell, [], [])

    def test_duplicate_source_ID_even_unselected_is_gate_failure(self):
        workbook(self.source, [(3, text("A", 3, "x")), (4, text("A", 4, "x"))])
        with self.assertRaisesRegex(ValueError, "duplicate source"):
            self.parse()

    def test_cohort_ID_absence_is_not_silent_completecase_exclusion(self):
        workbook(self.source, [(3, text("A", 3, "other"))])
        with self.assertRaisesRegex(ValueError, "cohort ID absent"):
            self.parse()

    def test_header_identity_sheet_and_cell_coordinates_are_exact(self):
        altered = dict(M.HEADERS, AN="log2Skew in HepG2")
        workbook(self.source, headers=altered)
        with self.assertRaisesRegex(ValueError, "header differs"):
            self.parse()
        workbook(self.source, name="wrong sheet")
        with self.assertRaisesRegex(ValueError, "sheet identity"):
            self.parse()
        workbook(self.source, [(3, text("A", 4, "v1"))])
        with self.assertRaisesRegex(ValueError, "coordinate/row"):
            self.parse()

    def test_duplicate_coordinate_and_duplicate_physical_row_rejected(self):
        for rows in [[(3, text("A", 3, "v1") + text("A", 3, "x"))],
                     [(3, text("A", 3, "x")), (3, text("A", 3, "y"))]]:
            workbook(self.source, rows)
            with self.assertRaisesRegex(ValueError, "duplicate"):
                self.parse()

    def test_sealed_gate_and_binding_hash_checked_before_source_access(self):
        value = binding(self.source, self.cohort)
        value["source_cohort_gate"]["passed"] = False
        file = self.root / "binding.json"
        file.write_bytes(M.encoded(value))
        with mock.patch.object(M.zipfile, "ZipFile", side_effect=AssertionError("source ZIP opened")):
            with self.assertRaisesRegex(ValueError, "gate failed"):
                M.validate_binding(file, M.sha256(file), Path("never-open-protocol"), self.source, self.cohort)
            with self.assertRaisesRegex(ValueError, "binding SHA256"):
                M.validate_binding(file, "0" * 64, Path("never-open-protocol"), self.source, self.cohort)

    def test_unsubstantiated_NA_policy_and_tampered_cohort_rejected(self):
        value = binding(self.source, self.cohort, tokens=["NA"], protocol_id="variant-transfer-v2")
        value["missing_policy"]["evidence"] = ""
        file = self.root / "binding.json"
        file.write_bytes(M.encoded(value))
        protocol = Path(__file__).resolve().parents[2] / "experiments/variant-transfer-v2.json"
        with self.assertRaisesRegex(ValueError, "requires exposure evidence"):
            M.validate_binding(file, M.sha256(file), protocol, self.source, self.cohort)
        value["missing_policy"]["evidence"] = "fixture source metadata"
        value["cohort"]["sha256"] = "0" * 64
        file.write_bytes(M.encoded(value))
        with self.assertRaisesRegex(ValueError, "cohort SHA256"):
            M.validate_binding(file, M.sha256(file), protocol, self.source, self.cohort)

    def test_v1_cannot_authorize_NA_through_modified_binding(self):
        value = binding(self.source, self.cohort, tokens=["NA"])
        file = self.root / "binding.json"
        file.write_bytes(M.encoded(value))
        protocol = Path(__file__).resolve().parents[2] / "experiments/variant-transfer-v1.json"
        with mock.patch.object(M.zipfile, "ZipFile", side_effect=AssertionError("source ZIP opened")):
            with self.assertRaisesRegex(ValueError, "v1 has no authorized missing token"):
                M.validate_binding(file, M.sha256(file), protocol, self.source, self.cohort)

    def test_v1_NA_response_remains_gate_failure_with_original_provenance(self):
        workbook(self.source, [(3, text("A", 3, "v1") + number("AN", 3, "1") + text("AO", 3, "NA"))])
        file = self.root / "binding.json"
        file.write_bytes(M.encoded(binding(self.source, self.cohort)))
        args = types.SimpleNamespace(binding=file, expected_binding_sha256=M.sha256(file), source=self.source,
            cohort=self.cohort, protocol=Path(__file__).resolve().parents[2] / "experiments/variant-transfer-v1.json", out=self.root / "strict-v1")
        with self.assertRaises(M.ObservationCellError):
            M.build(args)
        self.assertFalse((args.out / "observations.json.gz").exists())
        for name in ["failure.json", "manifest.json"]:
            self.assertEqual(M.read_json(args.out / name)["protocol_sha256"], M.PROTOCOL_SHA256)

    def test_v1_numeric_success_keeps_original_protocol_identity(self):
        file = self.root / "binding.json"
        file.write_bytes(M.encoded(binding(self.source, self.cohort)))
        args = types.SimpleNamespace(binding=file, expected_binding_sha256=M.sha256(file), source=self.source,
            cohort=self.cohort, protocol=Path(__file__).resolve().parents[2] / "experiments/variant-transfer-v1.json", out=self.root / "numeric-v1")
        M.build(args)
        for name in ["summary.json", "manifest.json"]:
            self.assertEqual(M.read_json(args.out / name)["protocol_sha256"], M.PROTOCOL_SHA256)
        payload = json.loads(gzip.decompress((args.out / "observations.json.gz").read_bytes()))
        self.assertEqual(payload["protocol_sha256"], M.PROTOCOL_SHA256)
        self.assertEqual(payload["protocol_id"], "variant-transfer-v1")
        self.assertIsNone(payload["operational_missingness"])

    def test_v2_requires_exposure_and_unknown_reason_operational_rule(self):
        original = binding(self.source, self.cohort, tokens=["NA"], protocol_id="variant-transfer-v2")
        protocol = Path(__file__).resolve().parents[2] / "experiments/variant-transfer-v2.json"
        changes = [("exposure", {}), ("exposure", {"v1_failure_observed": False}),
                   ("operational_missingness", {}),
                   ("operational_missingness", dict(M.OPERATIONAL_MISSINGNESS, source_semantics_established=True))]
        for key, changed in changes:
            with self.subTest(key=key, changed=changed):
                value = dict(original, **{key: changed})
                file = self.root / "binding.json"
                file.write_bytes(M.encoded(value))
                with self.assertRaisesRegex(ValueError, "operational missingness|prior-v1-failure"):
                    M.validate_binding(file, M.sha256(file), protocol, self.source, self.cohort)

    def test_protocol_ID_SHA_cross_binding_and_unknown_version_are_rejected(self):
        value = binding(self.source, self.cohort, tokens=["NA"], protocol_id="variant-transfer-v2")
        value["protocol_sha256"] = M.PROTOCOL_SHA256
        file = self.root / "binding.json"
        protocol = Path(__file__).resolve().parents[2] / "experiments/variant-transfer-v2.json"
        for version in ["variant-transfer-v2", "variant-transfer-v3"]:
            value["protocol_id"] = version
            file.write_bytes(M.encoded(value))
            with self.assertRaisesRegex(ValueError, "unexpected semantic protocol"):
                M.validate_binding(file, M.sha256(file), protocol, self.source, self.cohort)

    def test_full_build_records_missing_cell_and_excludes_whole_variant(self):
        self.cohort.write_bytes(M.encoded({"schema_version": 1, "pairs": [metadata(), metadata("v2", "group2")]}))
        workbook(self.source, [(3, text("A", 3, "v1") + number("AN", 3, "1") + number("AO", 3, "2")),
                               (4, text("A", 4, "v2") + text("AN", 4, "NA") + number("AO", 4, "3")),
                               (5, text("A", 5, "unselected") + number("AN", 5, "garbage", "b"))])
        file = self.root / "binding.json"
        file.write_bytes(M.encoded(binding(self.source, self.cohort, count=2, tokens=["NA"], protocol_id="variant-transfer-v2")))
        args = types.SimpleNamespace(binding=file, expected_binding_sha256=M.sha256(file), source=self.source,
            cohort=self.cohort, protocol=Path(__file__).resolve().parents[2] / "experiments/variant-transfer-v2.json", out=self.root / "out")
        M.build(args)
        summary = M.read_json(args.out / "summary.json")
        self.assertEqual((summary["metadata_eligible_variants"], summary["complete_case_variants"], summary["missing_cells"]), (2, 1, 1))
        self.assertEqual(summary["protocol_sha256"], M.PROTOCOLS["variant-transfer-v2"])
        self.assertEqual(summary["analysis_status"], "exploratory_after_v1_missingness_failure")
        self.assertFalse(summary["operational_missingness"]["source_semantics_established"])
        self.assertEqual(M.read_json(args.out / "manifest.json")["protocol_sha256"], M.PROTOCOLS["variant-transfer-v2"])
        pairs = M.read_json(args.out / "completecasepairs.json")["pairs"]
        self.assertEqual([r["id"] for r in pairs], ["v1::K562", "v1::HepG2"])
        observations = json.loads(gzip.decompress((args.out / "observations.json.gz").read_bytes()))
        self.assertEqual(observations["protocol_sha256"], M.PROTOCOLS["variant-transfer-v2"])
        self.assertEqual(len(observations["variants"]), 2)
        self.assertEqual(observations["variants"][1]["missing_contexts"], ["HepG2"])
        with gzip.open(args.out / "source-type-float-ledger.jsonl.gz", "rt") as handle:
            self.assertEqual(len(list(handle)), 4)
        for name, commitment in M.read_json(args.out / "manifest.json")["files"].items():
            self.assertEqual(M.sha256(args.out / name), commitment["sha256"])
        with self.assertRaisesRegex(ValueError, "output path exists"):
            M.build(args)

    def test_invalid_cell_failure_retains_typed_record_without_predictions(self):
        workbook(self.source, [(3, text("A", 3, "v1") + number("AN", 3, "1") + number("AO", 3, "1", "b"))])
        file = self.root / "binding.json"
        file.write_bytes(M.encoded(binding(self.source, self.cohort)))
        args = types.SimpleNamespace(binding=file, expected_binding_sha256=M.sha256(file), source=self.source,
            cohort=self.cohort, protocol=Path(__file__).resolve().parents[2] / "experiments/variant-transfer-v1.json", out=self.root / "failed-out")
        with self.assertRaises(M.ObservationCellError):
            M.build(args)
        self.assertTrue((args.out / "failure.json").exists())
        self.assertFalse((args.out / "completecasepairs.json").exists())
        with gzip.open(args.out / "source-type-float-ledger.jsonl.gz", "rt") as handle:
            record = json.loads(next(handle))
        self.assertEqual((record["context"], record["xml_type"], record["status"]), ("K562", "b", "gate_failure"))

    def test_cohort_metadata_is_validated_without_DNA_repairs(self):
        value = binding(self.source, self.cohort)
        for change in [{"sequence_ref": "a" * 200}, {"position_0based": 98}, {"gc_preserving": True}, {"component": " group"}]:
            row = metadata() | change
            self.cohort.write_bytes(M.encoded({"schema_version": 1, "pairs": [row]}))
            with self.assertRaises(ValueError):
                M.load_cohort(self.cohort, value)

    def test_zip_path_traversal_and_sheet_relationship_escape_rejected(self):
        workbook(self.source, extra_entries={"../outside": "bad"})
        with self.assertRaisesRegex(ValueError, "unsafe ZIP"):
            self.parse()


if __name__ == "__main__":
    unittest.main()
