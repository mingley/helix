"""Synthetic archive-integrity tests. No network or empirical workbook inputs."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import stat
import struct
import tempfile
import unittest
import warnings
import zipfile
import zlib

SPEC=importlib.util.spec_from_file_location("context_archive_integrity",Path(__file__).with_name("verify_archive.py"))
V=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


class ArchiveIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root=Path(self.temporary.name)
        self.path=self.root/"synthetic.zip"
        self.make_zip()

    def make_zip(self,entries=None,compression=zipfile.ZIP_DEFLATED):
        with zipfile.ZipFile(self.path,"w",compression=compression)as archive:
            for name,payload in entries or [("tables/",b""),("tables/one.xlsx",b"synthetic one"),("tables/two.xlsx",b"synthetic two")]:
                archive.writestr(name,payload)

    def bound_inputs(self,mime="application/zip"):
        contract={"authoritative_source":{"primary_archive":{
            "url":"synthetic://no-network-source","expected_size_bytes":self.path.stat().st_size}}}
        protocol=self.root/"protocol-input.json"
        protocol.write_text(json.dumps(contract))
        protocol_sha=V.sha256(protocol)
        acquired={"protocol_sha256":protocol_sha,"url":contract["authoritative_source"]["primary_archive"]["url"],
            "complete":True,"status":200,"bytes":self.path.stat().st_size,"sha256":V.sha256(self.path),
            "headers":{"content-type":mime}}
        acquisition=self.root/"acquisition-input.json"
        acquisition.write_text(json.dumps(acquired))
        return protocol,protocol_sha,acquisition

    def test_all_members_crc_and_hash_are_verified_including_unselected_members(self):
        report=V.inspect_and_scan(self.path,["tables/one.xlsx"])
        self.assertTrue(report["all_members_crc_verified"])
        self.assertEqual(report["member_count"],3)
        self.assertEqual(report["total_uncompressed_bytes"],len(b"synthetic onesynthetic two"))
        self.assertTrue(all(member["crc_verified"]for member in report["members"]))
        self.assertEqual(sum(member["selected_for_authorized_schema_inspection"]for member in report["members"]),1)
        self.assertEqual(report["members"][2]["sha256"],hashlib.sha256(b"synthetic two").hexdigest())

    def test_traversal_absolute_duplicate_and_symlink_names_fail(self):
        for name in ["../escape.xlsx","/absolute.xlsx","drive:/x.xlsx","bad\\x.xlsx","bad//x.xlsx","bad/./x.xlsx"]:
            if name=="drive:/x.xlsx":
                name="C:/x.xlsx"
            with self.subTest(name=name):
                self.make_zip([(name,b"x")])
                with self.assertRaises(ValueError):V.inspect_and_scan(self.path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",UserWarning)
            self.make_zip([("same.xlsx",b"one"),("same.xlsx",b"two")])
        with self.assertRaisesRegex(ValueError,"duplicate"):V.inspect_and_scan(self.path)
        link=zipfile.ZipInfo("link.xlsx");link.create_system=3
        link.external_attr=(stat.S_IFLNK|0o777)<<16
        self.make_zip([(link,b"target")])
        with self.assertRaisesRegex(ValueError,"nonregular"):V.inspect_and_scan(self.path)

    def test_html_prefix_truncation_and_trailing_bytes_fail(self):
        raw=self.path.read_bytes()
        for changed in [b"<!doctype html>"+raw,raw[:-1],raw+b"unexpected"]:
            self.path.write_bytes(changed)
            with self.assertRaises((ValueError,zipfile.BadZipFile)):V.inspect_and_scan(self.path)

    def test_declared_budget_and_missing_selection_fail_before_crc_scan(self):
        with self.assertRaisesRegex(ValueError,"budget"):V.inspect_and_scan(self.path,max_member_bytes=2)
        with self.assertRaisesRegex(ValueError,"budget"):V.inspect_and_scan(self.path,max_total_bytes=20)
        with self.assertRaisesRegex(ValueError,"absent"):V.inspect_and_scan(self.path,["absent.xlsx"])
        with self.assertRaisesRegex(ValueError,"duplicate selected"):V.inspect_and_scan(self.path,["tables/one.xlsx"]*2)

    def test_out_of_bounds_local_header_is_rejected_without_payload_interpretation(self):
        raw=bytearray(self.path.read_bytes());central=raw.index(b"PK\x01\x02")
        struct.pack_into("<I",raw,central+42,len(raw)+999)
        self.path.write_bytes(raw)
        with self.assertRaisesRegex(ValueError,"outside member area"):V.inspect_and_scan(self.path)

    def test_central_local_filename_disagreement_is_rejected(self):
        raw=bytearray(self.path.read_bytes())
        with zipfile.ZipFile(self.path)as archive:offset=archive.infolist()[1].header_offset
        raw[offset+30]=ord("x");self.path.write_bytes(raw)
        with self.assertRaisesRegex(ValueError,"names disagree"):V.inspect_and_scan(self.path)

    def test_corrupt_unselected_member_crc_fails_full_scan(self):
        self.make_zip(compression=zipfile.ZIP_STORED)
        raw=self.path.read_bytes().replace(b"synthetic two",b"synthetic twX")
        self.path.write_bytes(raw)
        with self.assertRaises((ValueError,zipfile.BadZipFile,zlib.error)):
            V.inspect_and_scan(self.path,["tables/one.xlsx"])

    def test_data_descriptors_are_supported_and_checked(self):
        class NonSeeking(io.BytesIO):
            def seek(self,*_args):raise io.UnsupportedOperation("nonseekable")
        stream=NonSeeking()
        with zipfile.ZipFile(stream,"w",compression=zipfile.ZIP_DEFLATED)as archive:
            archive.writestr("one.xlsx",b"descriptor payload")
        self.path.write_bytes(stream.getvalue())
        self.assertTrue(V.inspect_and_scan(self.path)["all_members_crc_verified"])
        raw=bytearray(self.path.read_bytes());position=raw.index(b"PK\x07\x08");raw[position+4]^=1
        self.path.write_bytes(raw)
        with self.assertRaisesRegex(ValueError,"descriptor"):V.inspect_and_scan(self.path)

    def test_result_manifest_binds_receipts_inventory_and_implementation_snapshot(self):
        protocol,protocol_sha,acquisition=self.bound_inputs()
        out=self.root/"result"
        result=V.verify(self.path,V.sha256(self.path),protocol,protocol_sha,acquisition,out,
            ["tables/one.xlsx"],predicted_archive_sha256="0"*64)
        self.assertTrue(result["complete"])
        self.assertFalse(result["producer_checksum_asserted"])
        self.assertFalse(result["cross_source_code_prediction"]["matches_official_download"])
        manifest=json.loads((out/"manifest.json").read_bytes())
        self.assertEqual(set(manifest["files"]),{"implementation.py","verification.json","inventory.json","protocol.json","acquisition.json"})
        for name,expected in manifest["files"].items():self.assertEqual(V.sha256(out/name),expected)
        self.assertEqual(result["implementation_sha256"],V.sha256(out/"implementation.py"))
        with self.assertRaises(FileExistsError):
            V.verify(self.path,V.sha256(self.path),protocol,protocol_sha,acquisition,out)

    def test_wrong_mime_creates_bound_failure_receipt_and_no_success_inventory(self):
        protocol,protocol_sha,acquisition=self.bound_inputs(mime="text/html")
        out=self.root/"failure"
        with self.assertRaisesRegex(ValueError,"MIME"):
            V.verify(self.path,V.sha256(self.path),protocol,protocol_sha,acquisition,out)
        receipt=json.loads((out/"verification.json").read_bytes())
        self.assertFalse(receipt["complete"])
        self.assertEqual(receipt["status"],"failed")
        self.assertFalse((out/"inventory.json").exists())
        self.assertTrue((out/"manifest.json").exists())


if __name__=="__main__":unittest.main()
