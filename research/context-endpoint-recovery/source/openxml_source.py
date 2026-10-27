"""Bounded streaming source reader reusing the existing strict XML row parser."""
from pathlib import Path
import posixpath
import sys
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from fetch_context_panel import Workbook, NS, REL  # noqa: E402


class SourceWorkbook(Workbook):
    """Same Workbook.rows decoder, with explicit bounds for the larger source.

    Only shared strings and workbook relationships are loaded into memory;
    worksheet rows stream. This does not alter measurement parsing or filtering.
    """
    def __init__(self, path, maximum_uncompressed_bytes=1_000_000_000):
        self.zip = zipfile.ZipFile(path)
        try:
            entries = self.zip.infolist()
            names = [entry.filename for entry in entries]
            self.uncompressed_bytes = sum(entry.file_size for entry in entries)
            self.maximum_uncompressed_bytes = maximum_uncompressed_bytes
            if self.uncompressed_bytes > maximum_uncompressed_bytes:
                raise ValueError("Source workbook exceeds explicit uncompressed byte bound")
            if len(names) != len(set(names)):
                raise ValueError("Duplicate workbook ZIP member names")
            if any(entry.flag_bits & 1 for entry in entries):
                raise ValueError("Encrypted workbook unsupported")
            self.shared = []
            if "xl/sharedStrings.xml" in names:
                if self.zip.getinfo("xl/sharedStrings.xml").file_size > 128_000_000:
                    raise ValueError("Shared-string XML exceeds memory byte bound")
                tree = ET.fromstring(self.zip.read("xl/sharedStrings.xml"))
                self.shared = ["".join(node.itertext()) for node in tree.findall("m:si", NS)]
            book = ET.fromstring(self.zip.read("xl/workbook.xml"))
            relationships = ET.fromstring(self.zip.read("xl/_rels/workbook.xml.rels"))
            targets = {node.attrib["Id"]: node.attrib["Target"] for node in relationships if node.attrib.get("TargetMode") != "External"}
            self.sheets = {}
            for sheet in book.findall("m:sheets/m:sheet", NS):
                target = targets[sheet.attrib["{" + REL + "}id"]]
                member = posixpath.normpath(target.lstrip("/") if target.startswith("/") else "xl/" + target)
                if not member.startswith("xl/") or ".." in member.split("/") or "\\" in member:
                    raise ValueError("Workbook sheet relationship escapes archive")
                if sheet.attrib["name"] in self.sheets:
                    raise ValueError("Duplicate worksheet names")
                self.sheets[sheet.attrib["name"]] = member
        except BaseException:
            self.zip.close()
            raise
