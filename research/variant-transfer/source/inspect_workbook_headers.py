#!/usr/bin/env python3
"""CRC-validate an XLSX and retain only its first two physical rows as schema.

Non-string cell values are redacted without conversion or interpretation.
The workbook's data rows are never decoded. ZIP byte scanning is integrity
checking rather than authorization to interpret scientific values.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import posixpath
import sys
import xml.etree.ElementTree as E
import zipfile

ROOT = Path(__file__).resolve().parents[3]
S = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def inspect(path, out, sheet_prefix=None):
    protocol_path = ROOT / "experiments/variant-transfer-v1.json"
    protocol_sha = digest(protocol_path)
    if protocol_sha != "811930c2d3a5ee763014ea611acda4d2c84404b7d5ae399c7843c62d0a828149":
        raise ValueError("frozen semantic protocol hash mismatch")
    protocol = json.loads(protocol_path.read_text())
    limits = protocol["source"]["acquisition_limits"]
    helper = ROOT / "research/context-endpoint-recovery/source/verify_archive.py"
    spec = importlib.util.spec_from_file_location("source_archive_verifier", helper)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with zipfile.ZipFile(path) as z:
        if len(z.infolist()) > limits["max_zip_members"]:
            raise ValueError("ZIP member count exceeds frozen limit")
        xml_bytes = sum(i.file_size for i in z.infolist() if i.filename.endswith(".xml"))
        if xml_bytes > limits["max_total_xml_bytes"]:
            raise ValueError("XML bytes exceed frozen limit")
    integrity = module.inspect_and_scan(path, max_member_bytes=limits["max_total_uncompressed_bytes"],
                                       max_total_bytes=limits["max_total_uncompressed_bytes"])
    out.mkdir(parents=True, exist_ok=False)
    archive_sha = digest(path)
    (out / "archive-validation.json").write_text(json.dumps({
        "protocol_sha256": protocol_sha, "archive_sha256": archive_sha,
        "archive_bytes": path.stat().st_size, "helper_sha256": digest(helper),
        "total_xml_bytes": xml_bytes, "inspection": integrity,
        "scientific_outcome_values_decoded": False,
    }, indent=2) + "\n")
    with zipfile.ZipFile(path) as z:
        workbook = E.fromstring(z.read("xl/workbook.xml"))
        rels = E.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        targets = {}
        for rel in rels:
            if rel.attrib.get("TargetMode") == "External":
                raise ValueError("external workbook relationship")
            target = rel.attrib["Target"]
            target = posixpath.normpath(target.lstrip("/") if target.startswith("/")
                                         else posixpath.join("xl", target))
            if not target.startswith("xl/") or ".." in PurePosixPath(target).parts:
                raise ValueError("unsafe workbook relationship")
            targets[rel.attrib["Id"]] = target
        wanted = set()
        sheets = []
        selected_sheets = [sheet for sheet in workbook.find(S + "sheets")
                           if sheet_prefix is None or sheet.attrib["name"].startswith(sheet_prefix)]
        if sheet_prefix is not None and len(selected_sheets) != 1:
            raise ValueError("metadata sheet prefix does not identify exactly one worksheet")
        for sheet in selected_sheets:
            member = targets[sheet.attrib[R + "id"]]
            rows = []
            with z.open(member) as source:
                for _, row in E.iterparse(source, events=("end",)):
                    if row.tag != S + "row":
                        continue
                    number = int(row.attrib["r"])
                    if number > 2:
                        break
                    cells = []
                    for cell in row.findall(S + "c"):
                        typ = cell.attrib.get("t", "n")
                        item = {"coordinate": cell.attrib["r"], "type": typ}
                        if typ == "s":
                            index = int(cell.findtext(S + "v"))
                            wanted.add(index)
                            item["shared_string_index"] = index
                        elif typ == "inlineStr":
                            item["text"] = "".join(t.text or "" for t in cell.findall(".//" + S + "t"))
                        else:
                            item["value"] = "[redacted non-text cell; not interpreted]"
                        cells.append(item)
                    rows.append({"row_1based": number, "cells": cells})
                    row.clear()
            sheets.append({"sheet": sheet.attrib["name"], "member": member, "schema_rows": rows})
        strings = {}
        if wanted:
            with z.open("xl/sharedStrings.xml") as source:
                index = 0
                for _, node in E.iterparse(source, events=("end",)):
                    if node.tag != S + "si":
                        continue
                    if index in wanted:
                        strings[index] = "".join(t.text or "" for t in node.findall(".//" + S + "t"))
                    index += 1
                    node.clear()
                    if wanted <= strings.keys():
                        break
            if wanted != strings.keys():
                raise ValueError("missing selected shared string")
        for sheet in sheets:
            for row in sheet["schema_rows"]:
                for cell in row["cells"]:
                    if "shared_string_index" in cell:
                        cell["text"] = strings[cell["shared_string_index"]]
        properties = {i.filename: z.read(i).decode("utf-8") for i in z.infolist()
                      if i.filename.startswith("docProps/") and i.filename.endswith(".xml")}
    result = {"schema_version": 1, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
              "protocol_sha256": protocol_sha, "archive_sha256": archive_sha,
              "implementation_sha256": digest(__file__), "scope": "first two physical rows only; no numeric scientific values decoded",
              "sheets": sheets, "document_properties": properties}
    (out / "headers-before-values.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--sheet-prefix", help="inspect exactly one matching metadata sheet; no other sheet rows decoded")
    args = p.parse_args()
    inspect(args.source, args.out, args.sheet_prefix)


if __name__ == "__main__":
    main()
