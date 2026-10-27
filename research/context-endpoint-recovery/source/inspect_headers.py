#!/usr/bin/env python3
"""Inspect workbook names and first-row text headers, redacting numeric cells."""
import argparse
import datetime as dt
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from extract_members import sha256
from openxml_source import SourceWorkbook, NS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--workbook-sha256", required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if sha256(args.workbook) != args.workbook_sha256 or sha256(args.protocol) != args.protocol_sha256:
        raise ValueError("Protocol or source workbook hash mismatch")
    result = {"schema_version": 1, "inspected_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "protocol_sha256": args.protocol_sha256, "source_workbook": {"path": str(args.workbook), "sha256": args.workbook_sha256, "bytes": args.workbook.stat().st_size}, "scope": "Workbook metadata and first stored row text only. Numeric cell values redacted. No data rows/quantitative activity values interpreted.", "implementation_sha256": sha256(__file__), "reader_sha256": sha256(Path(__file__).with_name("openxml_source.py")), "sheets": {}}
    with SourceWorkbook(args.workbook) as workbook:
        result["workbook_uncompressed_bytes"] = workbook.uncompressed_bytes
        result["workbook_uncompressed_byte_bound"] = workbook.maximum_uncompressed_bytes
        for sheet, member in workbook.sheets.items():
            info = workbook.zip.getinfo(member)
            entry = {"xml_member": member, "xml_uncompressed_bytes": info.file_size, "sheet_dimension_not_a_record_count": None, "header_row": None, "headers": {}, "numeric_header_cells_redacted": []}
            with workbook.zip.open(member) as stream:
                for _, node in ET.iterparse(stream, events=("end",)):
                    if node.tag == "{" + NS["m"] + "}dimension":
                        entry["sheet_dimension_not_a_record_count"] = node.attrib.get("ref")
                    if node.tag != "{" + NS["m"] + "}row":
                        continue
                    row_number = int(node.attrib["r"])
                    entry["header_row"] = row_number
                    text_columns = set()
                    for cell in node.findall("m:c", NS):
                        column = re.fullmatch(r"([A-Z]+)[0-9]+", cell.attrib["r"]).group(1)
                        if cell.attrib.get("t", "n") in ("s", "str", "inlineStr"):
                            text_columns.add(column)
                        else:
                            entry["numeric_header_cells_redacted"].append(column)
                    # Existing strict row decoder is limited to text header columns.
                    entry["headers"] = next(cells for number, cells in workbook.rows(sheet, columns=text_columns, stop_after=row_number) if number == row_number)
                    break
            result["sheets"][sheet] = entry
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"workbook_uncompressed_bytes": result["workbook_uncompressed_bytes"], "sheets": {sheet: value["headers"] for sheet, value in result["sheets"].items()}}, indent=2))


if __name__ == "__main__":
    main()
