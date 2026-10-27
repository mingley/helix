#!/usr/bin/env python3
"""Export unchanged, schema-bound source lexemes; perform no activity arithmetic."""
import argparse
import collections
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from openxml_source import SourceWorkbook as Workbook, NS  # noqa: E402


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def encode(value):
    return json.dumps(value, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def raw_rows(workbook, sheet, columns):
    """Mirror Workbook.rows while retaining selected cell types/formula/cache."""
    with workbook.zip.open(workbook.sheets[sheet]) as stream:
        for _, node in ET.iterparse(stream, events=("end",)):
            if node.tag != "{" + NS["m"] + "}row":
                continue
            cells = {}
            for cell in node.findall("m:c", NS):
                reference = cell.attrib["r"]
                column = re.fullmatch(r"([A-Z]+)[0-9]+", reference).group(1)
                if column not in columns:
                    continue
                if column in cells:
                    raise ValueError("Duplicate selected XML cell reference")
                kind = cell.attrib.get("t", "n")
                raw = cell.findtext("m:v", None, NS)
                value = raw
                if kind == "s":
                    value = workbook.shared[int(raw)] if raw is not None else None
                elif kind == "inlineStr":
                    inline = cell.find("m:is", NS)
                    value = "".join(inline.itertext()) if inline is not None else None
                formula = cell.find("m:f", NS)
                cells[column] = {"reference": reference, "xml_type": kind, "raw_v": raw, "value": value, "formula": cell.findtext("m:f", None, NS), "formula_attributes": dict(formula.attrib) if formula is not None else None}
            yield int(node.attrib["r"]), cells
            node.clear()


def classification(cell):
    if cell is None:
        return "absent_cell"
    value = cell["value"]
    if value is None:
        return "empty_cell"
    if cell["xml_type"] == "e":
        return "excel_error"
    if value == "" or value.strip() in ("NA", "N/A", "#N/A", "NaN", "nan", "NULL", "null"):
        return "missing_token"
    try:
        return "finite_numeric" if Decimal(value).is_finite() else "nonfinite_numeric"
    except (InvalidOperation, TypeError):
        return "other_text"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binding", type=Path, required=True)
    parser.add_argument("--binding-sha256", required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if sha256(args.protocol) != args.protocol_sha256 or sha256(args.binding) != args.binding_sha256:
        raise ValueError("Protocol or schema binding hash mismatch")
    binding = json.loads(args.binding.read_text())
    if binding["status"] != "bound_before_quantitative_inspection" or not binding["semantic_binding_complete"] or binding["protocol_sha256"] != args.protocol_sha256:
        raise ValueError("A previously written endpoint binding is required")
    source = binding["source_workbook"]
    path = ROOT / source["path"]
    if sha256(path) != source["sha256"] or path.stat().st_size != source["bytes"]:
        raise ValueError("Bound source workbook differs")
    base = {"schema_version": 1, "protocol_sha256": args.protocol_sha256, "schema_binding_sha256": args.binding_sha256, "source": source, "strict_reader_source": {"path": "scripts/fetch_context_panel.py", "sha256": sha256(ROOT / "scripts/fetch_context_panel.py"), "larger_source_metadata_wrapper_sha256": sha256(Path(__file__).with_name("openxml_source.py"))}}
    totals = {}
    with Workbook(path) as workbook:
        # Recheck the exact header mappings before any selected data row reads.
        for sheet, table in binding["tables"].items():
            headers = next((cells for number, cells in workbook.rows(sheet, columns=set(table["headers"]), stop_after=table["header_row"]) if number == table["header_row"]), None)
            if headers != table["headers"]:
                raise ValueError(f"Bound exact source headers differ: {sheet}")
        with args.output.open("xb") as raw_output, args.ledger.open("xb") as raw_ledger:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, mtime=0) as output, gzip.GzipFile(filename="", mode="wb", fileobj=raw_ledger, mtime=0) as ledger:
                def write(value):
                    output.write(value.encode("utf8"))
                def log(value):
                    ledger.write((encode(value) + "\n").encode("utf8"))
                write(encode(base)[:-1] + ',"tables":{')
                log({**base, "record_type": "ledger_header", "scope": "All selected numeric cell types are counted; each absent/missing/nonfinite/text/non-default-type/formula cell is retained below. Default finite numeric cells have xml_type n and raw numeric lexeme in the export."})
                first_table = True
                for sheet, table in binding["tables"].items():
                    if not first_table:
                        write(",")
                    first_table = False
                    id_column = table["id_column"]
                    value_columns = table["value_columns"]
                    key_columns = table.get("empirical_key_columns", [id_column])
                    if not set(key_columns).issubset({id_column, *value_columns}):
                        raise ValueError("Empirical key columns must be included in the export")
                    value_headers = [table["headers"][column] for column in value_columns]
                    if len(value_headers) != len(set(value_headers)):
                        raise ValueError("Duplicate selected value headers")
                    metadata = {"header_row": table["header_row"], "id_header": table["headers"][id_column], "headers": table["headers"], "selected_value_headers": value_headers, "empirical_key_headers": [table["headers"][column] for column in key_columns], "role": table["role"]}
                    write(encode(sheet) + ":" + encode(metadata)[:-1] + ',"rows":[')
                    ids = collections.Counter()
                    empirical_keys = collections.Counter()
                    counts = {"stored_data_row_elements": 0, "nonempty_id_rows": 0, "missing_id_rows": 0, "missing_empirical_key_rows": 0, "value_classifications": {header: collections.Counter() for header in value_headers}, "xml_types": {header: collections.Counter() for header in value_headers}, "formula_cells": 0}
                    first_row = True
                    for number, cells in raw_rows(workbook, sheet, {id_column, *value_columns}):
                        if number <= table["header_row"]:
                            continue
                        counts["stored_data_row_elements"] += 1
                        id_cell = cells.get(id_column)
                        identity = id_cell["value"] if id_cell is not None else None
                        if identity is not None and identity != "":
                            ids[identity] += 1
                            counts["nonempty_id_rows"] += 1
                        else:
                            counts["missing_id_rows"] += 1
                            log({"record_type": "missing_id", "sheet": sheet, "source_row": number, "cell": id_cell})
                        key = tuple(cells.get(column, {}).get("value") for column in key_columns)
                        if any(value is None or value == "" for value in key):
                            counts["missing_empirical_key_rows"] += 1
                        else:
                            empirical_keys[key] += 1
                        values = {}
                        for column in value_columns:
                            header = table["headers"][column]
                            cell = cells.get(column)
                            kind = classification(cell)
                            values[header] = cell["value"] if cell is not None else None
                            counts["value_classifications"][header][kind] += 1
                            counts["xml_types"][header][cell["xml_type"] if cell is not None else "absent"] += 1
                            formula = cell is not None and cell["formula"] is not None
                            counts["formula_cells"] += int(formula)
                            if kind != "finite_numeric" or cell["xml_type"] != "n" or formula:
                                log({"record_type": "cell_exception", "sheet": sheet, "source_row": number, "id": identity, "header": header, "column": column, "classification": kind, "cell": cell})
                        if not first_row:
                            write(",")
                        first_row = False
                        write(encode({"id": identity, "source_row": number, "values": values}))
                    counts["unique_nonempty_ids"] = len(ids)
                    counts["duplicate_id_rows_beyond_first"] = sum(number - 1 for number in ids.values())
                    counts["unique_complete_empirical_keys"] = len(empirical_keys)
                    counts["duplicate_empirical_key_rows_beyond_first"] = sum(number - 1 for number in empirical_keys.values())
                    for identity, number in ids.items():
                        if number > 1:
                            log({"record_type": "duplicate_id", "sheet": sheet, "id": identity, "count": number})
                    for key, number in empirical_keys.items():
                        if number > 1:
                            log({"record_type": "duplicate_empirical_key", "sheet": sheet, "key_headers": metadata["empirical_key_headers"], "key_values": key, "count": number})
                    write('],"counts":' + encode(counts) + "}")
                    totals[sheet] = counts
                write("}}\n")
                log({"record_type": "ledger_summary", "tables": totals})
    receipt = {**base, "implementation_sha256": sha256(__file__), "export": {"path": str(args.output), "bytes": args.output.stat().st_size, "sha256": sha256(args.output)}, "ledger": {"path": str(args.ledger), "bytes": args.ledger.stat().st_size, "sha256": sha256(args.ledger)}, "tables": totals, "activity_arithmetic_performed": False, "numeric_classification_only": "Decimal.is_finite plus missing-token/type/formula and row/ID counts; no means, variances, correlations or fitting"}
    with args.receipt.open("x") as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(encode({"export": receipt["export"], "ledger": receipt["ledger"], "table_count": len(totals)}))


if __name__ == "__main__":
    main()
