#!/usr/bin/env python3
"""Strict XML observations from a separately sealed Siraj S4/cohort binding.

No model operations. Source responses are decoded only for exact metadata-frozen
cohort IDs, after the caller-pinned binding and all prerequisite hashes pass.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import posixpath
import re
import shutil
import stat
import struct
import sys
import xml.etree.ElementTree as ET
import zipfile

PROTOCOL_SHA256 = "811930c2d3a5ee763014ea611acda4d2c84404b7d5ae399c7843c62d0a828149"
MODEL_SHA256 = "2ea476bc5413dc13486dec9f3d7e9162a50a0e47f1d83ce29efabc5f5cc5af9d"
CONTEXTS = ["K562", "HepG2"]
SHEET = "S4.MPRA results"
HEADERS = {"A": "Variant", "AN": "log2Skew in HEPG2", "AO": "log2Skew in K562"}
EFFECTS = {"HepG2": "AN", "K562": "AO"}
MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS = {"m": MAIN}
HEX = re.compile(r"[0-9a-f]{64}\Z")
CELL = re.compile(r"([A-Z]+)([1-9][0-9]*)\Z")
NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")
MAX_ARCHIVE_BYTES = 400_000_000
MAX_UNCOMPRESSED_BYTES = 2_000_000_000
MAX_MEMBERS = 10_000


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def unique_object(items):
    result = {}
    for key, value in items:
        require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def decode_json(text):
    return json.loads(text, object_pairs_hook=unique_object,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def read_json(path):
    return decode_json(Path(path).read_text(encoding="utf-8"))


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def write_gzip(path, value):
    with Path(path).open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, filename="", mode="wb", mtime=0) as handle:
            handle.write(encoded(value))


def write_jsonl(path, rows):
    with Path(path).open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, filename="", mode="wb", mtime=0) as handle:
            for row in rows:
                handle.write((json.dumps(row, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode())


def label(value):
    return isinstance(value, str) and bool(value) and value.strip() == value and "\n" not in value


def valid_digest(value):
    return isinstance(value, str) and bool(HEX.fullmatch(value))


def validate_binding(binding_path, expected_sha256, protocol_path, source_path, cohort_path):
    """Never open source ZIP or cohort payload before the sealed gate passes."""
    require(valid_digest(expected_sha256) and sha256(binding_path) == expected_sha256,
            "binding SHA256 mismatch")
    binding = read_json(binding_path)
    require(type(binding.get("schema_version")) is int and binding["schema_version"] == 1,
            "unexpected binding schema")
    require(binding.get("protocol_id") == "variant-transfer-v1"
            and binding.get("protocol_sha256") == PROTOCOL_SHA256, "unexpected semantic protocol binding")
    require(isinstance(binding.get("source_cohort_gate"), dict)
            and binding["source_cohort_gate"].get("passed") is True, "sealed source/cohort gate failed")
    require(binding.get("model_artifacts_sha256") == MODEL_SHA256, "frozen model binding mismatch")
    require(sha256(protocol_path) == PROTOCOL_SHA256, "semantic protocol bytes changed")
    sheet = binding.get("sheet", {})
    require(sheet.get("name") == SHEET and sheet.get("header_row_1based") == 2
            and type(sheet.get("header_row_1based")) is int and sheet.get("headers") == HEADERS,
            "exact source sheet/header binding differs")
    require(binding.get("effects") == EFFECTS, "cell-to-column binding differs")
    policy = binding.get("missing_policy", {})
    tokens = policy.get("allowed_tokens")
    require(tokens in ([], ["NA"]), "only a source-documented exact NA token may be allowed")
    require(not tokens or bool(policy.get("evidence")), "missing token requires source evidence")
    source = binding.get("source", {})
    require(valid_digest(source.get("sha256")) and type(source.get("bytes")) is int
            and 0 < source["bytes"] <= MAX_ARCHIVE_BYTES, "invalid source byte commitment")
    cohort = binding.get("cohort", {})
    require(valid_digest(cohort.get("sha256")) and cohort.get("format") in {"json", "jsonl.gz"}
            and type(cohort.get("variant_count")) is int and cohort["variant_count"] > 0,
            "invalid metadata cohort commitment")
    require(not Path(source_path).is_symlink() and Path(source_path).stat().st_size == source["bytes"]
            and sha256(source_path) == source["sha256"], "raw S4 source byte commitment failed")
    require(not Path(cohort_path).is_symlink() and sha256(cohort_path) == cohort["sha256"],
            "metadata cohort SHA256 mismatch")
    return binding


def load_cohort(path, binding):
    if binding["cohort"]["format"] == "jsonl.gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            rows = [decode_json(line) for line in handle]
    else:
        payload = read_json(path)
        require(type(payload.get("schema_version")) is int and payload["schema_version"] == 1,
                "invalid cohort schema")
        require(set(payload) == {"schema_version", "pairs"}, "expected schema_version/pairs cohort object")
        rows = payload["pairs"]
    require(isinstance(rows, list) and len(rows) == binding["cohort"]["variant_count"],
            "cohort count differs from second freeze")
    result = {}
    for row in rows:
        require(isinstance(row, dict) and all(label(row.get(key)) for key in ("id", "component", "locus")),
                "invalid cohort identity/component/locus")
        require(row["id"] not in result, "duplicate metadata cohort ID")
        require("::" not in row["id"], "variant ID conflicts with context-pair delimiter")
        for name in ("sequence_ref", "sequence_alt"):
            sequence = row.get(name)
            require(isinstance(sequence, str) and len(sequence) == 200 and not set(sequence) - set("ACGT"),
                    "invalid fixed 200nt cohort DNA")
        changes = [index for index, (a, b) in enumerate(zip(row["sequence_ref"], row["sequence_alt"])) if a != b]
        require(changes == [99] and type(row.get("position_0based", 99)) is int
                and row.get("position_0based", 99) == 99, "cohort is not a centered same-length SNV")
        preserving = (row["sequence_ref"][99] in "GC") == (row["sequence_alt"][99] in "GC")
        require(type(row.get("gc_preserving")) is bool and row["gc_preserving"] == preserving,
                "cohort GC metadata contradicts DNA")
        result[row["id"]] = {"id": row["id"], "component": row["component"], "locus": row["locus"],
                              "reference": row["sequence_ref"], "alternate": row["sequence_alt"],
                              "position_0based": 99, "gc_preserving": preserving}
    return result


def validate_archive(archive):
    entries = archive.infolist()
    names = [entry.filename for entry in entries]
    require(len(entries) <= MAX_MEMBERS and len(names) == len(set(names)), "ZIP member limit or duplicate name")
    require(sum(entry.file_size for entry in entries) <= MAX_UNCOMPRESSED_BYTES,
            "ZIP total uncompressed byte budget exceeded")
    observed = 0
    for entry in entries:
        part = PurePosixPath(entry.filename)
        require(not part.is_absolute() and ".." not in part.parts and "\\" not in entry.filename
                and str(part) == entry.filename.rstrip("/"), "unsafe ZIP member name")
        require(not entry.flag_bits & 1 and not stat.S_ISLNK(entry.external_attr >> 16),
                "encrypted or symlink ZIP member")
        count = 0
        with archive.open(entry) as handle:
            while block := handle.read(1024 * 1024):
                count += len(block)
                observed += len(block)
                require(count <= entry.file_size and observed <= MAX_UNCOMPRESSED_BYTES,
                        "actual ZIP data exceeded byte budget")
        require(count == entry.file_size, "ZIP member decompressed length mismatch")
    return {"members": len(entries), "uncompressed_bytes": observed, "all_member_crc_checked": True}


def source_sheet(archive):
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {}
    for relation in relationships:
        require(relation.attrib["Id"] not in targets, "duplicate workbook relationship ID")
        if relation.attrib.get("TargetMode") != "External":
            targets[relation.attrib["Id"]] = relation.attrib["Target"]
    sheets = workbook.findall("m:sheets/m:sheet", NS)
    require(len(sheets) == 1 and sheets[0].attrib.get("name") == SHEET, "workbook sheet identity differs")
    target = targets[sheets[0].attrib["{" + REL + "}id"]]
    member = posixpath.normpath(target.lstrip("/") if target.startswith("/") else "xl/" + target)
    require(member.startswith("xl/") and ".." not in member.split("/") and "\\" not in member,
            "unsafe sheet relationship target")
    require(member in archive.namelist(), "worksheet target absent")
    return member


def shared_strings(archive):
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    result = []
    with archive.open("xl/sharedStrings.xml") as handle:
        for event, element in ET.iterparse(handle, events=("end",)):
            if element.tag == "{" + MAIN + "}si":
                result.append("".join(node.text or "" for node in element.iter("{" + MAIN + "}t")))
                element.clear()
    return result


def text_cell(cell, shared):
    require(cell is not None and cell.find("m:f", NS) is None, "ID/header formula or absent cell")
    require(len(cell.findall("m:v", NS)) <= 1 and len(cell.findall("m:is", NS)) <= 1,
            "duplicate ID/header cell payload")
    kind = cell.attrib.get("t", "n")
    if kind == "s":
        raw = cell.findtext("m:v", default=None, namespaces=NS)
        require(raw is not None and raw.isdecimal() and int(raw) < len(shared), "invalid shared-string reference")
        return shared[int(raw)]
    if kind == "inlineStr":
        inline = cell.find("m:is", NS)
        require(inline is not None, "missing inline text")
        return "".join(node.text or "" for node in inline.iter("{" + MAIN + "}t"))
    if kind == "str":
        raw = cell.findtext("m:v", default=None, namespaces=NS)
        require(raw is not None, "missing string value")
        return raw
    raise ValueError("ID/header is not an ordinary text cell")


class ObservationCellError(ValueError):
    def __init__(self, message, record):
        super().__init__(message)
        self.record = record


def effect_cell(cell, shared, allowed_tokens):
    """Returns raw type/lexeme and exact binary64 encoding; no coercion repairs."""
    kind = cell.attrib.get("t", "n") if cell is not None else None
    raw = cell.findtext("m:v", default=None, namespaces=NS) if cell is not None else None
    record = {"coordinate": cell.attrib.get("r") if cell is not None else None,
              "xml_type": kind, "raw_lexeme": raw,
              "formula": cell.findtext("m:f", default=None, namespaces=NS) if cell is not None else None}
    try:
        require(cell is not None, "unexplained absent effect cell")
        require(len(cell.findall("m:v", NS)) <= 1 and len(cell.findall("m:is", NS)) <= 1,
                "duplicate effect cell payload")
        require(cell.find("m:f", NS) is None, "effect formula is not an observation")
        if kind in {"s", "inlineStr", "str"}:
            text = text_cell(cell, shared)
            record["text"] = text
            require(text in allowed_tokens, "unexplained text or numeric string in effect cell")
            record.update(status="source_documented_missing", value=None, float64_le_hex=None)
            return record
        require(kind == "n", "effect is an error, boolean, date or unsupported type")
        require(raw is not None and bool(NUMBER.fullmatch(raw)), "unexplained blank or invalid numeric literal")
        value = float(raw)
        require(math.isfinite(value), "effect numeric literal is nonfinite")
        record.update(status="ordinary_finite_numeric", value=value,
                      float64_le_hex=struct.pack("<d", value).hex())
        return record
    except ValueError as error:
        record.update(status="gate_failure", reason=str(error))
        raise ObservationCellError(str(error), record) from error


def physical_rows(archive, member):
    with archive.open(member) as handle:
        stack = []
        for event, element in ET.iterparse(handle, events=("start", "end")):
            if event == "start":
                stack.append(element)
            else:
                if element.tag == "{" + MAIN + "}row":
                    yield element
                    element.clear()
                    if len(stack) >= 2:
                        stack[-2].remove(element)
                stack.pop()


def row_cells(row, row_number):
    cells = {}
    for cell in row.findall("m:c", NS):
        coordinate = CELL.fullmatch(cell.attrib.get("r", ""))
        require(coordinate is not None and int(coordinate[2]) == row_number, "cell coordinate/row mismatch")
        column = coordinate[1]
        require(column not in cells, "duplicate cell coordinate")
        cells[column] = cell
    return cells


def parse_source(path, cohort, allowed_tokens, id_ledger, effect_ledger):
    selected, seen, row_numbers = {}, set(), set()
    header_seen = False
    with zipfile.ZipFile(path) as archive:
        integrity = validate_archive(archive)
        member = source_sheet(archive)
        shared = shared_strings(archive)
        for row in physical_rows(archive, member):
            raw_number = row.attrib.get("r", "")
            require(raw_number.isdecimal() and int(raw_number) > 0, "invalid physical row number")
            number = int(raw_number)
            require(number not in row_numbers, "duplicate physical row number")
            row_numbers.add(number)
            if number < 2:
                continue
            cells = row_cells(row, number)
            if number == 2:
                for column, expected in HEADERS.items():
                    require(text_cell(cells.get(column), shared) == expected, "source header differs: " + column)
                header_seen = True
                continue
            require(header_seen, "source data precedes frozen header")
            if not cells:
                id_ledger.append({"row_1based": number, "id": None, "status": "physically_empty_row"})
                continue
            identity = text_cell(cells.get("A"), shared)
            require(label(identity), "invalid source Variant ID")
            require(identity not in seen, "duplicate source Variant ID")
            seen.add(identity)
            id_ledger.append({"row_1based": number, "id": identity,
                              "selected_by_frozen_metadata": identity in cohort})
            if identity not in cohort:
                continue
            observations = {}
            for context in CONTEXTS:
                try:
                    record = effect_cell(cells.get(EFFECTS[context]), shared, allowed_tokens)
                except ObservationCellError as error:
                    error.record.update(id=identity, context=context, row_1based=number)
                    effect_ledger.append(error.record)
                    raise
                record.update(id=identity, context=context, row_1based=number)
                effect_ledger.append(record)
                observations[context] = record
            selected[identity] = observations
    require(header_seen, "frozen header row absent")
    require(set(selected) == set(cohort), "metadata cohort ID absent from source table")
    return selected, integrity, {"source_unique_IDs": len(seen), "physical_rows": len(row_numbers),
                                  "metadata_selected_IDs": len(selected)}


def build(args):
    require(not args.out.exists(), "output path exists; preserve all attempts")
    binding = validate_binding(args.binding, args.expected_binding_sha256, args.protocol, args.source, args.cohort)
    cohort = load_cohort(args.cohort, binding)
    args.out.mkdir(parents=True)
    shutil.copyfile(Path(__file__), args.out / "implementation.py")
    shutil.copyfile(args.binding, args.out / "binding.json")
    shutil.copyfile(args.protocol, args.out / "protocol.json")
    ids, cells = [], []
    try:
        selected, integrity, inventory = parse_source(args.source, cohort,
            binding["missing_policy"]["allowed_tokens"], ids, cells)
        variants, pairs, missing = [], [], []
        for identity in sorted(cohort):
            metadata, observations = cohort[identity], selected[identity]
            unavailable = [context for context in CONTEXTS if observations[context]["value"] is None]
            complete = not unavailable
            variant = dict(metadata, complete_case=complete,
                observed={context: observations[context]["value"] for context in CONTEXTS},
                source_float64_le_hex={context: observations[context]["float64_le_hex"] for context in CONTEXTS},
                missing_contexts=unavailable)
            variants.append(variant)
            for context in unavailable:
                missing.append({"id": identity, "component": metadata["component"], "context": context,
                                **observations[context]})
            if complete:
                for context in CONTEXTS:
                    pairs.append({"id": identity + "::" + context, "context": context,
                                  "position_0based": 99, "reference": metadata["reference"],
                                  "alternate": metadata["alternate"]})
        payload = {"schema_version": 1, "protocol_sha256": PROTOCOL_SHA256,
                   "binding_sha256": args.expected_binding_sha256, "variants": variants}
        write_gzip(args.out / "observations.json.gz", payload)
        (args.out / "completecasepairs.json").write_bytes(encoded({"schema_version": 1, "pairs": pairs}))
        write_jsonl(args.out / "missingness.jsonl.gz", missing)
        component_counts = Counter(row["component"] for row in variants if row["complete_case"])
        summary = {"schema_version": 1, "protocol_sha256": PROTOCOL_SHA256,
                   "binding_sha256": args.expected_binding_sha256, "model_artifacts_sha256": MODEL_SHA256,
                   "source_sha256": binding["source"]["sha256"], "source_bytes": binding["source"]["bytes"],
                   "cohort_sha256": binding["cohort"]["sha256"], "source_cohort_gate": binding["source_cohort_gate"],
                   "source_observation_gate_passed": True, "inventory": inventory, "archive_integrity": integrity,
                   "metadata_eligible_variants": len(variants), "complete_case_variants": len(pairs) // 2,
                   "complete_case_components": len(component_counts), "missing_variants": sum(not row["complete_case"] for row in variants),
                   "missing_cells": len(missing), "missing_cells_by_context": dict(Counter(row["context"] for row in missing)),
                   "missing_cells_by_component": dict(Counter(row["component"] for row in missing)),
                   "minimum_support_gate_passed": len(pairs) // 2 >= 100 and len(component_counts) >= 10,
                   "models_fitted": False, "predictions_computed": False,
                   "limits": "Literal published allelic-effect cells only; complete-case rule changes no metadata eligibility or grouping."}
        (args.out / "summary.json").write_bytes(encoded(summary))
        print(json.dumps({key: summary[key] for key in ("metadata_eligible_variants", "complete_case_variants", "missing_cells")}, sort_keys=True))
    except BaseException as error:
        (args.out / "failure.json").write_bytes(encoded({"error_type": type(error).__name__, "error": str(error),
            "binding_sha256": args.expected_binding_sha256, "source_observation_gate_passed": False}))
        raise
    finally:
        write_jsonl(args.out / "source-IDs.jsonl.gz", ids)
        write_jsonl(args.out / "source-type-float-ledger.jsonl.gz", cells)
        (args.out / "manifest.json").write_bytes(encoded({"schema_version": 1,
            "protocol_sha256": PROTOCOL_SHA256, "binding_sha256": args.expected_binding_sha256,
            "source_sha256": binding["source"]["sha256"], "cohort_sha256": binding["cohort"]["sha256"],
            "files": {path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
                      for path in sorted(args.out.iterdir()) if path.name != "manifest.json"}}))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binding", type=Path, required=True)
    parser.add_argument("--expected-binding-sha256", required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, default=Path("experiments/variant-transfer-v1.json"))
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args(argv))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, TypeError, OSError, zipfile.BadZipFile, ET.ParseError) as error:
        print(f"variant observation import failed: {error}", file=sys.stderr)
        raise SystemExit(2)
