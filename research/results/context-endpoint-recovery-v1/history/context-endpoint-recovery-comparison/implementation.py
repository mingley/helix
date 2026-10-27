#!/usr/bin/env python3
"""Independent openpyxl inventory and endpoint helpers for a future frozen source.

The command line currently inventories explicitly checksum-bound local XLSX
files only. It never downloads quantitative payloads, fits a model, changes an
old target or selects an endpoint. Empirical endpoint comparisons await the
centrally frozen payload, source schema, equations and ID-mapping contract.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import datetime
import gzip
from decimal import Decimal
from fractions import Fraction
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import re
import shutil
import traceback
import zipfile

import numpy as np
import openpyxl
from openpyxl.worksheet._reader import WorkSheetParser
from openpyxl.utils.cell import coordinate_to_tuple
from openpyxl.utils.cell import column_index_from_string


PINNED_OPENPYXL = "3.1.5"
FROZEN_RECOVERY_SHA256 = "00175e41cf9fa0cbf177ce8fc7ef143eb508aec8b11d058f7acb4458284dc868"
VALUE_TAG = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}v"
FORMULA_TAG = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}f"


def hashes(path: Path) -> dict:
    sha, md5 = hashlib.sha256(), hashlib.md5(usedforsecurity=False)
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            size += len(block)
            sha.update(block)
            md5.update(block)
    return {"bytes": size, "sha256": sha.hexdigest(), "md5": md5.hexdigest()}


@dataclass(frozen=True)
class SourceCell:
    row: int
    column: int
    kind: str
    value: object
    xml_type: str
    raw_v: str | None
    formula_present: bool
    formula: str | None

    @property
    def meaningful(self) -> bool:
        return self.formula_present or (self.value is not None and self.value != "")

    def numeric_fraction(self, *, allow_cached_formula: bool = False) -> Fraction:
        """Read a source literal exactly; never treat booleans/errors as numbers."""
        if self.formula_present:
            if not allow_cached_formula:
                raise ValueError("formula outcome requires explicit frozen cached-value permission")
            if self.xml_type not in {"n", ""}:
                raise ValueError("formula cache is not a numeric source cell")
        elif self.kind != "n" or self.xml_type not in {"n", ""}:
            raise ValueError("source cell is not an ordinary numeric measurement")
        if self.raw_v is None:
            raise ValueError("missing numeric source literal")
        try:
            decimal = Decimal(self.raw_v)
            if not decimal.is_finite():
                raise ValueError("nonfinite source decimal")
            value = Fraction(decimal)
        except (ValueError, ArithmeticError) as error:
            raise ValueError("invalid or nonfinite numeric source literal") from error
        try:
            finite_float = float(value)
        except OverflowError as error:
            raise ValueError("numeric source literal is outside finite float64") from error
        if not math.isfinite(finite_float):
            raise ValueError("numeric source literal is outside finite float64")
        return value

    def header_json(self) -> dict:
        """Only explicitly declared header cells are exposed in inventories."""
        if self.formula_present:
            value = None
        elif self.kind == "n":
            value = None
        elif self.kind != "s":
            value = None
        else:
            value = self.value
        return {"column_1based": self.column, "kind": self.kind, "value": value,
                "xml_type": self.xml_type, "formula_present": self.formula_present,
                "nontext_header_value_redacted": self.kind != "s" or self.formula_present}


class LiteralParser(WorkSheetParser):
    """Use a different XLSX library, retaining literals before its float casting."""
    def __init__(self, *args, selected_columns=None, metadata_only=False, header_rows=(), **kwargs):
        super().__init__(*args, **kwargs)
        self.selected_columns = selected_columns
        self.metadata_only = metadata_only
        self.header_rows = set(header_rows)

    def parse_cell(self, element):
        coordinate = element.get("r")
        row, column = coordinate_to_tuple(coordinate) if coordinate else (self.row_counter, self.col_counter + 1)
        if row != self.row_counter:
            raise ValueError("cell address disagrees with physically stored row")
        selected = self.selected_columns is None or column in self.selected_columns
        xml_type = element.attrib.get("t", "n")
        raw_v = element.findtext(VALUE_TAG)
        raw_formula = element.find(FORMULA_TAG)
        text_header = self.metadata_only and row in self.header_rows and xml_type in {"s", "str", "inlineStr"} and raw_formula is None
        if not selected or (self.metadata_only and not text_header):
            self.col_counter = column
            nonempty = raw_formula is not None or (raw_v is not None and raw_v != "")
            if xml_type == "inlineStr":
                nonempty = any((node.text or "") != "" for node in element.iter()
                               if node.tag.endswith("}t"))
            elif xml_type == "s" and raw_v is not None:
                nonempty = self.shared_strings[int(raw_v)] != ""
            kind = "f" if raw_formula is not None else ("s" if xml_type in {"str", "inlineStr"} else xml_type)
            return {"row": row, "column": column, "value": True if nonempty else None,
                    "data_type": kind, "source_raw_v": "" if raw_v is not None else None,
                    "source_xml_type": xml_type, "source_formula_present": raw_formula is not None,
                    "source_formula": None, "source_selected": selected}
        result = super().parse_cell(element)
        result["source_raw_v"] = raw_v
        result["source_xml_type"] = element.attrib.get("t", "n")
        result["source_formula_present"] = raw_formula is not None
        result["source_formula"] = None if raw_formula is None else raw_formula.text
        result["source_selected"] = True
        return result


class IndependentWorkbook:
    """Stream physical stored rows, including explicit blank/formula rows.

    Public iter_rows() fills worksheet gaps and can trust inflated dimensions.
    The pinned openpyxl parser instead yields actual stored <row> elements.
    This private-API dependency is explicit and version-checked, not inferred.
    """
    def __init__(self, path: Path):
        if openpyxl.__version__ != PINNED_OPENPYXL:
            raise ValueError("independent parser requires pinned openpyxl 3.1.5")
        with zipfile.ZipFile(path) as archive:
            entries = archive.infolist()
            if len({entry.filename for entry in entries}) != len(entries):
                raise ValueError("duplicate XLSX ZIP member names")
            if (any(entry.file_size > 512 * 1024 * 1024 for entry in entries)
                    or sum(entry.file_size for entry in entries) > 1024 * 1024 * 1024):
                raise ValueError("XLSX decompression exceeds fixed resource budget")
            if any(entry.flag_bits & 1 for entry in entries):
                raise ValueError("encrypted XLSX members are unsupported")
        self.workbook = openpyxl.load_workbook(path, read_only=True, data_only=False,
                                              keep_links=False, rich_text=False)

    @property
    def sheet_names(self) -> list[str]:
        return self.workbook.sheetnames

    def physical_rows(self, sheet: str, *, columns: set[int] | None = None,
                      metadata_only: bool = False, header_rows=()):
        if sheet not in self.sheet_names:
            raise ValueError(f"declared worksheet is absent: {sheet}")
        worksheet = self.workbook[sheet]
        previous = 0
        with worksheet._get_source() as source:
            parser = LiteralParser(source, worksheet._shared_strings, data_only=False,
                                   epoch=self.workbook.epoch, date_formats=self.workbook._date_formats,
                                   timedelta_formats=self.workbook._timedelta_formats, rich_text=False,
                                   selected_columns=columns, metadata_only=metadata_only, header_rows=header_rows)
            for row_index, parsed_cells in parser.parse():
                if row_index <= previous or row_index < 1:
                    raise ValueError("stored worksheet rows must increase strictly")
                previous = row_index
                cells = {}
                for parsed in parsed_cells:
                    if not parsed["source_selected"]:
                        continue
                    if parsed["row"] != row_index or parsed["column"] in cells:
                        raise ValueError("cell address disagrees with row or repeats a column")
                    cell = SourceCell(parsed["row"], parsed["column"], parsed["data_type"], parsed["value"],
                                      parsed["source_xml_type"], parsed["source_raw_v"],
                                      parsed["source_formula_present"], parsed["source_formula"])
                    cells[cell.column] = cell
                yield row_index, cells

    def inventory(self, declared_headers: dict[str, list[int]]) -> dict:
        if set(declared_headers) - set(self.sheet_names):
            raise ValueError("inventory declares an absent worksheet")
        result = []
        for name in self.sheet_names:
            worksheet = self.workbook[name]
            requested = set(declared_headers.get(name, []))
            if any(isinstance(row, bool) or not isinstance(row, int) or row < 1 for row in requested):
                raise ValueError("header rows must be positive integers")
            counters = Counter()
            kinds = Counter()
            header_cells = {}
            last_row, maximum_column = 0, 0
            for row_index, cells in self.physical_rows(name, metadata_only=True, header_rows=requested):
                counters["physical_stored_rows"] += 1
                counters["unrepresented_gaps_within_stored_extent"] += row_index - last_row - 1
                last_row = row_index
                counters["physical_stored_cells"] += len(cells)
                counters["meaningful_stored_rows"] += any(cell.meaningful for cell in cells.values())
                counters["explicit_stored_blank_rows"] += not any(cell.meaningful for cell in cells.values())
                if cells:
                    maximum_column = max(maximum_column, max(cells))
                for cell in cells.values():
                    kinds[cell.kind] += 1
                    counters["meaningful_stored_cells"] += cell.meaningful
                    counters["formula_cells"] += cell.formula_present
                    counters["formula_cells_with_cache_literal"] += cell.formula_present and cell.raw_v is not None
                    counters["blank_stored_cells"] += not cell.meaningful
                if row_index in requested:
                    header_cells[str(row_index)] = [cells[column].header_json() for column in sorted(cells)]
            if requested != {int(row) for row in header_cells}:
                raise ValueError("a declared header row is not physically stored")
            result.append({"sheet": name, "state": worksheet.sheet_state,
                           "declared_dimension_advisory": {"max_row": worksheet.max_row, "max_column": worksheet.max_column},
                           "largest_stored_row": last_row, "largest_stored_column": maximum_column,
                           "counts": dict(sorted(counters.items())), "cell_types": dict(sorted(kinds.items())),
                           "declared_headers_only": header_cells})
        return {"sheets": result, "numeric_outcome_summaries": "not computed by inventory",
                "row_accounting": "Physical stored rows; gaps and explicit blank rows are distinct; dimensions are advisory"}

    def close(self):
        self.workbook.close()

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        self.close()


def runtime() -> dict:
    return {"python": platform.python_version(), "openpyxl": openpyxl.__version__,
            "et_xmlfile": importlib.metadata.version("et-xmlfile"), "numpy": np.__version__,
            "recovery_protocol_sha256": FROZEN_RECOVERY_SHA256,
            "implementation_sha256": hashes(Path(__file__))["sha256"],
            "quantitative_source_downloads": "none; local checksum-bound inputs only"}


def extract_bound_table(workbook: IndependentWorkbook, *, sheet: str, header_row: int,
                         first_data_row: int, id_column: int, value_columns: dict[str, int],
                         expected_headers: dict[int, object], allow_cached_formula: bool = False) -> dict:
    """Extract only a future protocol's exact columns, retaining every stored row.

    No alias search, ID suffix transformation, missing-value imputation or
    outcome-dependent filtering occurs. Incomplete rows remain in the ledger.
    Header values and one-based column indices are explicit source contracts.
    """
    selected = {id_column, *value_columns.values()}
    if (first_data_row <= header_row or header_row < 1 or not value_columns
            or len(set(value_columns.values())) != len(value_columns)
            or id_column in value_columns.values() or not selected <= set(expected_headers)):
        raise ValueError("invalid explicit endpoint column/header binding")
    if any(isinstance(column, bool) or not isinstance(column, int) or column < 1 for column in selected):
        raise ValueError("bound column indices must be positive integers")
    records, excluded, seen = [], [], set()
    counts = Counter()
    missing_by_context = Counter()
    header_seen = False
    for number, cells in workbook.physical_rows(sheet, columns=selected | set(expected_headers)):
        if number == header_row:
            header_seen = True
            for column, expected in expected_headers.items():
                observed = None if column not in cells else cells[column].value
                if observed != expected:
                    raise ValueError(f"exact header mismatch at {sheet}/{number}/{column}")
        if number < first_data_row:
            continue
        counts["physical_data_rows"] += 1
        id_cell = cells.get(id_column)
        identifier = None if id_cell is None else id_cell.value
        if identifier is None or identifier == "" or (isinstance(identifier, str) and identifier.isspace()):
            counts["blank_id_rows"] += 1
            excluded.append({"source_row_1based": number, "reason": "blank_id"})
            continue
        if not isinstance(identifier, str) or id_cell.formula_present or "\n" in identifier or "\r" in identifier:
            raise ValueError("endpoint IDs must be explicit, exact, newline-free source strings")
        if identifier in seen:
            raise ValueError(f"duplicate endpoint ID: {identifier}")
        seen.add(identifier)
        record = {"id": identifier, "source_row_1based": number, "values": {}, "raw_literals": {},
                  "missing": {}, "formula_cache_used": []}
        for context, column in value_columns.items():
            cell = cells.get(column)
            try:
                if cell is None:
                    raise ValueError("absent_source_cell")
                exact = cell.numeric_fraction(allow_cached_formula=allow_cached_formula)
                record["values"][context] = float(exact)
                record["raw_literals"][context] = cell.raw_v
                if cell.formula_present:
                    record["formula_cache_used"].append(context)
            except ValueError as error:
                record["missing"][context] = {"reason": str(error),
                                                "kind": None if cell is None else cell.kind}
                missing_by_context[context] += 1
        record["complete"] = not record["missing"]
        counts["complete_ids"] += record["complete"]
        counts["incomplete_ids"] += not record["complete"]
        records.append(record)
    if not header_seen:
        raise ValueError("declared endpoint header row is absent")
    return {"sheet": sheet, "contexts": list(value_columns), "records": records,
            "accounting": {**dict(counts), "unique_ids": len(seen), "missing_by_context": dict(missing_by_context)},
            "excluded_stored_rows": excluded,
            "binding": {"header_row": header_row, "first_data_row": first_data_row, "id_column": id_column,
                        "value_columns": value_columns, "expected_headers": expected_headers,
                        "allow_cached_formula": allow_cached_formula}}


def exact_id_join(left: list[dict], right: list[dict], mapping: dict[str, str] | None = None) -> dict:
    """Exact joins or an explicitly supplied one-to-one mapping, never heuristics."""
    a, b = {row["id"]: row for row in left}, {row["id"]: row for row in right}
    if len(a) != len(left) or len(b) != len(right):
        raise ValueError("duplicate IDs cannot be silently joined")
    if mapping is None:
        pairs = [(identifier, identifier) for identifier in sorted(set(a) & set(b))]
        mapping_kind = "exact_case_sensitive_identity"
    else:
        if (not set(mapping) <= set(a) or not set(mapping.values()) <= set(b)
                or len(set(mapping.values())) != len(mapping)):
            raise ValueError("explicit ID mapping must name observed IDs one-to-one")
        pairs = sorted(mapping.items())
        mapping_kind = "explicit_external_mapping_requires_frozen_provenance"
    used_left, used_right = {x for x, _ in pairs}, {y for _, y in pairs}
    return {"pairs": pairs, "left_only_ids": sorted(set(a) - used_left), "right_only_ids": sorted(set(b) - used_right),
            "left_ids": len(a), "right_ids": len(b), "joined_ids": len(pairs), "mapping_kind": mapping_kind}


def compare_bound_tables(left: list[dict], right: list[dict], contexts: list[str], *,
                          absolute_tolerance: float, relative_tolerance: float = 0.0,
                          mapping: dict[str, str] | None = None) -> dict:
    """Compare every joined value, including exact source-literal arithmetic."""
    joined = exact_id_join(left, right, mapping)
    a, b = {row["id"]: row for row in left}, {row["id"]: row for row in right}
    deltas, signed, bounds, equal_rationals, equal_representations = [], [], [], 0, 0
    worst = None
    if any(not math.isfinite(value) or value < 0 for value in (absolute_tolerance, relative_tolerance)):
        raise ValueError("comparison tolerance must be explicit, finite and nonnegative")
    for left_id, right_id in joined["pairs"]:
        for context in contexts:
            if context not in a[left_id]["raw_literals"] or context not in b[right_id]["raw_literals"]:
                raise ValueError("joined endpoint values are incomplete; retain missingness before comparison")
            literal_a, literal_b = a[left_id]["raw_literals"][context], b[right_id]["raw_literals"][context]
            value_a, value_b = Fraction(literal_a), Fraction(literal_b)
            delta = float(abs(value_a - value_b))
            deltas.append(delta)
            signed.append(float(value_b - value_a))
            bounds.append(absolute_tolerance + relative_tolerance * abs(float(value_a)))
            equal_rationals += value_a == value_b
            equal_representations += literal_a == literal_b
            if worst is None or delta > worst["absolute_difference"]:
                worst = {"left_id": left_id, "right_id": right_id, "context": context,
                         "left_literal": literal_a, "right_literal": literal_b, "absolute_difference": delta}
    if not deltas:
        raise ValueError("endpoint comparison has no joined numeric cells")
    differences = np.array(deltas)
    return {"join": joined, "contexts": contexts, "values_compared": len(deltas),
            "equal_exact_numeric_literals": equal_rationals, "equal_source_representations": equal_representations,
            "absolute_tolerance": absolute_tolerance, "relative_tolerance": relative_tolerance,
            "reference_side": "left unchanged source", "exceeding_tolerance": int(np.count_nonzero(differences > np.array(bounds))),
            "maximum_absolute_difference": float(differences.max()), "mean_absolute_difference": float(differences.mean()),
            "signed_mean_error_candidate_minus_source": float(np.mean(signed)),
            "root_mean_squared_difference": float(np.sqrt(np.mean(differences ** 2))),
            "absolute_error_type7_quantiles_0_50_95_100": np.quantile(differences, [0, .5, .95, 1], method="linear").tolist(),
            "worst_source_cell": worst,
            "policy": "All joined source cells retained; no scaling, repair, endpoint selection or old-target mutation"}


def apply_fixed_equation(values, operations: tuple[str, ...], *, ddof: int,
                          population_rows: list[bool] | np.ndarray) -> tuple[np.ndarray, dict]:
    """Evaluate explicit fixed normalization operations; never infer a formula.

    Population membership, operation order and ddof must be supplied by the
    future frozen comparison. Complete finite rows are required. Unsupported
    or undefined operations fail instead of using nanmeans or zero scales.
    """
    matrix = np.asarray(values, dtype=np.float64)
    population = np.asarray(population_rows)
    if matrix.ndim != 2 or not len(matrix) or matrix.shape[1] < 2 or not np.isfinite(matrix).all():
        raise ValueError("equation input must be a complete finite element-by-context matrix")
    if population.dtype != bool or population.shape != (len(matrix),) or not population.any():
        raise ValueError("source normalization population must be an explicit nonempty boolean row mask")
    if isinstance(ddof, bool) or ddof not in {0, 1}:
        raise ValueError("ddof must be explicitly 0 or 1")
    if not operations or any(operation not in {"cell_zscore", "row_center"} for operation in operations):
        raise ValueError("only explicitly listed fixed cell_zscore/row_center operations are supported")
    result = matrix.copy()
    ledger = []
    for operation in operations:
        if operation == "cell_zscore":
            if population.sum() <= ddof:
                raise ValueError("normalization population is too small for specified ddof")
            mean = result[population].mean(axis=0)
            scale = result[population].std(axis=0, ddof=ddof)
            if not np.isfinite(scale).all() or np.any(scale <= 0):
                raise ValueError("constant or undefined source zscore column")
            result = (result - mean) / scale
            ledger.append({"operation": operation, "ddof": ddof, "source_population_rows": int(population.sum()),
                           "column_means": mean.tolist(), "column_scales": scale.tolist()})
        else:
            result -= result.mean(axis=1, keepdims=True)
            ledger.append({"operation": operation, "contexts_per_row": matrix.shape[1]})
    return result, {"operations": list(operations), "population_mask_sha256": hashlib.sha256(population.tobytes()).hexdigest(),
                    "ledger": ledger, "target_policy": "Diagnostic candidate only; no supplied or historical target is modified"}


def compare_fixed_matrix(observed, candidate, *, absolute_tolerance: float, relative_tolerance: float = 0.0) -> dict:
    actual = np.asarray(observed, dtype=np.float64)
    predicted = np.asarray(candidate, dtype=np.float64)
    if (actual.ndim != 2 or actual.shape != predicted.shape or not actual.size
            or not np.isfinite(actual).all() or not np.isfinite(predicted).all()
            or any(not math.isfinite(value) or value < 0 for value in (absolute_tolerance, relative_tolerance))):
        raise ValueError("fixed equation comparison needs aligned finite matrices and explicit nonnegative tolerance")
    residual = np.abs(actual - predicted)
    bounds = absolute_tolerance + relative_tolerance * np.abs(actual)
    by_context = []
    for context in range(actual.shape[1]):
        differences = residual[:, context]
        by_context.append({"context_index_0based": context, "values_compared": len(actual),
                           "exceeding_tolerance": int(np.count_nonzero(differences > bounds[:, context])),
                           "maximum_absolute_difference": float(differences.max()),
                           "mean_absolute_difference": float(differences.mean()),
                           "root_mean_squared_difference": float(np.sqrt(np.mean(differences ** 2))),
                           "signed_mean_error_candidate_minus_source": float(np.mean(predicted[:, context] - actual[:, context])),
                           "absolute_error_type7_quantiles_0_50_95_100": np.quantile(differences, [0, .5, .95, 1], method="linear").tolist()})
    return {"rows": len(actual), "contexts": actual.shape[1], "values_compared": int(actual.size),
            "absolute_tolerance": absolute_tolerance, "relative_tolerance": relative_tolerance,
            "exceeding_tolerance": int(np.count_nonzero(residual > bounds)),
            "maximum_absolute_difference": float(residual.max()), "mean_absolute_difference": float(residual.mean()),
            "root_mean_squared_difference": float(np.sqrt(np.mean(residual ** 2))),
            "signed_mean_error_candidate_minus_source": float(np.mean(predicted - actual)),
            "absolute_error_type7_quantiles_0_50_95_100": np.quantile(residual, [0, .5, .95, 1], method="linear").tolist(),
            "all_within_tolerance": bool(np.all(residual <= bounds)), "by_context": by_context}


def matrix_moments(values, contexts: list[str]) -> dict:
    """Descriptive source moments with explicitly counted finite populations."""
    matrix = np.asarray(values, dtype=float)
    if matrix.ndim != 2 or matrix.shape[1] != len(contexts) or np.isinf(matrix).any():
        raise ValueError("moment matrix must align contexts and use explicit NaN missingness")
    by_context = []
    for index, context in enumerate(contexts):
        finite = matrix[np.isfinite(matrix[:, index]), index]
        by_context.append({"context": context, "finite_count": len(finite),
                           "missing_count": len(matrix) - len(finite),
                           "mean": float(finite.mean()) if len(finite) else None,
                           "population_sd": float(finite.std(ddof=0)) if len(finite) else None,
                           "sample_sd": float(finite.std(ddof=1)) if len(finite) > 1 else None})
    complete = np.isfinite(matrix).all(axis=1)
    sums = matrix[complete].sum(axis=1)
    return {"source_rows": len(matrix), "complete_context_rows": int(complete.sum()), "by_context": by_context,
            "complete_row_sum": {"count": len(sums), "mean": float(sums.mean()) if len(sums) else None,
                                 "mean_absolute": float(np.abs(sums).mean()) if len(sums) else None,
                                 "root_mean_squared": float(np.sqrt(np.mean(sums ** 2))) if len(sums) else None,
                                 "maximum_absolute": float(np.abs(sums).max()) if len(sums) else None},
            "missingness_policy": "Finite populations counted explicitly; no values imputed and no nanmean used"}


def column_zscores(values, population_mask, *, ddof: int) -> tuple[np.ndarray, dict]:
    """Prespecified source populations with missingness retained, never imputed."""
    matrix = np.asarray(values, dtype=float)
    population = np.asarray(population_mask)
    if matrix.ndim != 2 or not matrix.size or np.isinf(matrix).any():
        raise ValueError("source matrix must be two-dimensional and contain only finite values or explicit missing NaNs")
    if population.dtype != bool or population.shape != matrix.shape:
        raise ValueError("per-cell normalization population must be an explicit aligned boolean matrix")
    if np.any(population & ~np.isfinite(matrix)):
        raise ValueError("normalization population contains missing values")
    if isinstance(ddof, bool) or ddof not in {0, 1}:
        raise ValueError("source zscore ddof must be explicitly 0 or 1")
    result = np.full_like(matrix, np.nan)
    moments = []
    for context in range(matrix.shape[1]):
        source = matrix[population[:, context], context]
        if len(source) <= ddof:
            raise ValueError("insufficient normalization population")
        mean = float(source.mean())
        scale = float(source.std(ddof=ddof))
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError("constant or undefined normalization column")
        finite = np.isfinite(matrix[:, context])
        result[finite, context] = (matrix[finite, context] - mean) / scale
        moments.append({"context_index_0based": context, "population_count": len(source),
                        "mean": mean, "standard_deviation": scale, "ddof": ddof})
    return result, {"moments": moments, "population_mask_sha256": hashlib.sha256(population.tobytes()).hexdigest()}


def replicate_activity_base(replicates, *, subtract_replicate_medians: bool) -> tuple[np.ndarray, dict]:
    """Mean exactly three source replicate values; incomplete triples stay missing."""
    values = np.asarray(replicates, dtype=float)
    if values.ndim != 3 or values.shape[2] != 3 or np.isinf(values).any():
        raise ValueError("replicate source must be element by context by exactly three finite-or-missing measurements")
    adjusted = values.copy()
    populations = []
    for context in range(values.shape[1]):
        for replicate in range(3):
            finite = np.isfinite(values[:, context, replicate])
            median = float(np.median(values[finite, context, replicate])) if finite.any() else None
            populations.append({"context_index_0based": context, "replicate_index_0based": replicate,
                                "finite_source_count": int(finite.sum()), "median": median})
            if subtract_replicate_medians and median is not None:
                adjusted[finite, context, replicate] -= median
    complete = np.isfinite(adjusted).all(axis=2)
    result = np.full(values.shape[:2], np.nan)
    result[complete] = adjusted[complete].mean(axis=1)
    return result, {"subtract_replicate_medians": subtract_replicate_medians,
                    "population": "All finite source scores independently per context/replicate, before old intersections",
                    "replicate_populations": populations, "complete_three_replicate_cells": int(complete.sum()),
                    "incomplete_three_replicate_cells": int((~complete).sum()), "imputation": False}


def fixed_grid_for_base(values, ids: list[str], base_id: str) -> list[dict]:
    """All twelve prespecified variants for one activity base, including failures."""
    matrix = np.asarray(values, dtype=float)
    if matrix.ndim != 2 or matrix.shape != (len(ids), 3) or len(set(ids)) != len(ids):
        raise ValueError("activity base must have unique IDs and exactly three context columns")
    complete = np.isfinite(matrix).all(axis=1)
    variants = []
    for population_name in ("complete3contextsourceintersection", "allfinite_independently_percell"):
        mask = np.repeat(complete[:, None], 3, axis=1) if population_name.startswith("complete3") else np.isfinite(matrix)
        for pre_ddof in (0, 1):
            for post_ddof in (None, 0, 1):
                key = f"{base_id}|{population_name}|pre_ddof{pre_ddof}|post_{'none' if post_ddof is None else 'ddof'+str(post_ddof)}"
                variant = {"id": key, "base_id": base_id, "pre_population": population_name,
                           "pre_ddof": pre_ddof, "post_ddof": post_ddof, "source_id_count": len(ids),
                           "complete3context_ids": [identifier for identifier, keep in zip(ids, complete) if keep],
                           "complete_id_sha256": hashlib.sha256("\n".join(identifier for identifier, keep in zip(ids, complete) if keep).encode()).hexdigest()}
                try:
                    if not complete.any():
                        raise ValueError("no complete three-context rows for this base")
                    standardized, pre_ledger = column_zscores(matrix, mask, ddof=pre_ddof)
                    centered = standardized[complete]
                    centered = centered - centered.mean(axis=1, keepdims=True)
                    variant["pre_moments"] = pre_ledger
                    if post_ddof is not None:
                        centered, post_ledger = column_zscores(centered, np.ones_like(centered, dtype=bool), ddof=post_ddof)
                        variant["post_moments"] = post_ledger
                    variant.update({"status": "estimable", "values": centered})
                except ValueError as error:
                    variant.update({"status": "nonestimable", "reason": str(error)})
                variants.append(variant)
    return variants


def published_ess_subset_test(source_records: list[dict], target_records: list[dict], contexts: list[str], *,
                               old_ids: set[str], old_admitted_ids: set[str]) -> dict:
    """The one frozen S11→all-S12-IDs sample-zscore test, without subset search."""
    source, target = {row["id"]: row for row in source_records}, {row["id"]: row for row in target_records}
    required = sorted(target)
    result = {"id": "pre_value_code_S11_ESS_to_all_exact_S12_ID_sample_zscore",
              "required_target_ids": len(target_records), "source_ids": len(source_records),
              "absolute_tolerance": 1e-12, "relative_tolerance": 0.0, "ddof": 1,
              "population_id_sha256": hashlib.sha256("\n".join(required).encode()).hexdigest(),
              "population_policy": "Entire exact S12 ID set; no missing ID/measurement may be dropped"}
    if len(source) != len(source_records) or len(target) != len(target_records):
        return dict(result, status="nonestimable", reason="duplicate source or target observation IDs")
    missing_ids = sorted(set(target) - set(source))
    if missing_ids:
        return dict(result, status="nonestimable", reason="S12 IDs lack exact S11 matches", missing_ids=missing_ids)
    incomplete = []
    for identifier in required:
        for side, rows in (("source", source), ("target", target)):
            values = rows[identifier].get("values", {})
            if any(context not in values or not math.isfinite(values[context]) for context in contexts):
                incomplete.append({"id": identifier, "side": side})
    if incomplete:
        return dict(result, status="nonestimable", reason="A required S12/source ID lacks a finite three-context measurement",
                    incomplete=incomplete)
    if len(required) < 2:
        return dict(result, status="nonestimable", reason="sample SD requires at least two exact S12 IDs")
    source_matrix = np.array([[source[identifier]["values"][context] for context in contexts] for identifier in required])
    observed = np.array([[target[identifier]["values"][context] for context in contexts] for identifier in required])
    try:
        candidate, moments = column_zscores(source_matrix, np.ones_like(source_matrix, dtype=bool), ddof=1)
    except ValueError as error:
        return dict(result, status="nonestimable", reason=str(error))
    comparisons = []
    cohorts = [("fullS12sourceIDset", set(required)), ("exactoldS9intersection", set(required) & old_ids),
               ("fixedoldadmittedintersection", set(required) & old_admitted_ids)]
    for name, included in cohorts:
        mask = np.array([identifier in included for identifier in required], dtype=bool)
        if not mask.any():
            comparisons.append({"cohort": name, "status": "nonestimable", "reason": "empty prescribed intersection"})
            continue
        comparison = compare_fixed_matrix(observed[mask], candidate[mask], absolute_tolerance=1e-12)
        comparison.update({"cohort": name, "status": "estimable", "id_count": int(mask.sum()),
                           "id_sha256": hashlib.sha256("\n".join(identifier for identifier in required if identifier in included).encode()).hexdigest(),
                           "normalization_recomputed_on_cohort": False})
        comparisons.append(comparison)
    return dict(result, status="estimable", normalization=moments, comparisons=comparisons,
                interpretation="Compatibility of two supplied releases only; no original value or historical pipeline is certified or changed")


def _write(path: Path, value) -> None:
    if path.suffix == ".gz":
        with path.open("xb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped:
            import io
            with io.TextIOWrapper(zipped, encoding="utf-8") as text:
                json.dump(value, text, allow_nan=False, separators=(",", ":"))
    else:
        with path.open("x") as text:
            json.dump(value, text, indent=2, allow_nan=False)
            text.write("\n")


def _json(path: Path):
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes())


def _parse_all_bound_sheets(binding: dict) -> dict:
    tables = {}
    with IndependentWorkbook(Path(binding["source_workbook"]["path"])) as workbook:
        for sheet, contract in binding["tables"].items():
            print("Independent openpyxl parsing: " + sheet, flush=True)
            headers = {column_index_from_string(column): label for column, label in contract["headers"].items()}
            header_stream = workbook.physical_rows(sheet, columns=set(headers))
            try:
                header_number, header_cells = next(header_stream)
            finally:
                header_stream.close()
            if header_number != contract["header_row"] or any(header_cells.get(column) is None or
                    header_cells[column].value != label for column, label in headers.items()):
                raise ValueError("source headers differ from frozen binding: " + sheet)
            id_column = column_index_from_string(contract["id_column"])
            selected = {column_index_from_string(column) for column in contract["value_columns"]}
            key_columns = [column_index_from_string(column) for column in contract["empirical_key_columns"]]
            keys, rows, missing, counts = set(), [], [], Counter()
            for number, cells in workbook.physical_rows(sheet, columns=selected | {id_column}):
                if number <= contract["header_row"]:
                    continue
                counts["stored_data_rows"] += 1
                id_cell = cells.get(id_column)
                identifier = None if id_cell is None else id_cell.value
                if not isinstance(identifier, str) or not identifier or id_cell.formula_present:
                    raise ValueError("invalid or missing bound source ID: " + sheet)
                values, numeric = {}, {}
                for column in sorted(selected):
                    label, cell = headers[column], cells.get(column)
                    raw = None if cell is None else (cell.raw_v if cell.kind == "n" else cell.value)
                    if cell is not None and cell.formula_present:
                        raise ValueError("unpermitted formula measurement: " + sheet)
                    values[label] = raw
                    try:
                        if cell is None:
                            raise ValueError("absent_source_cell")
                        cell.numeric_fraction()
                        numeric[label] = True
                    except ValueError as error:
                        numeric[label] = False
                        missing.append({"source_row": number, "id": identifier, "header": label,
                                        "xml_type": None if cell is None else cell.xml_type,
                                        "raw_token": raw, "reason": str(error)})
                empirical_key = tuple(identifier if column == id_column else values[headers[column]] for column in key_columns)
                if empirical_key in keys:
                    raise ValueError("duplicate empirical source key in " + sheet)
                keys.add(empirical_key)
                rows.append({"id": identifier, "source_row": number, "values": values, "finite_numeric": numeric})
            tables[sheet] = {"rows": rows, "counts": dict(counts), "missing_or_nonnumeric": missing,
                              "selected_headers": [headers[column] for column in sorted(selected)],
                              "unbound_columns_never_numeric_decoded": contract["unbound_first_row_columns_numeric_values_not_inspected"]}
    return tables


def _wide_records(table: dict, columns: dict[str, str]) -> list[dict]:
    result = []
    for row in table["rows"]:
        values, literals = {}, {}
        for context, header in columns.items():
            if row["finite_numeric"][header]:
                values[context] = float(Decimal(row["values"][header]))
                literals[context] = row["values"][header]
        result.append({"id": row["id"], "source_row_1based": row["source_row"], "values": values,
                       "raw_literals": literals, "complete": len(values) == len(columns)})
    if len({row["id"] for row in result}) != len(result):
        raise ValueError("duplicate wide endpoint IDs")
    return result


def _matrix_from_records(records: list[dict], contexts: list[str]) -> tuple[list[str], np.ndarray]:
    rows = sorted(records, key=lambda row: row["id"])
    return [row["id"] for row in rows], np.array([[row["values"].get(context, np.nan) for context in contexts] for row in rows])


def _cohort_comparisons(candidate_ids, candidate_values, target: list[dict], contexts, old_ids, admitted_ids) -> list[dict]:
    target_map = {row["id"]: row for row in target}
    lookup = {identifier: index for index, identifier in enumerate(candidate_ids)}
    complete_target = {identifier for identifier, row in target_map.items() if all(context in row["values"] for context in contexts)}
    common = set(candidate_ids) & complete_target
    result = []
    for name, cohort in (("full_appropriate_source_complete", common), ("exact_oldS9_intersection", common & old_ids),
                         ("fixed_old_admitted_intersection", common & admitted_ids)):
        ids = sorted(cohort)
        record = {"cohort": name, "matched_ids": len(ids), "id_sha256": hashlib.sha256("\n".join(ids).encode()).hexdigest(),
                  "candidate_ids_without_complete_target": len(set(candidate_ids) - complete_target),
                  "target_ids_without_candidate": len(complete_target - set(candidate_ids)), "normalization_refitted": False}
        if not ids:
            record.update(status="nonestimable", reason="empty prescribed source intersection")
        else:
            actual = np.array([[target_map[identifier]["values"][context] for context in contexts] for identifier in ids])
            predicted = np.array([candidate_values[lookup[identifier]] for identifier in ids])
            record.update(compare_fixed_matrix(actual, predicted, absolute_tolerance=1e-8, relative_tolerance=1e-8))
            record["status"] = "estimable"
        result.append(record)
    return result


def empirical_recovery(protocol_path: Path, binding_path: Path, output: Path, primary_export: Path | None) -> None:
    """Run only the exact centrally frozen source recovery; never a predictive fit."""
    if output.exists():
        raise FileExistsError("independent recovery output already exists; retain attempts")
    protocol_bytes, binding_bytes = protocol_path.read_bytes(), binding_path.read_bytes()
    if hashlib.sha256(protocol_bytes).hexdigest() != FROZEN_RECOVERY_SHA256:
        raise ValueError("recovery protocol differs from centrally frozen bytes")
    if hashlib.sha256(binding_bytes).hexdigest() != "a59592816f57f7254b780c83ba3e929e3b38f4fd06b2b196a0cddda1d9bfb98e":
        raise ValueError("exact pre-value schema binding differs")
    protocol, binding = json.loads(protocol_bytes), json.loads(binding_bytes)
    if binding["semantic_binding_complete"] is not True or binding["status"] != "bound_before_quantitative_inspection":
        raise ValueError("source schema gate is not complete")
    source = hashes(Path(binding["source_workbook"]["path"]))
    if source["sha256"] != binding["source_workbook"]["sha256"] or source["bytes"] != binding["source_workbook"]["bytes"]:
        raise ValueError("bound source workbook bytes changed")
    for evidence in binding["evidence"].values():
        if hashes(Path(evidence["path"]))["sha256"] != evidence["sha256"]:
            raise ValueError("schema provenance receipt changed")
    output.mkdir(parents=True)
    shutil.copyfile(__file__, output / "implementation.py")
    (output / "protocol.json").write_bytes(protocol_bytes)
    (output / "binding.json").write_bytes(binding_bytes)
    try:
        _write(output / "inputs.json", {"protocol_sha256": FROZEN_RECOVERY_SHA256,
               "schema_binding_sha256": hashlib.sha256(binding_bytes).hexdigest(), "source": source,
               "runtime": runtime(), "quantitative_open_started_at": datetime.datetime.now(datetime.timezone.utc).isoformat()})
        tables = _parse_all_bound_sheets(binding)
        _write(output / "independent-source-values.json.gz", tables)
        if primary_export is not None:
            expected = _json(primary_export)
            if (expected["protocol_sha256"] != FROZEN_RECOVERY_SHA256
                    or expected["schema_binding_sha256"] != hashlib.sha256(binding_bytes).hexdigest()
                    or expected["source"]["sha256"] != source["sha256"]):
                raise ValueError("primary source export binds a different protocol, schema or source")
            discrepancies = []
            for name, table in tables.items():
                projected = [{key: row[key] for key in ("id", "source_row", "values")} for row in table["rows"]]
                if projected != expected["tables"][name]["rows"]:
                    discrepancies.append(name)
            _write(output / "source-parser-comparison.json", {"all_source_rows_and_lexemes_exact": not discrepancies,
                   "differing_tables": discrepancies, "primary_export_sha256": hashes(primary_export)["sha256"],
                   "tables": len(tables), "stored_data_rows": sum(len(table["rows"]) for table in tables.values())})
            if discrepancies:
                raise ValueError("independent raw workbook parser differs from primary source export")
        contexts = protocol["endpoint_contract"]["contexts"]
        aggregate = _wide_records(tables[binding["roles"]["activity_aggregate"]["sheet"]], binding["roles"]["activity_aggregate"]["columns"])
        ess = _wide_records(tables[binding["roles"]["published_ess"]["sheet"]], binding["roles"]["published_ess"]["columns"])
        old_contract = protocol["old_inputs"]["s9_workbook"]
        if hashes(Path(old_contract["path"]))["sha256"] != old_contract["sha256"]:
            raise ValueError("old S9 source changed")
        with IndependentWorkbook(Path(old_contract["path"])) as old_book:
            old_table = extract_bound_table(old_book, sheet=old_contract["sheet"], header_row=1, first_data_row=2,
                id_column=2, value_columns={context: index + 3 for index, context in enumerate(contexts)},
                expected_headers={2: old_contract["id_column"], **{index + 3: old_contract["context_columns"][context]
                                  for index, context in enumerate(contexts)}})
        old = old_table["records"]
        old_ids = {row["id"] for row in old}
        if len(old_ids) != old_contract["nonempty_unique_ids"]:
            raise ValueError("independent old S9 parser lost a required source ID")
        admitted_contract = protocol["old_inputs"]["old_admitted_panel"]
        if hashes(Path(admitted_contract["path"]))["sha256"] != admitted_contract["sha256"]:
            raise ValueError("fixed old admitted panel changed")
        admitted_ids = {row["id"] for row in _json(Path(admitted_contract["path"]))["records"]}
        if len(admitted_ids) != admitted_contract["element_count"]:
            raise ValueError("fixed old admitted ID count differs")
        _write(output / "old-source-values.json.gz", old_table)
        ids, aggregate_values = _matrix_from_records(aggregate, contexts)
        ess_ids, ess_values = _matrix_from_records(ess, contexts)
        all_rep_ids = sorted({row["id"] for context in contexts for row in tables[binding["roles"]["replicates"][context]["sheet"]]["rows"]})
        rep_lookup = {identifier: index for index, identifier in enumerate(all_rep_ids)}
        replicate_cube = np.full((len(all_rep_ids), 3, 3), np.nan)
        for c, context in enumerate(contexts):
            role = binding["roles"]["replicates"][context]
            for row in tables[role["sheet"]]["rows"]:
                token = row["values"][role["replicate_header"]]
                if token not in role["replicate_values"].values():
                    raise ValueError("source replicate token differs from exact binding")
                replicate = int(token) - 1
                if row["finite_numeric"][role["activity_header"]]:
                    replicate_cube[rep_lookup[row["id"]], c, replicate] = float(Decimal(row["values"][role["activity_header"]]))
        mean_replicates, mean_ledger = replicate_activity_base(replicate_cube, subtract_replicate_medians=False)
        centered_replicates, centered_ledger = replicate_activity_base(replicate_cube, subtract_replicate_medians=True)
        bases = [("published_activity_aggregate", ids, aggregate_values, None),
                 ("mean_three_published_replicates", all_rep_ids, mean_replicates, mean_ledger),
                 ("mean_three_after_explicit_replicate_median_subtraction", all_rep_ids, centered_replicates, centered_ledger)]
        pipelines = []
        (output / "pipelines").mkdir()
        base_reports = []
        for base_id, base_ids, values, ledger in bases:
            finite = np.isfinite(values).all(axis=1)
            base_reports.append({"id": base_id, "moments": matrix_moments(values, contexts), "replicate_arithmetic": ledger,
                "against_published_aggregate": _cohort_comparisons([identifier for identifier, keep in zip(base_ids, finite) if keep],
                    values[finite], aggregate, contexts, old_ids, admitted_ids)})
            for variant in fixed_grid_for_base(values, base_ids, base_id):
                if variant["status"] == "estimable":
                    matrix = variant.pop("values")
                    name = f"pipeline-{len(pipelines):02d}.json.gz"
                    _write(output / "pipelines" / name, {"ids": variant["complete3context_ids"], "contexts": contexts, "values": matrix.tolist()})
                    variant["values_file"] = "pipelines/" + name
                    variant["moments"] = matrix_moments(matrix, contexts)
                    variant["against_final_published_ESS"] = _cohort_comparisons(variant["complete3context_ids"], matrix, ess, contexts, old_ids, admitted_ids)
                    variant["against_old_S9"] = _cohort_comparisons(variant["complete3context_ids"], matrix, old, contexts, old_ids, admitted_ids)
                pipelines.append(variant)
        code_test = published_ess_subset_test(ess, old, contexts, old_ids=old_ids, old_admitted_ids=admitted_ids)
        summary = {"schema_version": 1, "protocol_sha256": FROZEN_RECOVERY_SHA256, "schema_binding_sha256": hashlib.sha256(binding_bytes).hexdigest(),
                   "runtime": runtime(), "source_row_accounting": {name: {"rows": len(table["rows"]),
                       "missing_or_nonnumeric_selected_cells": len(table["missing_or_nonnumeric"])} for name, table in tables.items()},
                   "old_unique_ids": len(old_ids), "old_admitted_ids": len(admitted_ids),
                   "activity_bases": base_reports, "pipelines": pipelines,
                   "final_published_ESS_source_moments": matrix_moments(ess_values, contexts),
                   "published_ESS_vs_old_S9": compare_bound_tables(ess, old, contexts, absolute_tolerance=1e-8, relative_tolerance=1e-8),
                   "pre_value_code_hypothesis": code_test, "predictive_fitting_performed": False,
                   "historical_pipeline_identified": False, "source_values_or_old_targets_modified": False,
                   "interpretation": "All prescribed diagnostic routes retained. Compatibility is not producer history or independent biological validation."}
        _write(output / "summary.json", summary)
        print(json.dumps({"summary": str(output / "summary.json"), "variants": len(pipelines),
                          "estimable": sum(row["status"] == "estimable" for row in pipelines),
                          "code_test_status": code_test["status"]}), flush=True)
    except BaseException as error:
        _write(output / "failure.json", {"error_type": type(error).__name__, "error": str(error), "traceback": traceback.format_exc()})
        raise
    finally:
        files = {str(path.relative_to(output)): hashes(path)["sha256"] for path in sorted(output.rglob("*")) if path.is_file() and path.name != "manifest.json"}
        _write(output / "manifest.json", {"schema_version": 1, "protocol_sha256": FROZEN_RECOVERY_SHA256, "files": files})


def restored_source_arrays(root: Path) -> tuple[dict, list[dict], list[dict], set[str]]:
    """Replay saved source values, reusing the OLD archived replicate medians.

    No normalization population or median is estimated here. This also supplies
    immutable pre-existing activity bases for any separately authorized audit.
    """
    tables, binding, summary = _json(root / "independent-source-values.json.gz"), _json(root / "binding.json"), _json(root / "summary.json")
    contexts = ["HepG2", "K562", "WTC11"]
    aggregate = _wide_records(tables[binding["roles"]["activity_aggregate"]["sheet"]], binding["roles"]["activity_aggregate"]["columns"])
    ess = _wide_records(tables[binding["roles"]["published_ess"]["sheet"]], binding["roles"]["published_ess"]["columns"])
    ids, aggregate_values = _matrix_from_records(aggregate, contexts)
    rep_ids = sorted({row["id"] for context in contexts for row in tables[binding["roles"]["replicates"][context]["sheet"]]["rows"]})
    lookup = {identifier: index for index, identifier in enumerate(rep_ids)}
    cube = np.full((len(rep_ids), 3, 3), np.nan)
    for c, context in enumerate(contexts):
        role = binding["roles"]["replicates"][context]
        for row in tables[role["sheet"]]["rows"]:
            if row["finite_numeric"][role["activity_header"]]:
                cube[lookup[row["id"]], c, int(row["values"][role["replicate_header"]]) - 1] = float(Decimal(row["values"][role["activity_header"]]))
    complete = np.isfinite(cube).all(axis=2)
    mean = np.full(cube.shape[:2], np.nan)
    mean[complete] = cube[complete].mean(axis=1)
    centered = cube.copy()
    old_medians = summary["activity_bases"][2]["replicate_arithmetic"]["replicate_populations"]
    for item in old_medians:
        c, r, median = item["context_index_0based"], item["replicate_index_0based"], item["median"]
        if median is not None:
            finite = np.isfinite(centered[:, c, r])
            centered[finite, c, r] -= median
    centered_mean = np.full(cube.shape[:2], np.nan)
    centered_mean[complete] = centered[complete].mean(axis=1)
    bases = {"published_activity_aggregate": (ids, aggregate_values), "mean_three_published_replicates": (rep_ids, mean),
             "mean_three_after_explicit_replicate_median_subtraction": (rep_ids, centered_mean)}
    old = _json(root / "old-source-values.json.gz")["records"]
    protocol = _json(root / "protocol.json")
    admitted_path = Path(protocol["old_inputs"]["old_admitted_panel"]["path"])
    if hashes(admitted_path)["sha256"] != protocol["old_inputs"]["old_admitted_panel"]["sha256"]:
        raise ValueError("old admitted ID source changed")
    admitted = {row["id"] for row in _json(admitted_path)["records"]}
    return bases, ess, old, admitted


def verified_manifest(root: Path, filename: str) -> dict:
    directory = root.resolve()
    manifest = _json(directory / filename)
    for name, expected in manifest["files"].items():
        path = (directory / name).resolve()
        if directory not in path.parents or Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError("artifact path escapes recovery bundle")
        if hashes(path)["sha256"] != expected:
            raise ValueError("recovery artifact integrity mismatch: " + name)
    return {"path": str(directory / filename), "sha256": hashes(directory / filename)["sha256"]}


class RecoveryComparison:
    def __init__(self):
        self.checks, self.failures = [], []

    def exact(self, name, expected, observed):
        matched = expected == observed
        self.checks.append({"category": name, "kind": "exact", "matched": matched})
        if not matched:
            self.failures.append({"category": name, "kind": "exact"})

    def numeric(self, name, expected, observed, absolute=1e-8, relative=1e-8):
        a, b = np.asarray(expected, dtype=float), np.asarray(observed, dtype=float)
        if a.shape != b.shape:
            self.failures.append({"category": name, "kind": "shape", "expected": list(a.shape), "observed": list(b.shape)})
            return
        missing = np.isnan(a)
        finite = np.isfinite(a) & np.isfinite(b)
        invalid = ~(finite | (missing & np.isnan(b)))
        delta = np.abs(a[finite] - b[finite])
        failed = bool(invalid.any() or np.any(delta > absolute + relative * np.abs(a[finite])))
        item = {"category": name, "kind": "numeric", "matched": not failed, "finite_values": int(finite.sum()),
                "matching_missing": int((missing & np.isnan(b)).sum()), "maximum_absolute_difference": float(delta.max()) if len(delta) else 0.0,
                "absolute_tolerance": absolute, "relative_tolerance": relative}
        self.checks.append(item)
        if failed:
            self.failures.append(item)


def _json_id_hash(ids) -> str:
    return hashlib.sha256(json.dumps(sorted(ids), separators=(",", ":")).encode()).hexdigest()


def _compare_metric(checks: RecoveryComparison, prefix: str, own: dict, primary: dict):
    checks.exact(prefix + ":cells", own["values_compared"], primary["paired_finite_cells"])
    checks.exact(prefix + ":within_count", own["values_compared"] - own["exceeding_tolerance"], primary["within_fixed_tolerance_count"])
    for a, b in (("maximum_absolute_difference", "maximum_absolute_error"), ("mean_absolute_difference", "mae"),
                 ("root_mean_squared_difference", "rmse"), ("signed_mean_error_candidate_minus_source", "signed_mean_error")):
        checks.numeric(prefix + ":" + a, own[a], primary[b])
    checks.numeric(prefix + ":quantiles", own["absolute_error_type7_quantiles_0_50_95_100"],
                   [primary["absolute_error_type7_quantiles"][key] for key in ("0", "0.5", "0.95", "1")])


def _recovery_comparison_compute(independent_root: Path, primary_root: Path, output: Path) -> bool:
    own_manifest, primary_manifest = verified_manifest(independent_root, "manifest.json"), verified_manifest(primary_root, "artifacts.json")
    own, primary = _json(independent_root / "summary.json"), _json(primary_root / "summary.json")
    if own["protocol_sha256"] != FROZEN_RECOVERY_SHA256 or primary["protocol_sha256"] != FROZEN_RECOVERY_SHA256:
        raise ValueError("recovery protocol mismatch")
    bases, ess, old, admitted = restored_source_arrays(independent_root)
    contexts = ["HepG2", "K562", "WTC11"]
    old_ids = {row["id"] for row in old}
    checks = RecoveryComparison()
    pops = {"complete3contextsourceintersection": "complete_three_context_source", "allfinite_independently_percell": "all_finite_per_context"}
    def own_key(row):
        return row["base_id"], pops[row["pre_population"]], row["pre_ddof"], "none" if row["post_ddof"] is None else "column_zscore_ddof" + str(row["post_ddof"])
    def primary_key(row):
        return row["activity_base"], row["pre_zscore_population"], row["pre_zscore_ddof"], row["post_ESS_step"]
    a_map, b_map = {own_key(row): row for row in own["pipelines"]}, {primary_key(row): row for row in primary["variants"]}
    checks.exact("all36_parameter_tuples", sorted(a_map), sorted(b_map))
    for key in sorted(a_map):
        a, b = a_map[key], b_map[key]
        label = "|".join(map(str, key))
        checks.exact(label + ":estimability", a["status"] == "estimable", b["estimable"])
        if a["status"] != "estimable" or not b["estimable"]:
            continue
        a_values, b_values = _json(independent_root / a["values_file"]), _json(primary_root / b["derived_values_file"])
        checks.exact(label + ":all_source_IDs", a_values["ids"], b_values["ids"])
        checks.exact(label + ":complete_ID_commitment", _json_id_hash(a_values["ids"]), b["complete_source_population"]["ids_sha256"])
        checks.numeric(label + ":all_derived_cells", a_values["values"], b_values["values"])
        base_ids, base_values = bases[a["base_id"]]
        complete = np.isfinite(base_values).all(axis=1)
        for c, context in enumerate(contexts):
            params, expected = a["pre_moments"]["moments"][c], b["pre_parameters"][c]
            selected = complete if key[1] == "complete_three_context_source" else np.isfinite(base_values[:, c])
            checks.exact(label + ":preN:" + context, params["population_count"], expected["count"])
            checks.exact(label + ":pre_IDs:" + context, _json_id_hash([identifier for identifier, keep in zip(base_ids, selected) if keep]), expected["ids_sha256"])
            checks.numeric(label + ":pre_moments:" + context, [params["mean"], params["standard_deviation"]], [expected["mean"], expected["sd"]])
            if a["post_ddof"] is not None:
                params, expected = a["post_moments"]["moments"][c], b["post_parameters"][c]
                checks.exact(label + ":postN:" + context, params["population_count"], expected["count"])
                checks.exact(label + ":post_IDs:" + context, _json_id_hash(a_values["ids"]), expected["ids_sha256"])
                checks.numeric(label + ":post_moments:" + context, [params["mean"], params["standard_deviation"]], [expected["mean"], expected["sd"]])
        for own_name, primary_name in (("against_final_published_ESS", "final_S11_ESS"), ("against_old_S9", "unchanged_old_S9")):
            for a_cohort, name in zip(a[own_name], ("full_source_population", "exact_old_s9_intersection", "fixed_old_admitted_intersection")):
                b_cohort = b["comparisons"][primary_name]["cohorts"][name]
                prefix = label + ":" + primary_name + ":" + name
                checks.exact(prefix + ":paired_ids", a_cohort["matched_ids"], b_cohort["complete_three_context_paired_ids"])
                _compare_metric(checks, prefix + ":pooled", a_cohort, b_cohort["complete_three_context_pooled"])
                for c, context in enumerate(contexts):
                    _compare_metric(checks, prefix + ":" + context, a_cohort["by_context"][c], b_cohort["complete_three_context_metrics"][context])
    for base in own["activity_bases"]:
        if base["id"] == "published_activity_aggregate":
            continue
        for a, name in zip(base["against_published_aggregate"], ("full_source_population", "exact_old_s9_intersection", "fixed_old_admitted_intersection")):
            b = primary["aggregate_comparisons"][base["id"]]["cohorts"][name]
            checks.exact(base["id"] + ":aggregate:" + name + ":N", a["matched_ids"], b["complete_three_context_paired_ids"])
            _compare_metric(checks, base["id"] + ":aggregate:" + name, a, b["complete_three_context_pooled"])
    source_ids, source_matrix = _matrix_from_records(ess, contexts)
    source_lookup = dict(zip(source_ids, source_matrix))
    target_ids = sorted(old_ids)
    subset = np.array([source_lookup[identifier] for identifier in target_ids])
    code_values, _ = column_zscores(subset, np.ones_like(subset, dtype=bool), ddof=1)
    primary_code = _json(primary_root / "code-hypothesis-values.json.gz")
    checks.exact("code_test_entire_S12_IDs", target_ids, primary_code["ids"])
    checks.numeric("code_test_all_derived_cells", code_values, primary_code["values"], absolute=1e-12, relative=0.0)
    for a, name in zip(own["pre_value_code_hypothesis"]["comparisons"], ("full_s12_source_id_set", "exact_old_s9_intersection", "fixed_old_admitted_intersection")):
        b = primary["pre_value_code_hypothesis"]["fixed_1e-12_test_by_cohort"][name]
        checks.exact("code_test:" + name + ":N", a["id_count"], b["ids"])
        checks.exact("code_test:" + name + ":within", a["values_compared"] - a["exceeding_tolerance"], b["within_1e-12_cells"])
        checks.exact("code_test:" + name + ":pass", a["all_within_tolerance"], b["matched"])
    report = {"schema_version": 1, "protocol_sha256": FROZEN_RECOVERY_SHA256,
              "primary_artifacts_manifest": primary_manifest, "independent_manifest": own_manifest,
              "matched": not checks.failures, "failure_count": len(checks.failures), "failures": checks.failures, "checks": checks.checks,
              "ID_hash_conventions": {"primary": "SHA256 compact JSON of sorted IDs", "independent_original": "SHA256 newline-joined sorted IDs; exact lists compared and primary commitments recomputed"},
              "limits": "Numerical/source compatibility only; no old target mutation, predictive fit or producer-history certification"}
    _write(output / "comparison.json", report)
    print(json.dumps({"matched": report["matched"], "failure_count": report["failure_count"], "checks": len(report["checks"])}))
    return report["matched"]


def recovery_comparison(independent_root: Path, primary_root: Path, output: Path) -> bool:
    if output.exists():
        raise FileExistsError("comparison directory exists; retain attempts")
    output.mkdir(parents=True)
    shutil.copyfile(__file__, output / "implementation.py")
    try:
        return _recovery_comparison_compute(independent_root, primary_root, output)
    except BaseException as error:
        _write(output / "failure.json", {"error_type": type(error).__name__, "error": str(error), "traceback": traceback.format_exc()})
        raise
    finally:
        _write(output / "manifest.json", {"schema_version": 1, "protocol_sha256": FROZEN_RECOVERY_SHA256,
            "files": {path.name: hashes(path)["sha256"] for path in sorted(output.iterdir()) if path.is_file() and path.name != "manifest.json"}})


def inventory_file(path: Path, expected_sha256: str, headers: dict[str, list[int]], output: Path,
                    expected_md5: str | None = None) -> None:
    if output.exists():
        raise FileExistsError("inventory output already exists; retain attempts")
    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        raise ValueError("expected SHA256 must be explicit lowercase 64-digit hex")
    observed = hashes(path)
    if observed["sha256"] != expected_sha256 or (expected_md5 is not None and observed["md5"] != expected_md5):
        raise ValueError("local source workbook integrity differs from declared payload")
    with IndependentWorkbook(path) as workbook:
        inventory = workbook.inventory(headers)
    report = {"schema_version": 1, "status": "independent_workbook_inventory_only",
              "source": {"path": str(path), **observed}, "runtime": runtime(),
              "declared_header_rows": headers, **inventory,
              "limits": "No endpoint acceptance, source equation, target correction or model fit is implied"}
    output.mkdir(parents=True)
    shutil.copyfile(__file__, output / "implementation.py")
    (output / "inventory.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    manifest = {"schema_version": 1, "files": {name: hashes(output / name)["sha256"]
                for name in ("implementation.py", "inventory.json")}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"inventory": str(output / "inventory.json"), "worksheets": len(report["sheets"])}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xlsx", type=Path)
    parser.add_argument("--expected-sha256")
    parser.add_argument("--expected-md5")
    parser.add_argument("--header", action="append", default=[], metavar="SHEET:ROW1")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--describe", action="store_true")
    parser.add_argument("--empirical", action="store_true")
    parser.add_argument("--protocol", type=Path, default=Path("experiments/context-endpoint-recovery-v1.json"))
    parser.add_argument("--binding", type=Path, default=Path("research/context-endpoint-recovery/source/schema-binding.json"))
    parser.add_argument("--primary-export", type=Path)
    parser.add_argument("--compare-primary", type=Path)
    parser.add_argument("--independent-run", type=Path)
    args = parser.parse_args()
    if args.compare_primary is not None:
        if args.out is None or args.independent_run is None:
            parser.error("--compare-primary requires --independent-run and new --out")
        if not recovery_comparison(args.independent_run, args.compare_primary, args.out):
            raise SystemExit(2)
        return
    if args.empirical:
        if args.out is None:
            parser.error("--empirical requires a new --out directory")
        empirical_recovery(args.protocol, args.binding, args.out, args.primary_export)
        return
    if args.describe or args.xlsx is None:
        print(json.dumps(runtime(), indent=2))
        return
    if args.expected_sha256 is None or args.out is None:
        parser.error("inventory requires --xlsx, --expected-sha256 and new --out")
    headers = {}
    for specification in args.header:
        try:
            sheet, row = specification.rsplit(":", 1)
            integer = int(row)
        except (ValueError, AttributeError):
            parser.error("--header must be an exact SHEET:ROW1 pair")
        if integer in headers.setdefault(sheet, []):
            parser.error("duplicate declared header row")
        headers[sheet].append(integer)
    inventory_file(args.xlsx, args.expected_sha256, headers, args.out, args.expected_md5)


if __name__ == "__main__":
    main()
