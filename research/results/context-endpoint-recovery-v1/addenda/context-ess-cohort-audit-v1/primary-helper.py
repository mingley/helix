#!/usr/bin/env python3
"""Frozen 36-route source arithmetic audit; no fitting or endpoint replacement."""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation, localcontext
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import shutil

import numpy as np


PROTOCOL_SHA256 = "00175e41cf9fa0cbf177ce8fc7ef143eb508aec8b11d058f7acb4458284dc868"
CONTEXTS = ("HepG2", "K562", "WTC11")
BASES = ("published_activity_aggregate", "mean_three_published_replicates",
         "mean_three_after_explicit_replicate_median_subtraction")
POPULATIONS = ("complete_three_context_source", "all_finite_per_context")
POST = ("none", "column_zscore_ddof0", "column_zscore_ddof1")
ABS_TOL, REL_TOL = 1e-8, 1e-8


def sha(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def ids_hash(ids):
    return hashlib.sha256(json.dumps(sorted(ids), ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def write_gzip(path, value):
    with Path(path).open("xb") as raw:
        with gzip.GzipFile(fileobj=raw, filename="", mode="wb", mtime=0) as handle:
            handle.write(json.dumps(value, allow_nan=False, separators=(",", ":")).encode())


def source_number(raw):
    """Preserve tokens elsewhere; no text imputation, cached formula or float-to-decimal claims."""
    if raw is None:
        return math.nan, "absent_cell"
    if not isinstance(raw, str):
        raise ValueError("exported source values must be raw strings or null")
    try:
        value = Decimal(raw)
    except InvalidOperation:
        return math.nan, "missing_or_nonnumeric_token:" + raw
    if not value.is_finite():
        return math.nan, "nonfinite_decimal:" + raw
    result = float(value)
    return (result, None) if math.isfinite(result) else (math.nan, "outside_finite_float64")


@dataclass
class Endpoint:
    name: str
    ids: list[str]
    values: np.ndarray
    units: str
    raw: dict[str, list[str | None]] | None = None
    blocked: str | None = None
    ledger: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.values.shape != (len(self.ids), 3) or len(set(self.ids)) != len(self.ids):
            raise ValueError("endpoint matrix or identifier uniqueness differs")
        self.index = {identifier: i for i, identifier in enumerate(self.ids)}

    def subset(self, ids):
        result = np.full((len(ids), 3), np.nan)
        present = [(i, self.index[identifier]) for i, identifier in enumerate(ids) if identifier in self.index]
        if present:
            left, right = zip(*present)
            result[list(left)] = self.values[list(right)]
        return result


def table_role(tables, role, fields):
    if not role:
        return None
    table = tables[role["sheet"]]
    if table["id_header"] != role["id_header"]:
        raise ValueError("bound table identifier header differs")
    headers = set(table["headers"].values())
    if not set(fields) <= headers or not set(fields) <= set(table["selected_value_headers"]):
        raise ValueError("bound numeric field is absent from selected export headers")
    return table


def direct_endpoint(name, tables, role, units):
    if not role:
        return Endpoint(name, [], np.empty((0, 3)), units, blocked="required source role absent")
    columns = [role["columns"][context] for context in CONTEXTS]
    table = table_role(tables, role, columns)
    by_id, ledger, duplicates = {}, [], []
    for row in table["rows"]:
        identifier = row["id"]
        raw = [row["values"].get(header) for header in columns]
        parsed = [source_number(value) for value in raw]
        reason = [item[1] for item in parsed]
        entry = {"id": identifier, "source_row": row["source_row"], "missing_by_context": dict(zip(CONTEXTS, reason))}
        if not isinstance(identifier, str) or not identifier:
            entry["identifier_status"] = "missing_or_invalid"
        elif identifier in by_id:
            entry["identifier_status"] = "duplicate_blocks_endpoint"
            duplicates.append(identifier)
        else:
            entry["identifier_status"] = "unique"
            by_id[identifier] = (raw, [item[0] for item in parsed])
        ledger.append(entry)
    identifiers = sorted(by_id)
    values = np.asarray([by_id[i][1] for i in identifiers], dtype=np.float64).reshape((-1, 3))
    return Endpoint(name, identifiers, values, units, raw={i: by_id[i][0] for i in identifiers},
                    blocked="duplicate source IDs block this endpoint" if duplicates else None,
                    ledger={"sheet": role["sheet"], "rows": ledger, "duplicate_ids": sorted(set(duplicates)),
                            "source_columns": dict(zip(CONTEXTS, columns))})


def replicate_bases(tables, roles):
    """Long empirical keys are (exact ID, declared replicate), never ID alone."""
    observations, ledger, identifiers, blocked = {}, {}, set(), []
    for context in CONTEXTS:
        role = (roles or {}).get(context)
        if not role:
            blocked.append(context + ": required replicate role absent")
            continue
        wide = "replicate_columns" in role
        columns = ([role["replicate_columns"][str(k)] for k in (1, 2, 3)] if wide else
                   [role["replicate_header"], role["activity_header"]])
        table = table_role(tables, role, columns)
        entries, seen = [], set()
        for row in table["rows"]:
            identifier = row["id"]
            entry = {"id": identifier, "source_row": row["source_row"], "observations": []}
            if not isinstance(identifier, str) or not identifier:
                entry["identifier_status"] = "missing_or_invalid"
                entries.append(entry)
                continue
            identifiers.add(identifier)
            if wide:
                observed = [(k, row["values"].get(role["replicate_columns"][str(k)])) for k in (1, 2, 3)]
            else:
                token = row["values"].get(role["replicate_header"])
                mapping = role.get("replicate_values", {str(k): str(k) for k in (1, 2, 3)})
                matches = [k for k in (1, 2, 3) if token == mapping[str(k)]]
                if len(matches) != 1:
                    entry["replicate_status"] = "unknown_or_missing_source_replicate_token"
                    entry["source_replicate_token"] = token
                    blocked.append(context + ": undeclared replicate token")
                    entries.append(entry)
                    continue
                observed = [(matches[0], row["values"].get(role["activity_header"]))]
            for replicate, raw in observed:
                key = (identifier, context, replicate)
                value, reason = source_number(raw)
                duplicate = key in seen
                if duplicate:
                    blocked.append(context + ": duplicate (ID,replicate) empirical key")
                else:
                    seen.add(key)
                    observations[key] = value
                entry["observations"].append({"replicate": replicate, "missing_reason": reason,
                                               "duplicate_empirical_key": duplicate})
            entries.append(entry)
        ledger[context] = {"sheet": role["sheet"], "format": "wide" if wide else "long", "rows": entries}
    identifiers = sorted(identifiers)
    cube = np.asarray([[[observations.get((i, c, r), math.nan) for r in (1, 2, 3)]
                        for c in CONTEXTS] for i in identifiers], dtype=np.float64).reshape((-1, 3, 3))
    reason = "; ".join(sorted(set(blocked))) if blocked else None
    medians, populations = np.full((3, 3), np.nan), []
    if not reason:
        for c, context in enumerate(CONTEXTS):
            for r in range(3):
                valid = np.isfinite(cube[:, c, r])
                if np.any(valid):
                    medians[c, r] = np.median(cube[valid, c, r])
                populations.append({"context": context, "replicate": r + 1, "count": int(valid.sum()),
                                    "ids_sha256": ids_hash(np.asarray(identifiers)[valid].tolist()),
                                    "median": float(medians[c, r]) if np.any(valid) else None})
    uncentered = np.mean(cube, axis=2) if len(identifiers) else np.empty((0, 3))
    centered = np.mean(cube - medians[None, :, :], axis=2) if len(identifiers) else np.empty((0, 3))
    common_ledger = {"source_rows_by_context": ledger, "replicate_median_populations": populations,
                     "source_identifier_count": len(identifiers), "blocked_reasons": sorted(set(blocked)),
                     "missing_required_empirical_keys": [
                         {"id": identifier, "context": context, "replicate": replicate,
                          "reason": "unobserved_empirical_key" if (identifier, context, replicate) not in observations else "nonfinite_or_missing_source_value"}
                         for identifier in identifiers for context in CONTEXTS for replicate in (1, 2, 3)
                         if not math.isfinite(observations.get((identifier, context, replicate), math.nan))]}
    return (
        Endpoint(BASES[1], identifiers, uncentered, "log2(RNA/DNA)", blocked=reason, ledger=common_ledger),
        Endpoint(BASES[2], identifiers, centered, "log2(RNA/DNA)", blocked=reason, ledger=common_ledger),
    )


def moments(values):
    result = {}
    for column, context in enumerate(CONTEXTS):
        finite = values[np.isfinite(values[:, column]), column]
        n = len(finite)
        result[context] = {"finite_count": n, "missing_count": len(values) - n,
                           "mean": float(np.mean(finite)) if n else None,
                           "population_sd_ddof0": float(np.std(finite, ddof=0)) if n else None,
                           "sample_sd_ddof1": float(np.std(finite, ddof=1)) if n > 1 else None}
    return result


def row_sum_summary(endpoint):
    mask = np.isfinite(endpoint.values).all(axis=1)
    sums = endpoint.values[mask].sum(axis=1)
    result = {"complete_three_context_rows": len(sums), "mean": float(np.mean(sums)) if len(sums) else None,
              "maximum_absolute": float(np.max(np.abs(sums))) if len(sums) else None,
              "count_absolute_gt_1e-8": int(np.sum(np.abs(sums) > ABS_TOL))}
    if endpoint.raw is not None:
        exact = []
        with localcontext() as context:
            context.prec = 100
            for identifier, complete in zip(endpoint.ids, mask):
                if complete:
                    exact.append((identifier, sum((Decimal(v) for v in endpoint.raw[identifier]), Decimal(0))))
        witness = max(exact, key=lambda item: abs(item[1])) if exact else None
        result["exact_decimal"] = {"zero_count": sum(value == 0 for _, value in exact),
                                   "maximum_absolute": str(abs(witness[1])) if witness else None,
                                   "witness_id": witness[0] if witness else None}
    else:
        result["exact_decimal"] = {"applicable": False, "reason": "derived binary-float arithmetic"}
    return result


def metric(candidate, target, decimal_pairs=None):
    finite = np.isfinite(candidate) & np.isfinite(target)
    left, right = candidate[finite], target[finite]
    errors = left - right
    absolute = np.abs(errors)
    n = len(errors)
    if decimal_pairs is not None and len(decimal_pairs) != n:
        raise ValueError("exact source-cell count differs from paired finite cell count")
    equal = (sum(Decimal(a) == Decimal(b) for a, b in decimal_pairs) if decimal_pairs is not None else None)
    return {"estimable": bool(n), "nonestimable_reason": None if n else "no paired finite cells",
            "paired_finite_cells": n, "within_fixed_tolerance_count": int(np.sum(absolute <= ABS_TOL + REL_TOL * np.abs(right))),
            "signed_mean_error": float(np.mean(errors)) if n else None,
            "mae": float(np.mean(absolute)) if n else None,
            "rmse": float(np.sqrt(np.mean(errors ** 2))) if n else None,
            "maximum_absolute_error": float(np.max(absolute)) if n else None,
            "absolute_error_type7_quantiles": dict(zip(("0", "0.5", "0.95", "1"),
                np.quantile(absolute, [0, 0.5, 0.95, 1], method="linear").tolist())) if n else None,
            "exact_decimal": {"applicable": decimal_pairs is not None, "equal_count": equal,
                              "reason": None if decimal_pairs is not None else "not applicable to derived arithmetic"}}


def compare(candidate, target, old_ids, admitted_ids):
    if candidate.blocked or target.blocked:
        return {"estimable": False, "reason": candidate.blocked or target.blocked}
    all_ids = set(candidate.ids) | set(target.ids)
    selectors = {"full_source_population": all_ids,
                 "exact_old_s9_intersection": set(candidate.ids) & set(target.ids) & set(old_ids),
                 "fixed_old_admitted_intersection": set(candidate.ids) & set(target.ids) & set(admitted_ids)}
    reports = {}
    for label, selected in selectors.items():
        identifiers = sorted(selected)
        left, right = candidate.subset(identifiers), target.subset(identifiers)
        paired = np.isfinite(left) & np.isfinite(right)
        complete = paired.all(axis=1)
        complete_ids = [identifier for identifier, valid in zip(identifiers, complete) if valid]
        direct = candidate.raw is not None and target.raw is not None

        def pairs(mask, column=None):
            if not direct:
                return None
            return [(candidate.raw[identifier][c], target.raw[identifier][c])
                    for i, identifier in enumerate(identifiers) for c in range(3)
                    if mask[i, c] and (column is None or c == column)]

        available = {c: metric(left[:, j], right[:, j], pairs(paired, j)) for j, c in enumerate(CONTEXTS)}
        complete_mask = np.repeat(complete[:, None], 3, axis=1)
        complete_metrics = {c: metric(left[complete, j], right[complete, j], pairs(complete_mask, j))
                            for j, c in enumerate(CONTEXTS)}
        reports[label] = {"estimable": bool(paired.any()), "complete_three_context_estimable": bool(complete.any()),
                          "nonestimable_reason": None if paired.any() else "no paired finite source cells in this fixed cohort",
                          "universe_ids": len(identifiers), "universe_ids_sha256": ids_hash(identifiers),
                          "candidate_present_ids": sum(i in candidate.index for i in identifiers),
                          "target_present_ids": sum(i in target.index for i in identifiers),
                          "complete_three_context_paired_ids": len(complete_ids),
                          "complete_three_context_ids_sha256": ids_hash(complete_ids),
                          "incomplete_or_unmatched_ids": len(identifiers) - len(complete_ids),
                          "available_context_pair_metrics": available,
                          "available_context_pairs_pooled": metric(left.ravel(), right.ravel(), pairs(paired)),
                          "complete_three_context_metrics": complete_metrics,
                          "complete_three_context_pooled": metric(left[complete].ravel(), right[complete].ravel(), pairs(complete_mask)),
                          "candidate_moments_on_complete_comparison": moments(left[complete]),
                          "source_moments_on_complete_comparison": moments(right[complete])}
    return {"estimable": True, "error_direction": "candidate minus unchanged target", "cohorts": reports}


def pipeline(base, population, ddof, post):
    if population not in POPULATIONS or ddof not in (0, 1) or post not in POST:
        raise ValueError("operation is outside the frozen Cartesian grid")
    identifier = f"{base.name}__pre_{population}__ddof{ddof}__post_{post}"
    record = {"id": identifier, "activity_base": base.name, "pre_zscore_population": population,
              "pre_zscore_ddof": ddof, "post_ESS_step": post, "estimable": False}
    if base.blocked:
        record["nonestimable_reason"] = base.blocked
        return record, None
    complete = np.isfinite(base.values).all(axis=1)
    identifiers = [i for i, keep in zip(base.ids, complete) if keep]
    record["complete_source_population"] = {"count": len(identifiers), "ids_sha256": ids_hash(identifiers)}
    if not len(identifiers):
        record["nonestimable_reason"] = "no complete three-context source rows"
        return record, None
    means, sds, parameters = [], [], []
    for column, context in enumerate(CONTEXTS):
        mask = complete if population == POPULATIONS[0] else np.isfinite(base.values[:, column])
        values = base.values[mask, column]
        if len(values) <= ddof:
            record["nonestimable_reason"] = context + ": insufficient pre-zscore population"
            return record, None
        if np.all(values == values[0]):
            record["nonestimable_reason"] = context + ": nonpositive pre-zscore SD (exactly constant values)"
            return record, None
        mean, sd = float(np.mean(values)), float(np.std(values, ddof=ddof))
        if not math.isfinite(sd) or sd <= 0:
            record["nonestimable_reason"] = context + ": nonpositive pre-zscore SD"
            return record, None
        means.append(mean)
        sds.append(sd)
        parameters.append({"context": context, "count": len(values), "mean": mean, "sd": sd,
                           "ids_sha256": ids_hash([i for i, keep in zip(base.ids, mask) if keep])})
    record["pre_parameters"] = parameters
    standardized = (base.values[complete] - np.asarray(means)) / np.asarray(sds)
    values = standardized - np.mean(standardized, axis=1, keepdims=True)
    record["post_parameters"] = None
    if post != "none":
        post_ddof = 0 if post == POST[1] else 1
        if len(values) <= post_ddof:
            record["nonestimable_reason"] = "insufficient post-zscore source population"
            return record, None
        if np.any(np.all(values == values[0], axis=0)):
            record["nonestimable_reason"] = "nonpositive post-ESS column SD (exactly constant values)"
            return record, None
        post_mean, post_sd = np.mean(values, axis=0), np.std(values, axis=0, ddof=post_ddof)
        if not np.all(np.isfinite(post_sd) & (post_sd > 0)):
            record["nonestimable_reason"] = "nonpositive post-ESS column SD"
            return record, None
        record["post_parameters"] = {c: {"count": len(values), "mean": float(post_mean[j]),
                                          "sd": float(post_sd[j]), "ddof": post_ddof}
                                     for j, c in enumerate(CONTEXTS)}
        values = (values - post_mean) / post_sd
    record["estimable"] = True
    result = Endpoint(identifier, identifiers, values, "dimensionless derived compatibility score")
    record["moments"] = moments(values)
    record["row_sums"] = row_sum_summary(result)
    return record, result


def code_hypothesis(source_ess, s12, old_ids, admitted_ids):
    record = {"estimable": False, "operation": "S11 ESS restricted to every exact S12 ID, column sample zscore ddof1",
              "fixed_maximum_absolute_tolerance": 1e-12, "historical_pipeline_certified": False}
    if s12 is None or source_ess.blocked or s12.blocked:
        record["nonestimable_reason"] = source_ess.blocked or (s12.blocked if s12 else "final S12 workbook unavailable or version mismatch")
        return record, None
    identifiers = sorted(s12.ids)
    source = source_ess.subset(identifiers)
    record["normalization_population"] = {"count": len(identifiers), "ids_sha256": ids_hash(identifiers),
                                           "identifiers": identifiers}
    missing = [identifier for identifier, row in zip(identifiers, source) if not np.isfinite(row).all()]
    if missing or not np.isfinite(s12.subset(identifiers)).all() or len(identifiers) < 2:
        record["nonestimable_reason"] = "every unique S12 ID must have finite three-context source and target; no subset rescue"
        record["missing_or_incomplete_S11_ids"] = missing
        return record, None
    if np.any(np.all(source == source[0], axis=0)):
        record["nonestimable_reason"] = "nonpositive source SD (exactly constant values)"
        return record, None
    mean, sd = np.mean(source, axis=0), np.std(source, axis=0, ddof=1)
    if not np.all(np.isfinite(sd) & (sd > 0)):
        record["nonestimable_reason"] = "nonpositive or nonfinite source SD"
        return record, None
    result = Endpoint("pre_value_code_hypothesis", identifiers, (source - mean) / sd, "dimensionless derived compatibility score")
    record["estimable"] = True
    record["moments"] = moments(result.values)
    record["row_sums"] = row_sum_summary(result)
    record["normalization_parameters"] = {c: {"mean": float(mean[j]), "sample_sd_ddof1": float(sd[j])}
                                          for j, c in enumerate(CONTEXTS)}
    comparison = compare(result, s12, old_ids, admitted_ids)
    record["comparison"] = comparison
    record["fixed_1e-12_test_by_cohort"] = {}
    cohorts = {"full_s12_source_id_set": identifiers,
               "exact_old_s9_intersection": sorted(set(identifiers) & set(old_ids)),
               "fixed_old_admitted_intersection": sorted(set(identifiers) & set(admitted_ids))}
    for name, selected in cohorts.items():
        absolute = np.abs(result.subset(selected) - s12.subset(selected))
        maximum = float(np.max(absolute)) if absolute.size else None
        record["fixed_1e-12_test_by_cohort"][name] = {"ids": len(selected), "ids_sha256": ids_hash(selected),
            "maximum_absolute_error": maximum, "matched": maximum <= 1e-12 if maximum is not None else None,
            "within_1e-12_cells": int(np.sum(absolute <= 1e-12)), "cells": absolute.size}
    return record, result


def old_endpoint(path, name):
    # This existing stdlib XML reader preserves source decimal lexemes.
    from audit_context_source_arithmetic import source_rows
    rows, ledger = source_rows(path.read_bytes())
    identifiers = [row["id"] for row in rows]
    raw = {row["id"]: row["raw"] for row in rows}
    return Endpoint(name, identifiers, np.asarray([[float(v) for v in raw[i]] for i in identifiers]),
                    "dimensionless published source score", raw=raw, ledger=ledger)


def save_endpoint(path, endpoint):
    write_gzip(path, {"name": endpoint.name, "contexts": list(CONTEXTS), "units": endpoint.units,
                      "derived_values_not_source_replacements": endpoint.raw is None,
                      "blocked_reason": endpoint.blocked,
                      "ids": endpoint.ids,
                      "values": [[float(v) if math.isfinite(v) else None for v in row] for row in endpoint.values]})


def run(args):
    if args.out.exists():
        raise FileExistsError("retain earlier attempts; audit output must be new")
    if sha(args.protocol) != PROTOCOL_SHA256:
        raise ValueError("frozen recovery protocol SHA differs")
    protocol = json.loads(args.protocol.read_bytes())
    binding = json.loads(args.binding.read_bytes())
    if (binding.get("protocol_sha256") != PROTOCOL_SHA256
            or binding.get("status") != "bound_before_quantitative_inspection"
            or binding.get("semantic_binding_complete") is not True):
        raise ValueError("source schema gate is not bound")
    exported = json.loads(gzip.decompress(args.source.read_bytes()))
    if (exported.get("schema_version") != 1 or exported.get("protocol_sha256") != PROTOCOL_SHA256
            or exported.get("schema_binding_sha256") != sha(args.binding)):
        raise ValueError("source export does not bind this frozen schema/protocol")
    extra_inputs = []
    evidence_inputs = {}
    for label, receipt in binding.get("evidence", {}).items():
        path = Path(receipt["path"])
        if sha(path) != receipt["sha256"]:
            raise ValueError("source binding evidence receipt differs: " + label)
        extra_inputs.append(path)
        evidence_inputs[label] = path
    export_receipt = None
    if args.export_receipt is not None:
        export_receipt = json.loads(args.export_receipt.read_bytes())
        if (export_receipt.get("protocol_sha256") != PROTOCOL_SHA256
                or export_receipt.get("schema_binding_sha256") != sha(args.binding)
                or export_receipt["export"].get("sha256") != sha(args.source)
                or export_receipt["export"].get("bytes") != args.source.stat().st_size):
            raise ValueError("source export receipt binds different source bytes")
        extra_inputs.append(args.export_receipt)
        source_ledger = Path(export_receipt["ledger"]["path"])
        if (sha(source_ledger) != export_receipt["ledger"]["sha256"]
                or source_ledger.stat().st_size != export_receipt["ledger"]["bytes"]):
            raise ValueError("source cell-type/missing/formula ledger hash differs")
        extra_inputs.append(source_ledger)
    source = exported["source"]
    source_path = Path(source["path"])
    expected_s11 = next(v for v in protocol["pre_value_code_hypothesis"]["downstream_predicted_members"] if "Table S11" in v["path"])
    if (sha(source_path) != source["sha256"] or source["sha256"] != expected_s11["sha256"]
            or source_path.stat().st_size != source["bytes"] or source["bytes"] != expected_s11["bytes"]
            or any(binding["source_workbook"].get(k) != source[k] for k in ("sha256", "bytes"))):
        raise ValueError("official Table11 source/export/schema byte commitments differ")
    old_path = Path(protocol["old_inputs"]["s9_workbook"]["path"])
    panel_path = Path(protocol["old_inputs"]["old_admitted_panel"]["path"])
    if sha(old_path) != protocol["old_inputs"]["s9_workbook"]["sha256"] or sha(panel_path) != protocol["old_inputs"]["old_admitted_panel"]["sha256"]:
        raise ValueError("frozen old source/panel commitments differ")
    old = old_endpoint(old_path, "unchanged_old_S9")
    panel = json.loads(panel_path.read_bytes())
    admitted = [row["id"] for row in panel["records"]]
    if len(admitted) != len(set(admitted)) or len(admitted) != protocol["old_inputs"]["old_admitted_panel"]["element_count"]:
        raise ValueError("frozen admitted ID set differs")
    tables, roles = exported["tables"], binding["roles"]
    aggregate = direct_endpoint(BASES[0], tables, roles.get("activity_aggregate"), "log2(RNA/DNA)")
    published_ess = direct_endpoint("unchanged_final_S11_ESS", tables, roles.get("published_ess"), "dimensionless published ESS")
    means, median_means = replicate_bases(tables, roles.get("replicates"))
    bases = (aggregate, means, median_means)
    earlier_identity = binding.get("earlier_source_identity", {})
    s12_path = args.s12_workbook
    if s12_path is None and earlier_identity.get("table12_retained_path"):
        s12_path = Path(earlier_identity["table12_retained_path"])
    s12, s12_identity = None, {"provided_or_bound_path": s12_path is not None, "byte_identical_to_old_s9": None,
                              "official_member_identity_from_binding": earlier_identity.get("table12_folds_byte_identical_old_s9") is True,
                              "official_design_identity_from_binding": earlier_identity.get("table10_design_byte_identical_old_s7") is True}
    if s12_path is not None:
        s12_identity.update({"file": str(s12_path), "bytes": s12_path.stat().st_size,
                             "sha256": sha(s12_path)})
        s12_identity["byte_identical_to_old_s9"] = (s12_identity["sha256"] == sha(old_path)
            and earlier_identity.get("table12_sha256") == s12_identity["sha256"]
            and s12_identity["official_member_identity_from_binding"])
        if s12_identity["byte_identical_to_old_s9"]:
            s12 = old_endpoint(s12_path, "unchanged_final_S12")
        else:
            s12_identity["status"] = "source-version discrepancy; affected S11-to-S12 test blocked"

    args.out.mkdir(parents=True)
    (args.out / "pipelines").mkdir()
    bindings = {str(p): {"bytes": p.stat().st_size, "sha256": sha(p)}
                for p in [args.protocol, args.binding, args.source, source_path, old_path, panel_path, *extra_inputs]}
    if s12_path is not None:
        bindings[str(s12_path)] = {"bytes": s12_path.stat().st_size, "sha256": sha(s12_path)}
    inputs = {"protocol_sha256": PROTOCOL_SHA256, "files": bindings,
              "implementation_sha256": sha(Path(__file__)),
              "old_s9_reader_sha256": sha(Path(__file__).with_name("audit_context_source_arithmetic.py")),
              "runtime": {"python": platform.python_version(), "numpy": np.__version__},
              "source_reader_provenance": exported.get("strict_reader_source"),
              "source_exporter_implementation_sha256": export_receipt.get("implementation_sha256") if export_receipt else None,
              "source_schema_bound_before_numeric_checks": True,
              "exact_decimal_policy": "Only direct preserved source-cell comparisons; derived means/zscores are not exact decimal claims"}
    write_json(args.out / "inputs.json", inputs)
    shutil.copyfile(__file__, args.out / "implementation.py")
    shutil.copyfile(Path(__file__).with_name("audit_context_source_arithmetic.py"), args.out / "old-s9-reader.py")
    for source_file, target in ((args.protocol, "protocol.json"), (args.binding, "schema-binding.json"),
                                (args.source, "source-values.json.gz")):
        shutil.copyfile(source_file, args.out / target)
    if export_receipt is not None:
        shutil.copyfile(args.export_receipt, args.out / "source-export-receipt.json")
        shutil.copyfile(source_ledger, args.out / "source-cell-type-missingness-ledger.jsonl.gz")
    (args.out / "source-evidence").mkdir()
    for label, path in evidence_inputs.items():
        shutil.copyfile(path, args.out / "source-evidence" / (label + path.suffix))
    ledgers = {"raw_columns_and_all_stored_rows": "source-values.json.gz",
               "endpoint_rows": {e.name: e.ledger for e in (aggregate, published_ess, means)},
               "old_S9": old.ledger, "final_S12_identity": s12_identity,
               "fixed_old_admitted_ids": sorted(admitted), "fixed_old_admitted_ids_sha256": ids_hash(admitted),
               "all_exported_table_counts": {name: table.get("counts") for name, table in tables.items()}}
    write_gzip(args.out / "id-missingness-ledger.json.gz", ledgers)
    write_json(args.out / "code-hypothesis-population.json", {
        "committed_before_diagnostic_arithmetic": True,
        "definition": "Every exact unchanged final S12 source ID; no source-value-based subset selection",
        "identifiers": sorted(s12.ids) if s12 is not None else None,
        "count": len(s12.ids) if s12 is not None else None,
        "ids_sha256": ids_hash(s12.ids) if s12 is not None else None,
        "source_version_discrepancy": s12_identity})
    for endpoint in bases:
        save_endpoint(args.out / (endpoint.name + ".json.gz"), endpoint)
    variants = []
    for base, population, ddof, post in itertools.product(bases, POPULATIONS, (0, 1), POST):
        record, derived = pipeline(base, population, ddof, post)
        if derived is not None:
            relative = "pipelines/" + record["id"] + ".json.gz"
            save_endpoint(args.out / relative, derived)
            record["derived_values_file"] = relative
            record["comparisons"] = {"final_S11_ESS": compare(derived, published_ess, old.ids, admitted),
                                     "unchanged_old_S9": compare(derived, old, old.ids, admitted)}
        variants.append(record)
    if len(variants) != protocol["fixed_pipeline_grid"]["variant_count_if_all_fields_present"]:
        raise ValueError("frozen Cartesian grid cardinality differs")
    hypothesis, hypothesis_values = code_hypothesis(published_ess, s12, old.ids, admitted)
    if hypothesis_values is not None:
        save_endpoint(args.out / "code-hypothesis-values.json.gz", hypothesis_values)
    summary = {"schema_version": 1, "protocol_id": protocol["protocol_id"], "protocol_sha256": PROTOCOL_SHA256,
               "contexts": list(CONTEXTS), "source": source, "floating_tolerance": {"absolute": ABS_TOL, "relative": REL_TOL},
               "variant_count": len(variants), "estimable_variant_count": sum(v["estimable"] for v in variants),
               "variants": variants, "published_ESS_row_sums": row_sum_summary(published_ess) if not published_ess.blocked else {"nonestimable_reason": published_ess.blocked},
               "endpoint_populations": {e.name: {"source_ids": len(e.ids), "complete_three_context_ids": int(np.isfinite(e.values).all(axis=1).sum()),
                   "ids_sha256": ids_hash(e.ids), "blocked_reason": e.blocked, "moments": moments(e.values) if not e.blocked else None,
                   "ESS_row_sums": row_sum_summary(e) if not e.blocked and "dimensionless" in e.units else None}
                   for e in (*bases, published_ess, old)},
               "replicate_median_populations": means.ledger["replicate_median_populations"],
               "aggregate_comparisons": {base.name: compare(base, aggregate, old.ids, admitted) for base in bases[1:]},
               "unchanged_final_ESS_vs_old_S9": compare(published_ess, old, old.ids, admitted),
               "pre_value_code_hypothesis": hypothesis,
               "flags": {"empirical_raw_endpoint_binding": binding["semantic_binding_complete"],
                         "official_source_authenticity": "official-distribution source/schema receipts bound; SHA integrity is not publisher-signed authenticity",
                         "old_S9_release_relationship": s12_identity,
                         "replicate_aggregate_arithmetic": {"status": "all fixed comparisons reported; no aggregate selected by agreement",
                             "report": "summary.json#/aggregate_comparisons"},
                         "publishedESS_formula": {"status": "all 36 fixed routes enumerated; compatibility is not unique historical identification",
                             "report": "summary.json#/variants", "separate_code_hypothesis_report": "summary.json#/pre_value_code_hypothesis"},
                         "historical_pipeline_identified": False, "endpoint_replaced": False, "predictive_models_fitted": False,
                         "independent_experimental_validation": False},
               "limits": ["Every frozen arithmetic route is retained; no best endpoint is selected.",
                          "Numerical identity can establish file compatibility, not original producer operations.",
                          "Source normalization uses the whole source population, including future held-out outcomes.",
                          "Final release and old S9 are versions of one experiment, not independent biological replication.",
                          "Missing/duplicate empirical keys and every raw source token remain in the source and missingness ledgers."]}
    write_json(args.out / "summary.json", summary)
    files = {p.relative_to(args.out).as_posix(): {"bytes": p.stat().st_size, "sha256": sha(p)}
             for p in sorted(args.out.rglob("*")) if p.is_file()}
    write_json(args.out / "artifacts.json", {"schema_version": 1, "protocol_sha256": PROTOCOL_SHA256, "files": files})
    print(json.dumps({"output": str(args.out), "variants": len(variants), "estimable": summary["estimable_variant_count"],
                      "historical_pipeline_identified": False, "endpoint_replaced": False}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=Path("experiments/context-endpoint-recovery-v1.json"))
    parser.add_argument("--binding", type=Path, default=Path("research/context-endpoint-recovery/source/schema-binding.json"))
    parser.add_argument("--source", type=Path, default=Path("research/context-endpoint-recovery/source/jointdata-source-values.json.gz"))
    parser.add_argument("--s12-workbook", type=Path)
    parser.add_argument("--export-receipt", type=Path)
    parser.add_argument("--out", type=Path, default=Path("runs/recovered-context-primary"))
    run(parser.parse_args())


if __name__ == "__main__":
    main()
