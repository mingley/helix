#!/usr/bin/env python3
"""Independent NumPy replay of frozen published-context-specificity study v2.

This module neither imports the regulatory replication nor calls Rust. Before
fitting, it verifies the frozen protocol, every source file, raw XLSX cells,
declared reference intervals, complete-panel eligibility, families and guards.
Published ESS is used unchanged; it is not raw reporter activity.
"""
from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import datetime
import gzip
import hashlib
import itertools
import json
import math
import pathlib
import platform
import posixpath
import re
import shutil
import traceback
from typing import Mapping, Sequence
import xml.etree.ElementTree as ET
import zipfile

import numpy as np


BASES = "ACGT"
COMPLEMENT = str.maketrans(BASES, "TGCA")
MODELS = (
    "global_mean", "context_mean", "shared_dna", "additive_context_dna",
    "independent_context_dna",
)


def reverse_complement(sequence: str) -> str:
    """Return an exact reverse complement, rejecting ambiguous/nonuppercase DNA."""
    if not isinstance(sequence, str) or not sequence or set(sequence) - set(BASES):
        raise ValueError("sequence must be nonempty uppercase A/C/G/T")
    return sequence.translate(COMPLEMENT)[::-1]


def canonical_kmers(k: int = 3) -> tuple[str, ...]:
    if isinstance(k, bool) or not isinstance(k, int) or k <= 0:
        raise ValueError("k must be a positive integer")
    representatives = {
        min(word, reverse_complement(word))
        for bases in itertools.product(BASES, repeat=k)
        for word in ["".join(bases)]
    }
    return tuple(sorted(representatives))


KMER_ORDER = canonical_kmers(3)
FEATURE_NAMES = tuple(f"rc3mer_{word}" for word in KMER_ORDER) + ("gc_fraction",)
FROZEN_V2_SHA256 = "8f0576845fecb5cfcf9bfa354fab446d849244448a3489905a6ce33bdeb21cec"
ENDPOINT = "published_element_specificity_score"


@dataclass(frozen=True)
class ExcelNumber:
    literal: str

    def finite_float(self) -> float:
        value = float(self.literal)
        if not math.isfinite(value):
            raise ValueError("nonfinite source numeric cell")
        return value


@dataclass(frozen=True)
class ExcelError:
    text: str


class IndependentWorkbook:
    """Read XLSX cell types directly, without root importer or openpyxl.

    Strings retain case/whitespace. Numeric cells retain their XML literal;
    Excel errors and booleans stay distinct. Formulas are rejected rather than
    evaluated or silently accepted as original numeric measurements.
    """
    _main = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    _rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

    def __init__(self, path: pathlib.Path):
        self.archive = zipfile.ZipFile(path)
        try:
            members = self.archive.infolist()
            if len({item.filename for item in members}) != len(members):
                raise ValueError("XLSX has duplicate ZIP member names")
            if any(item.file_size > 512 * 1024 * 1024 for item in members):
                raise ValueError("XLSX member exceeds the 512MiB parser budget")
            self.strings = []
            if "xl/sharedStrings.xml" in self.archive.namelist():
                with self.archive.open("xl/sharedStrings.xml") as stream:
                    for event, node in ET.iterparse(stream, events=("end",)):
                        if node.tag == "{" + self._main + "}si":
                            self.strings.append("".join(part.text or "" for part in node.iter()
                                                        if part.tag == "{" + self._main + "}t"))
                            node.clear()
            relations = ET.fromstring(self.archive.read("xl/_rels/workbook.xml.rels"))
            targets = {}
            for relation in relations:
                if relation.attrib.get("TargetMode") == "External":
                    raise ValueError("external workbook relationship is unsupported")
                identity = relation.attrib["Id"]
                if identity in targets:
                    raise ValueError("duplicate workbook relationship")
                target = relation.attrib["Target"]
                if "\\" in target:
                    raise ValueError("noncanonical workbook relationship path")
                resolved = posixpath.normpath(target.lstrip("/") if target.startswith("/")
                                              else posixpath.join("xl", target))
                if not resolved.startswith("xl/") or resolved not in self.archive.namelist():
                    raise ValueError("workbook relationship escapes or misses ZIP member")
                targets[identity] = resolved
            workbook = ET.fromstring(self.archive.read("xl/workbook.xml"))
            self.sheets = {}
            for sheet in workbook.findall("{" + self._main + "}sheets/{" + self._main + "}sheet"):
                name = sheet.attrib["name"]
                if name in self.sheets:
                    raise ValueError("duplicate worksheet name")
                self.sheets[name] = targets[sheet.attrib["{" + self._rel + "}id"]]
        except BaseException:
            self.close()
            raise

    @staticmethod
    def column_index(address: str, row: int) -> int:
        match = re.fullmatch(r"([A-Z]+)([1-9][0-9]*)", address)
        if not match or int(match.group(2)) != row:
            raise ValueError("cell address disagrees with worksheet row")
        column = 0
        for letter in match.group(1):
            column = column * 26 + ord(letter) - ord("A") + 1
        return column - 1

    def rows(self, sheet_name: str):
        if sheet_name not in self.sheets:
            raise ValueError(f"expected source sheet is absent: {sheet_name}")
        prefix = "{" + self._main + "}"
        previous_row = 0
        with self.archive.open(self.sheets[sheet_name]) as stream:
            for event, node in ET.iterparse(stream, events=("end",)):
                if node.tag != prefix + "row":
                    continue
                row = int(node.attrib["r"])
                if row <= previous_row:
                    raise ValueError("worksheet rows must have increasing positive indices")
                previous_row = row
                cells = {}
                for cell in node.findall(prefix + "c"):
                    column = self.column_index(cell.attrib["r"], row)
                    if column in cells:
                        raise ValueError("duplicate cell address")
                    if cell.find(prefix + "f") is not None:
                        raise ValueError("formula cell requires an explicit source semantic decision")
                    kind = cell.attrib.get("t", "n")
                    value = cell.find(prefix + "v")
                    literal = None if value is None else value.text
                    if kind == "inlineStr":
                        parsed = "".join(part.text or "" for part in cell.iter()
                                         if part.tag == prefix + "t")
                    elif kind == "s":
                        if literal is None or not re.fullmatch(r"[0-9]+", literal):
                            raise ValueError("invalid shared-string index")
                        index = int(literal)
                        if index >= len(self.strings):
                            raise ValueError("shared-string index outside table")
                        parsed = self.strings[index]
                    elif kind == "n":
                        parsed = None if literal is None else ExcelNumber(literal)
                    elif kind == "b":
                        if literal not in {"0", "1"}:
                            raise ValueError("invalid boolean source cell")
                        parsed = literal == "1"
                    elif kind == "e":
                        parsed = ExcelError(literal or "")
                    elif kind == "str":
                        parsed = literal or ""
                    else:
                        raise ValueError(f"unsupported source cell type: {kind}")
                    cells[column] = parsed
                node.clear()
                yield row, cells

    def close(self):
        self.archive.close()

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        self.close()


def integer_sequence_features(sequence: str) -> dict:
    """Lossless source-derived counts before a single fixed division per feature."""
    reverse_complement(sequence)  # Validation is intentionally shared with RC.
    if len(sequence) < 3:
        raise ValueError("at least three unambiguous bases are required")
    counts = dict.fromkeys(KMER_ORDER, 0)
    for start in range(len(sequence) - 2):
        word = sequence[start:start + 3]
        counts[min(word, reverse_complement(word))] += 1
    return {"length": len(sequence), "gc_count": sequence.count("G") + sequence.count("C"),
            "canonical_trimer_counts": [counts[word] for word in KMER_ORDER]}


def sequence_features(sequence: str) -> np.ndarray:
    """32 RC-canonical overlapping 3-mer frequencies followed by GC fraction.

    Counts divide by L-2; GC divides by L. Complete A/C/G/T windows only are
    allowed. No orientation, locus, context, outcome or sequence length is a
    separate covariate. The 32 frequencies sum to one before floating roundoff.
    """
    counts = integer_sequence_features(sequence)
    return np.array(
        [count / (counts["length"] - 2) for count in counts["canonical_trimer_counts"]]
        + [counts["gc_count"] / counts["length"]],
        dtype=np.float64,
    )


def feature_matrix(sequences: Sequence[str]) -> np.ndarray:
    if not len(sequences):
        raise ValueError("at least one sequence is required")
    return np.vstack([sequence_features(sequence) for sequence in sequences])


def _matrix(values, name: str, columns: int | None = None) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if result.ndim != 2 or not np.isfinite(result).all():
        raise ValueError(f"{name} must be a finite two-dimensional array")
    if columns is not None and result.shape[1] != columns:
        raise ValueError(f"{name} has an incompatible feature count")
    return result


def _labels(values: Sequence[str], count: int) -> tuple[str, ...]:
    result = tuple(values)
    if len(result) != count or any(not isinstance(label, str) or not label for label in result):
        raise ValueError("context labels must be nonempty strings aligned to rows")
    return result


def _moments(x: np.ndarray, threshold: float) -> tuple[np.ndarray, np.ndarray]:
    if not math.isfinite(threshold) or threshold < 0:
        raise ValueError("scale threshold must be finite and nonnegative")
    mean = x.mean(axis=0)
    scale = np.sqrt(np.mean((x - mean) ** 2, axis=0))
    scale[scale <= threshold] = 1.0
    return mean, scale


def _penalized_solution(x: np.ndarray, y: np.ndarray, alpha: float) -> np.ndarray:
    if not x.shape[1]:
        return np.empty(0, dtype=np.float64)
    # Independent dense solve, rather than the Rust or previous SciPy cache.
    return np.linalg.solve(x.T @ x + alpha * np.eye(x.shape[1]), x.T @ y)


@dataclass
class RidgeFit:
    """A standardized linear fit; intercepts and feature scaling are unpenalized."""
    feature_means: np.ndarray
    feature_scales: np.ndarray
    standardized_coefficients: np.ndarray
    intercept: float
    context_intercepts: dict[str, float]
    training_rows: int
    alpha: float

    def predict(self, features, contexts: Sequence[str]) -> np.ndarray:
        x = _matrix(features, "prediction features", len(self.feature_means))
        labels = _labels(contexts, len(x))
        unknown = set(labels) - set(self.context_intercepts)
        if self.context_intercepts and unknown:
            raise ValueError(f"unseen contexts have no fitted intercept: {sorted(unknown)}")
        # Evaluate identical raw rows once. This is a function of input values,
        # not BLAS batch alignment; no rank rounding or epsilon is introduced.
        unique, inverse = np.unique(x, axis=0, return_inverse=True)
        linear = np.array([
            float(np.dot((row - self.feature_means) / self.feature_scales,
                         self.standardized_coefficients))
            for row in unique
        ])
        intercepts = np.array([
            self.context_intercepts[label] if self.context_intercepts else self.intercept
            for label in labels
        ])
        return linear[inverse] + intercepts

    def diagnostics(self) -> dict:
        return {
            "feature_means": self.feature_means.tolist(),
            "feature_scales": self.feature_scales.tolist(),
            "standardized_coefficients": self.standardized_coefficients.tolist(),
            "intercept": self.intercept,
            "context_intercepts": self.context_intercepts,
            "training_rows": self.training_rows,
            "alpha_summed_sse": self.alpha,
        }


@dataclass
class ContextModel:
    model: str
    shared: RidgeFit | None
    by_context: dict[str, RidgeFit]
    training_context_counts: dict[str, int]

    def predict(self, features, contexts: Sequence[str]) -> np.ndarray:
        if self.shared is not None:
            return self.shared.predict(features, contexts)
        x = _matrix(features, "prediction features")
        labels = _labels(contexts, len(x))
        unknown = set(labels) - set(self.by_context)
        if unknown:
            raise ValueError(f"unseen contexts have no fitted sequence model: {sorted(unknown)}")
        prediction = np.empty(len(x), dtype=np.float64)
        for context in sorted(set(labels)):
            indices = np.array([i for i, label in enumerate(labels) if label == context])
            prediction[indices] = self.by_context[context].predict(
                x[indices], [context] * len(indices))
        return prediction

    def diagnostics(self) -> dict:
        return {
            "model": self.model,
            "training_context_counts": self.training_context_counts,
            "shared": self.shared.diagnostics() if self.shared else None,
            "by_context": {label: fit.diagnostics() for label, fit in self.by_context.items()},
        }


def fit_context_model(features, targets, contexts: Sequence[str], model: str,
                      alpha: float = 10.0, scale_threshold: float = 1e-12) -> ContextModel:
    """Fit one generic row-SSE model, with alpha explicitly supplied by caller.

    The objective is sum_i residual_i^2 + alpha*sum_j beta_j^2. The additive
    model has a shared DNA slope and a free intercept per context. Both its X
    and y are residualized within context before the ridge solve. Separate fits
    each use their own training-only scaling and alpha. A future frozen study
    adapter must state whether pooled alpha is multiplied by the context count
    to give shared and separate fits comparable per-sequence regularization.
    """
    if model not in MODELS:
        raise ValueError(f"unknown context model: {model}")
    x = _matrix(features, "training features")
    y = np.asarray(targets, dtype=np.float64)
    if y.ndim != 1 or len(y) != len(x) or not len(y) or not np.isfinite(y).all():
        raise ValueError("training targets must be a nonempty finite aligned vector")
    if not math.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be finite and positive")
    labels = _labels(contexts, len(x))
    context_indices = {
        label: np.array([i for i, value in enumerate(labels) if value == label], dtype=np.int64)
        for label in sorted(set(labels))
    }
    counts = {label: len(indices) for label, indices in context_indices.items()}
    if model == "independent_context_dna":
        by_context = {
            label: fit_context_model(x[indices], y[indices], [label] * len(indices),
                                     "shared_dna", alpha, scale_threshold).shared
            for label, indices in context_indices.items()
        }
        return ContextModel(model, None, by_context, counts)
    mean, scale = _moments(x, scale_threshold)
    z = (x - mean) / scale
    beta = np.zeros(x.shape[1], dtype=np.float64)
    intercept = float(y.mean())
    intercepts = {}
    if model == "shared_dna":
        beta = _penalized_solution(z, y - intercept, alpha)
    elif model == "context_mean":
        intercepts = {label: float(y[indices].mean()) for label, indices in context_indices.items()}
    elif model == "additive_context_dna":
        residual_x = np.empty_like(z)
        residual_y = np.empty_like(y)
        context_moments = {}
        for label, indices in context_indices.items():
            mean_x = z[indices].mean(axis=0)
            mean_y = float(y[indices].mean())
            context_moments[label] = (mean_x, mean_y)
            residual_x[indices] = z[indices] - mean_x
            residual_y[indices] = y[indices] - mean_y
        beta = _penalized_solution(residual_x, residual_y, alpha)
        intercepts = {
            label: float(mean_y - np.dot(mean_x, beta))
            for label, (mean_x, mean_y) in context_moments.items()
        }
    fit = RidgeFit(mean, scale, beta, intercept, intercepts, len(x), alpha)
    return ContextModel(model, fit, {}, counts)


def fit_balanced_panel(features, target_matrix, context_order: Sequence[str],
                       alpha: float = 10.0, scale_threshold: float = 1e-12) -> dict[str, ContextModel]:
    """Fit provisional aligned-penalty models to a complete element x context panel.

    Each element occurs once in X and once in each column of Y. Shared slopes
    are fitted to the per-element context-mean target with alpha; separate
    slopes each use their context's target with the same alpha. This equals
    pooled shared ridge with context_count*alpha, not pooled alpha. S and A
    reuse exactly the same slope object, with global versus context intercepts.
    This helper is not permission to parse real outcomes before protocol freeze.
    """
    x = _matrix(features, "panel features")
    y = _matrix(target_matrix, "panel targets")
    labels = tuple(context_order)
    if not len(x) or len(x) != len(y) or len(labels) != y.shape[1] or not labels:
        raise ValueError("panel features/targets/contexts must be nonempty and aligned")
    if len(set(labels)) != len(labels) or any(not isinstance(label, str) or not label for label in labels):
        raise ValueError("panel context labels must be unique nonempty strings")
    mean_target = y.mean(axis=1)
    row_contexts = [labels[0]] * len(x)
    shared_fit = fit_context_model(x, mean_target, row_contexts, "shared_dna", alpha,
                                   scale_threshold).shared
    context_intercepts = {label: float(y[:, index].mean()) for index, label in enumerate(labels)}
    counts = {label: len(x) for label in labels}

    def copy_fit(beta, intercept, intercepts):
        return RidgeFit(shared_fit.feature_means.copy(), shared_fit.feature_scales.copy(),
                        beta, intercept, intercepts, len(x), alpha)

    global_intercept = float(mean_target.mean())
    zero = np.zeros(x.shape[1], dtype=np.float64)
    global_fit = copy_fit(zero.copy(), global_intercept, {})
    context_fit = copy_fit(zero.copy(), global_intercept, context_intercepts.copy())
    additive_fit = copy_fit(shared_fit.standardized_coefficients, global_intercept,
                            context_intercepts.copy())
    separate = {
        label: fit_context_model(x, y[:, index], [label] * len(x), "shared_dna", alpha,
                                 scale_threshold).shared
        for index, label in enumerate(labels)
    }
    return {
        "global_mean": ContextModel("global_mean", global_fit, {}, counts),
        "context_mean": ContextModel("context_mean", context_fit, {}, counts),
        "shared_dna": ContextModel("shared_dna", shared_fit, {}, counts),
        "additive_context_dna": ContextModel("additive_context_dna", additive_fit, {}, counts),
        "independent_context_dna": ContextModel("independent_context_dna", None, separate, counts),
    }


def paired_panel_indices(element_ids: Sequence[str], contexts: Sequence[str],
                         context_order: Sequence[str]) -> dict[str, tuple[int, ...]]:
    """Return only complete, unique context panels without imputing missing rows.

    Unexpected contexts and duplicate element/context cells are errors. Missing
    cells are omitted from the complete-panel result, never filled with zero.
    Caller retains and accounts for all source rows independently.
    """
    elements = tuple(element_ids)
    labels = _labels(contexts, len(elements))
    expected = _labels(context_order, len(context_order))
    if not expected or len(set(expected)) != len(expected):
        raise ValueError("context order must contain unique nonempty labels")
    cells: dict[str, dict[str, int]] = {}
    for index, (element, label) in enumerate(zip(elements, labels)):
        if not isinstance(element, str) or not element:
            raise ValueError("element IDs must be nonempty strings")
        if label not in expected:
            raise ValueError(f"unexpected context: {label}")
        if label in cells.setdefault(element, {}):
            raise ValueError(f"duplicate element/context cell: {element}/{label}")
        cells[element][label] = index
    return {element: tuple(by_context[label] for label in expected)
            for element, by_context in sorted(cells.items())
            if set(by_context) == set(expected)}


@dataclass(frozen=True)
class SequenceElement:
    element_id: str
    sequence: str
    assembly: str | None = None
    chromosome: str | None = None
    start0: int | None = None
    end0: int | None = None


def identity_overlap_families(elements: Sequence[SequenceElement]) -> dict[str, str]:
    """Transitive families from RC identity or positive genomic interval overlap.

    Genomic coordinates are zero-based half-open; touching intervals do not
    overlap. Coordinates must be all present or all absent. Identity can link
    different loci, while overlap requires matching assembly and chromosome.
    A family's stable ID is its lexicographically smallest source element ID.
    No outcome or sequence-feature similarity is consulted.
    """
    ids = [element.element_id for element in elements]
    if len(set(ids)) != len(ids) or any(not isinstance(value, str) or not value for value in ids):
        raise ValueError("source element IDs must be unique nonempty strings")
    parent = {value: value for value in ids}

    def find(value):
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def merge(left, right):
        a, b = find(left), find(right)
        parent[max(a, b)] = min(a, b)

    by_identity = {}
    intervals = {}
    for element in elements:
        rc = reverse_complement(element.sequence)
        canonical = min(element.sequence, rc)
        if canonical in by_identity:
            merge(element.element_id, by_identity[canonical])
        else:
            by_identity[canonical] = element.element_id
        coordinates = (element.assembly, element.chromosome, element.start0, element.end0)
        if all(value is None for value in coordinates):
            continue
        if any(value is None for value in coordinates):
            raise ValueError("genomic assembly/chromosome/start/end must be all present or absent")
        if not isinstance(element.assembly, str) or not element.assembly:
            raise ValueError("assembly must be nonempty")
        if not isinstance(element.chromosome, str) or not element.chromosome:
            raise ValueError("chromosome must be nonempty")
        if any(isinstance(value, bool) or not isinstance(value, int)
               for value in (element.start0, element.end0)):
            raise ValueError("genomic interval bounds must be integers")
        if not 0 <= element.start0 < element.end0:
            raise ValueError("genomic intervals must be nonempty and nonnegative")
        intervals.setdefault((element.assembly, element.chromosome), []).append(
            (element.start0, element.end0, element.element_id))
    for chromosome_intervals in intervals.values():
        component_id = None
        maximum_end = -1
        for start, end, element_id in sorted(chromosome_intervals):
            if start < maximum_end:
                merge(component_id, element_id)
                maximum_end = max(maximum_end, end)
            else:
                component_id, maximum_end = element_id, end
    return {element_id: find(element_id) for element_id in sorted(ids)}


def grouped_folds(row_families: Sequence[str], family_fold: Mapping[str, int],
                  fold_count: int) -> list[dict]:
    """Apply an externally specified family assignment; all contexts stay together."""
    families = tuple(row_families)
    if not families or any(not isinstance(value, str) or not value for value in families):
        raise ValueError("row families must be nonempty strings")
    if isinstance(fold_count, bool) or not isinstance(fold_count, int) or fold_count < 2:
        raise ValueError("at least two folds are required")
    if set(families) != set(family_fold):
        raise ValueError("family assignment must exactly cover observed families")
    if any(isinstance(value, bool) or not isinstance(value, int) or not 0 <= value < fold_count
           for value in family_fold.values()):
        raise ValueError("family fold index is invalid")
    folds = []
    for fold in range(fold_count):
        test = [i for i, family in enumerate(families) if family_fold[family] == fold]
        train = [i for i, family in enumerate(families) if family_fold[family] != fold]
        if not train or not test:
            raise ValueError("every fold must have train and test observations")
        train_families = sorted({families[i] for i in train})
        test_families = sorted({families[i] for i in test})
        if set(train_families) & set(test_families):
            raise ValueError("family leakage across train and test")
        folds.append({"fold_index_0based": fold, "train_indices": train, "test_indices": test,
                      "train_families": train_families, "test_families": test_families})
    return folds


def reconstruct_grouping(records: Sequence[dict], *, bin_bp: int = 1_000_000,
                          candidate_budget: int = 100_000_000, fold_count: int = 5,
                          guard_bp: int = 10_000, minimum_components: int = 100,
                          largest_component_fraction: float = 0.4,
                          minimum_training_fraction: float = 0.5) -> dict:
    """Independent full relation/guard reconstruction, without outcome access.

    Hamming verification uses NumPy byte comparisons, not root's packed-integer
    code. Guard lookup merges expanded test intervals and searches their ends,
    rather than root's starts/prefix-maximum algorithm. Candidate generation is
    the frozen exhaustive pigeonhole relation; all source IDs count in budget.
    """
    rows = sorted(records, key=lambda row: row["id"])
    ids = [row["id"] for row in rows]
    if not rows or len(set(ids)) != len(ids) or any(not isinstance(value, str) or not value or "\n" in value for value in ids):
        raise ValueError("grouping needs unique nonempty newline-free IDs")
    if bin_bp <= 0 or candidate_budget < 0 or fold_count < 2 or guard_bp < 0:
        raise ValueError("invalid grouping constants")
    parent = list(range(len(rows)))

    def representative(index):
        trail = []
        while index != parent[index]:
            trail.append(index)
            index = parent[index]
        for visited in trail:
            parent[visited] = index
        return index

    def connect(left, right):
        left, right = representative(left), representative(right)
        parent[max(left, right)] = min(left, right)

    first_parent, first_bin, first_exact = {}, {}, {}
    sequences, reverses, arrays, reverse_arrays = [], [], [], []
    for index, row in enumerate(rows):
        sequence = row["sequence"]
        reverse = reverse_complement(sequence)
        left, right = row["start"], row["end"]
        if (any(isinstance(bound, bool) or not isinstance(bound, int) for bound in (left, right))
                or left < 0 or right - left != len(sequence) or len(sequence) < 3
                or not isinstance(row["chromosome"], str) or not row["chromosome"]):
            raise ValueError("invalid exact genomic grouping record")
        sequences.append(sequence)
        reverses.append(reverse)
        arrays.append(np.frombuffer(sequence.encode("ascii"), dtype=np.uint8))
        reverse_arrays.append(np.frombuffer(reverse.encode("ascii"), dtype=np.uint8))
        documented_parent = row.get("parent_id") or row["id"]
        if documented_parent in first_parent:
            connect(index, first_parent[documented_parent])
        first_parent.setdefault(documented_parent, index)
        for genomic_bin in range(left // bin_bp, (right - 1) // bin_bp + 1):
            key = (row["chromosome"], genomic_bin)
            if key in first_bin:
                connect(index, first_bin[key])
            first_bin.setdefault(key, index)
        canonical = min(sequence, reverse)
        if canonical in first_exact:
            connect(index, first_exact[canonical])
        first_exact.setdefault(canonical, index)
    buckets = defaultdict(set)
    candidate_count = 0
    for index, sequence in enumerate(sequences):
        length = len(sequence)
        k = length // 10
        boundaries = [part * length // (k + 1) for part in range(k + 2)]
        candidates = set()
        for part, (start, end) in enumerate(zip(boundaries, boundaries[1:])):
            candidates.update(buckets.get((length, part, sequence[start:end]), ()))
        candidate_count += len(candidates)
        if candidate_count > candidate_budget:
            raise ValueError("independent exhaustive Hamming candidate budget exceeded")
        for peer in sorted(candidates):
            if representative(index) == representative(peer):
                continue
            forward_distance = int(np.count_nonzero(arrays[index] != arrays[peer]))
            if forward_distance <= k:
                connect(index, peer)
            elif int(np.count_nonzero(arrays[index] != reverse_arrays[peer])) <= k:
                connect(index, peer)
        for oriented in {sequence, reverses[index]}:
            for part, (start, end) in enumerate(zip(boundaries, boundaries[1:])):
                buckets[(length, part, oriented[start:end])].add(index)
    members = defaultdict(list)
    for index, identifier in enumerate(ids):
        members[representative(index)].append(identifier)
    components = [{"group": hashlib.sha256("\n".join(member_ids).encode()).hexdigest(),
                   "ids": member_ids, "size": len(member_ids)} for member_ids in members.values()]
    components.sort(key=lambda item: (-item["size"], item["group"], item["ids"]))
    if len(components) < minimum_components or components[0]["size"] / len(rows) > largest_component_fraction:
        raise ValueError("independent component feasibility gate failed")
    loads = [0] * fold_count
    assignment = {}
    group_by_id = {}
    for component in components:
        selected = min(range(fold_count), key=lambda fold: (loads[fold], fold))
        loads[selected] += component["size"]
        assignment[component["group"]] = selected
        group_by_id.update({identifier: component["group"] for identifier in component["ids"]})
    folds = []
    for fold in range(fold_count):
        test_groups = {group for group, selected in assignment.items() if selected == fold}
        intervals = defaultdict(list)
        for row in rows:
            if group_by_id[row["id"]] in test_groups:
                intervals[row["chromosome"]].append((row["start"] - guard_bp, row["end"] + guard_bp))
        merged = {}
        for chromosome, spans in intervals.items():
            joined = []
            for left, right in sorted(spans):
                if joined and left <= joined[-1][1]:
                    joined[-1] = (joined[-1][0], max(right, joined[-1][1]))
                else:
                    joined.append((left, right))
            merged[chromosome] = ([left for left, _ in joined], [right for _, right in joined])
        guard_groups = set()
        for row in rows:
            group = group_by_id[row["id"]]
            if group in test_groups or row["chromosome"] not in merged:
                continue
            starts, ends = merged[row["chromosome"]]
            first = bisect_left(ends, row["start"])
            if first < len(starts) and starts[first] <= row["end"]:
                guard_groups.add(group)
        train_ids = [identifier for identifier in ids if group_by_id[identifier] not in test_groups | guard_groups]
        test_ids = [identifier for identifier in ids if group_by_id[identifier] in test_groups]
        guard_ids = [identifier for identifier in ids if group_by_id[identifier] in guard_groups]
        if not test_ids or len(train_ids) / len(ids) < minimum_training_fraction:
            raise ValueError("independent guarded training-support gate failed")
        folds.append({"fold": fold, "test_groups": sorted(test_groups), "guard_groups": sorted(guard_groups),
                      "train_ids": train_ids, "test_ids": test_ids, "guard_ids": guard_ids})
    return {"records": [{"id": identifier, "group": group_by_id[identifier]} for identifier in ids],
            "components": components, "folds": folds, "candidate_pairs": candidate_count,
            "algorithm": "frozen pigeonhole candidates; NumPy exact Hamming; merged-expanded-interval guard"}


def errors(actual, predicted) -> dict:
    """Primary scalar losses only; exact-rank correlation is deliberately omitted."""
    y = np.asarray(actual, dtype=np.float64)
    prediction = np.asarray(predicted, dtype=np.float64)
    if y.ndim != 1 or y.shape != prediction.shape or not len(y):
        raise ValueError("loss vectors must be nonempty and aligned")
    if not np.isfinite(y).all() or not np.isfinite(prediction).all():
        raise ValueError("loss vectors must be finite")
    residual = y - prediction
    return {"count": len(y), "mse": float(np.mean(residual ** 2)),
            "mae": float(np.mean(np.abs(residual)))}


class SplitMix64:
    def __init__(self, seed: int):
        if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**64:
            raise ValueError("SplitMix64 seed must be uint64")
        self.state = seed

    def next_u64(self) -> int:
        mask = (1 << 64) - 1
        self.state = (self.state + 0x9E3779B97F4A7C15) & mask
        mixed = self.state
        mixed = ((mixed ^ (mixed >> 30)) * 0xBF58476D1CE4E5B9) & mask
        mixed = ((mixed ^ (mixed >> 27)) * 0x94D049BB133111EB) & mask
        return mixed ^ (mixed >> 31)


def paired_component_bootstrap(group_ids: Sequence[str], sequence_mse, sequence_mae,
                               comparisons: Sequence[tuple[str, str]], replicates: int,
                               seed: int) -> dict:
    """Common component draws; element-weighted ratio primary, group sensitivity.

    Each loss row already averages the observed contexts of one complete-panel
    element. Sampling a component repeats all its elements, with multiplicity.
    No outcomes are transformed and no models are refitted during resampling.
    """
    groups = tuple(group_ids)
    mse = _matrix(sequence_mse, "per-element MSE", len(MODELS))
    mae = _matrix(sequence_mae, "per-element MAE", len(MODELS))
    if len(groups) != len(mse) or mse.shape != mae.shape or not len(mse):
        raise ValueError("component IDs and per-element losses must align")
    if any(not isinstance(group, str) or not group for group in groups):
        raise ValueError("component IDs must be nonempty strings")
    if np.any(mse < 0) or np.any(mae < 0):
        raise ValueError("losses must be nonnegative")
    if isinstance(replicates, bool) or not isinstance(replicates, int) or replicates < 1:
        raise ValueError("bootstrap replicate count must be positive")
    order = sorted(set(groups))
    if len(order) * replicates > 50_000_000:
        raise ValueError("bootstrap index table exceeds 50-million-index budget")
    pairs = tuple(comparisons)
    if not pairs or len(set(pairs)) != len(pairs):
        raise ValueError("model comparisons must be nonempty and unique")
    for baseline, candidate in pairs:
        if baseline not in MODELS or candidate not in MODELS or baseline == candidate:
            raise ValueError("invalid model contrast")
    members = defaultdict(list)
    for index, group in enumerate(groups):
        members[group].append(index)
    indices = {group: np.array(members[group], dtype=np.int64) for group in order}
    counts = np.array([len(indices[group]) for group in order], dtype=np.float64)
    mse_sums = np.vstack([mse[indices[group]].sum(axis=0) for group in order])
    mae_sums = np.vstack([mae[indices[group]].sum(axis=0) for group in order])
    mean_mse = mse_sums / counts[:, None]
    mean_mae = mae_sums / counts[:, None]
    rng = SplitMix64(seed)
    draws = np.array([[rng.next_u64() % len(order) for _ in order]
                      for _ in range(replicates)], dtype=np.int32)
    sampled_sequence_mse = np.empty((replicates, len(MODELS)))
    sampled_sequence_mae = np.empty_like(sampled_sequence_mse)
    sampled_group_mse = np.empty_like(sampled_sequence_mse)
    sampled_group_mae = np.empty_like(sampled_sequence_mse)
    for replicate, selected in enumerate(draws):
        denominator = counts[selected].sum()
        sampled_sequence_mse[replicate] = mse_sums[selected].sum(axis=0) / denominator
        sampled_sequence_mae[replicate] = mae_sums[selected].sum(axis=0) / denominator
        sampled_group_mse[replicate] = mean_mse[selected].mean(axis=0)
        sampled_group_mae[replicate] = mean_mae[selected].mean(axis=0)

    def summary(model_mse, model_mae, original_mse, original_mae):
        contrast_mse = np.column_stack([model_mse[:, MODELS.index(baseline)]
                                       - model_mse[:, MODELS.index(candidate)]
                                       for baseline, candidate in pairs])
        contrast_mae = np.column_stack([model_mae[:, MODELS.index(baseline)]
                                       - model_mae[:, MODELS.index(candidate)]
                                       for baseline, candidate in pairs])
        summaries = []
        for j, (baseline, candidate) in enumerate(pairs):
            left, right = MODELS.index(baseline), MODELS.index(candidate)
            summaries.append({
                "comparison": {"baseline": baseline, "candidate": candidate},
                "observed_mse_improvement": float(original_mse[left] - original_mse[right]),
                "observed_mae_improvement": float(original_mae[left] - original_mae[right]),
                "mse_interval": np.quantile(contrast_mse[:, j], [0.025, 0.975], method="linear").tolist(),
                "mae_interval": np.quantile(contrast_mae[:, j], [0.025, 0.975], method="linear").tolist(),
            })
        return {"sampled_model_mses": model_mse.tolist(), "sampled_model_maes": model_mae.tolist(),
                "sampled_mse_improvements": contrast_mse.tolist(),
                "sampled_mae_improvements": contrast_mae.tolist(), "comparisons": summaries}

    return {
        "group_order": order, "model_order": list(MODELS), "group_sequence_counts": counts.astype(int).tolist(),
        "comparison_order": [{"baseline": baseline, "candidate": candidate} for baseline, candidate in pairs],
        "replicates": replicates, "seed": seed, "interval_probabilities": [0.025, 0.975],
        "resampling_indices": draws.tolist(),
        "equal_sequences": summary(sampled_sequence_mse, sampled_sequence_mae, mse.mean(axis=0), mae.mean(axis=0)),
        "equal_groups": summary(sampled_group_mse, sampled_group_mae, mean_mse.mean(axis=0), mean_mae.mean(axis=0)),
        "quantile_method": "Type7 linear interpolation; no significance or novelty decision",
        "interpretation_limit": "Conditional fixed-OOF-loss resampling of selected components; no model-refit or independent-study uncertainty",
    }


def evaluate_complete_panel(records: Sequence[dict], contexts: Sequence[str], folds: Sequence[dict],
                            alpha: float = 10.0) -> dict:
    """Evaluate a supplied synthetic/authorized panel without source loading.

    Folds must partition each element into train/test/whole-family guard, and
    every element must be tested exactly once. All contexts share these folds.
    """
    rows = sorted(records, key=lambda row: row["id"])
    ids = [row["id"] for row in rows]
    if not rows or len(set(ids)) != len(ids):
        raise ValueError("panel IDs must be nonempty and unique")
    labels = tuple(contexts)
    if not labels or len(set(labels)) != len(labels) or any(not isinstance(label, str) or not label for label in labels):
        raise ValueError("contexts must be unique nonempty strings")
    lookup = {identifier: index for index, identifier in enumerate(ids)}
    group_by_id = {row["id"]: row["group"] for row in rows}
    if any(not isinstance(group, str) or not group for group in group_by_id.values()):
        raise ValueError("panel groups must be nonempty strings")
    x = feature_matrix([row["sequence"] for row in rows])
    for row in rows:
        if set(row["observations"]) != set(labels):
            raise ValueError("every retained element must have exactly the frozen contexts")
        if any(isinstance(row["observations"][context], bool) for context in labels):
            raise ValueError("boolean source values cannot be observations")
    y = _matrix([[row["observations"][context] for context in labels] for row in rows], "outcome panel")
    predictions = np.full((len(rows), len(labels), len(MODELS)), np.nan)
    coverage = np.zeros(len(rows), dtype=np.int64)
    fold_index = np.full(len(rows), -1, dtype=np.int64)
    fits = []
    max_beta_identity = 0.0
    if len(folds) < 2:
        raise ValueError("OOF evaluation requires at least two folds")
    for number, fold in enumerate(folds):
        partitions = [fold[key] for key in ("train_ids", "test_ids", "guard_ids")]
        if any(len(values) != len(set(values)) for values in partitions):
            raise ValueError("duplicate element in a fold partition")
        if set().union(*map(set, partitions)) != set(ids) or sum(map(len, partitions)) != len(ids):
            raise ValueError("train/test/guard partitions must exactly cover panel IDs once")
        if not partitions[0] or not partitions[1]:
            raise ValueError("each fold requires nonempty training and test sets")
        family_sets = [{group_by_id[identifier] for identifier in values} for values in partitions]
        if any(family_sets[left] & family_sets[right] for left in range(3) for right in range(left + 1, 3)):
            raise ValueError("sequence family leaks across train/test/guard partitions")
        train = np.array(sorted(lookup[identifier] for identifier in partitions[0]), dtype=np.int64)
        test = np.array(sorted(lookup[identifier] for identifier in partitions[1]), dtype=np.int64)
        fitted = fit_balanced_panel(x[train], y[train], labels, alpha)
        for model_index, model in enumerate(MODELS):
            for context_index, context in enumerate(labels):
                predictions[test, context_index, model_index] = fitted[model].predict(
                    x[test], [context] * len(test))
        independent_beta = np.vstack([
            fitted["independent_context_dna"].by_context[context].standardized_coefficients
            for context in labels])
        max_beta_identity = max(max_beta_identity, float(np.max(np.abs(
            independent_beta.mean(axis=0) - fitted["shared_dna"].shared.standardized_coefficients))))
        coverage[test] += 1
        fold_index[test] = number
        fits.append({"fold_index_0based": number, "train_ids": [ids[i] for i in train],
                     "test_ids": [ids[i] for i in test], "guard_ids": sorted(partitions[2]),
                     "models": {model: fitted[model].diagnostics() for model in MODELS}})
    if not np.all(coverage == 1) or not np.isfinite(predictions).all():
        raise ValueError("test folds must cover every panel element exactly once")
    squared = (y[:, :, None] - predictions) ** 2
    absolute = np.abs(y[:, :, None] - predictions)
    per_sequence_mse = squared.mean(axis=1)
    per_sequence_mae = absolute.mean(axis=1)
    groups = sorted(set(group_by_id.values()))
    group_indices = defaultdict(list)
    for index, row in enumerate(rows):
        group_indices[row["group"]].append(index)
    group_indices = {group: group_indices[group] for group in groups}
    per_group = [{"group": group, "sequences": len(indices),
                  "model_losses": {"mse": per_sequence_mse[indices].mean(axis=0).tolist(),
                                   "mae": per_sequence_mae[indices].mean(axis=0).tolist()}}
                 for group, indices in group_indices.items()]
    summaries = []
    for index, model in enumerate(MODELS):
        summaries.append({"model": model, "sequence_weighted_mse": float(per_sequence_mse[:, index].mean()),
                          "sequence_weighted_mae": float(per_sequence_mae[:, index].mean()),
                          "macro_group_mse": float(np.mean([group["model_losses"]["mse"][index] for group in per_group])),
                          "macro_group_mae": float(np.mean([group["model_losses"]["mae"][index] for group in per_group])),
                          "by_context": [{"context": context, "mse": float(squared[:, j, index].mean()),
                                          "mae": float(absolute[:, j, index].mean())}
                                         for j, context in enumerate(labels)]})
    contrast_rows = []
    pair_improvement = []
    for left in range(len(labels)):
        for right in range(left + 1, len(labels)):
            actual = y[:, right] - y[:, left]
            predicted = predictions[:, right] - predictions[:, left]
            contrast_sq = (actual[:, None] - predicted) ** 2
            contrast_abs = np.abs(actual[:, None] - predicted)
            pair_improvement.append(contrast_sq[:, 3] - contrast_sq[:, 4])
            contrast_rows.append({"context_left": labels[left], "context_right": labels[right],
                                  "observed": actual.tolist(), "predictions": predicted.tolist(),
                                  "sequence_weighted_mse": contrast_sq.mean(axis=0).tolist(),
                                  "sequence_weighted_mae": contrast_abs.mean(axis=0).tolist()})
    factor = (len(labels) - 1) / (2 * len(labels))
    contrast_residual = 0.0 if not pair_improvement else float(np.max(np.abs(
        per_sequence_mse[:, 3] - per_sequence_mse[:, 4]
        - factor * np.vstack(pair_improvement).mean(axis=0))))
    prediction_rows = [{"fold_index_0based": int(fold_index[i]), "sequence_id": row["id"],
                        "group": row["group"], "context": context, "observed": float(y[i, j]),
                        "predictions": {model: float(predictions[i, j, k]) for k, model in enumerate(MODELS)}}
                       for i, row in enumerate(rows) for j, context in enumerate(labels)]
    return {
        "schema_version": 1, "endpoint": ENDPOINT, "contexts": list(labels), "model_order": list(MODELS),
        "alpha": alpha, "sequence_count": len(rows), "group_count": len(groups),
        "measurement_count": len(rows) * len(labels), "sequence_ids": ids,
        "feature_names": list(FEATURE_NAMES), "features": x.tolist(), "folds": fits,
        "predictions": prediction_rows, "per_sequence_errors": [
            {"sequence_id": row["id"], "group": row["group"],
             "model_losses": {"mse": per_sequence_mse[i].tolist(), "mae": per_sequence_mae[i].tolist()}}
            for i, row in enumerate(rows)], "per_group_errors": per_group, "models": summaries,
        "context_contrasts": contrast_rows,
        "balanced_design_checks": {
            "max_abs_mean_independent_beta_minus_shared_beta": max_beta_identity,
            "max_abs_mean_independent_prediction_minus_shared_prediction": float(np.max(np.abs(
                predictions[:, :, 4].mean(axis=1) - predictions[:, 0, 2]))),
            "max_abs_per_sequence_mse_contrast_identity_residual": contrast_residual},
        "fixed_zero_prediction_audit": {
            "sequence_weighted_mse": float(np.mean(y ** 2)), "sequence_weighted_mae": float(np.mean(np.abs(y))),
            "max_abs_per_sequence_context_mean": float(np.max(np.abs(y.mean(axis=1)))),
            "mean_by_context": {context: float(y[:, index].mean()) for index, context in enumerate(labels)},
            "scope": "Published ESS target descriptive audit; no zero-sum constraint is imposed"},
    }


def file_hashes(path: pathlib.Path) -> dict:
    sha, md5 = hashlib.sha256(), hashlib.md5(usedforsecurity=False)
    size = 0
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            size += len(chunk)
            sha.update(chunk)
            md5.update(chunk)
    return {"bytes": size, "sha256": sha.hexdigest(), "md5": md5.hexdigest()}


def checked_protocol(path: pathlib.Path) -> tuple[dict, str]:
    raw = path.read_bytes()
    observed = hashlib.sha256(raw).hexdigest()
    if observed != FROZEN_V2_SHA256:
        raise ValueError("only the separately frozen context-specificity-v2 protocol is accepted")
    protocol = json.loads(raw)
    if (protocol["status"] != "frozen" or protocol["protocol_id"] != "context-specificity-v2"
            or protocol["features"]["names"] != list(FEATURE_NAMES)
            or not protocol["frozen_before_ess_numeric_value_inspection_and_model_fitting"]):
        raise ValueError("frozen protocol identity/features/gate mismatch")
    return protocol, observed


def _source_integer(value) -> int:
    literal = value.literal if isinstance(value, ExcelNumber) else value
    if not isinstance(literal, str):
        raise ValueError("missing_or_invalid_coordinate")
    try:
        decimal = Decimal(literal)
        if not decimal.is_finite() or decimal != decimal.to_integral_value():
            raise ValueError("missing_or_invalid_coordinate")
        return int(decimal)
    except InvalidOperation as error:
        raise ValueError("missing_or_invalid_coordinate") from error


def verify_raw_workbooks(panel_directory: pathlib.Path, protocol: dict, normalized_records: Sequence[dict],
                         reference_path: pathlib.Path) -> tuple[list[dict], dict]:
    """Independently reconstruct all eligible raw S7/S9 rows, then compare panel.

    Exact IDs, source cell types/literals, source orientation and every reference
    interval are verified. Membership is reconstructed from all source rows;
    checking only the retained normalized rows would miss favorable exclusions.
    """
    from twobit_reference import TwoBitReader

    contexts = protocol["source"]["contexts"]
    source_paths = {}
    commitments = {}
    for key in ("design_workbook", "folds_workbook"):
        declared = protocol["source"][key]
        path = panel_directory / "raw" / pathlib.Path(declared["local_path"]).name
        observed = file_hashes(path)
        if observed != {"bytes": declared["size_bytes"], "sha256": declared["sha256"], "md5": declared["md5"]}:
            raise ValueError(f"frozen source workbook integrity mismatch: {key}")
        source_paths[key] = path
        commitments[key] = {"file": str(path), **observed}
    expected_reference = {
        "bytes": 835393456, "sha256": "1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0",
        "md5": "dcc3ea27079aa6dc3f9deccd7275e0f8"}
    if file_hashes(reference_path) != expected_reference:
        raise ValueError("independently pinned hg38 twoBit integrity mismatch")
    binding = protocol["schema_binding"]
    design_binding, score_binding = binding["design"], binding["scores"]
    expected_design_header = design_binding["all_header_columns"]
    expected_score_header = ["Fold", "Sequence ID"] + [score_binding["context_columns"][context] for context in contexts]
    designs, scores, eligible = {}, {}, {}
    exclusions = []
    chromosome_aliases = {}
    orientation_counts = defaultdict(int)
    source_rows = {"design": 0, "scores": 0, "blank_design_id": 0, "blank_score_id": 0}
    with IndependentWorkbook(source_paths["design_workbook"]) as book, TwoBitReader(reference_path) as genome:
        names = genome.sequence_sizes
        header_seen = False
        for number, cells in book.rows(design_binding["sheet"]):
            if number == design_binding["header_row_1based"]:
                if [cells.get(i) for i in range(len(expected_design_header))] != expected_design_header:
                    raise ValueError("S7 header differs from frozen binding")
                header_seen = True
            if number < design_binding["first_data_row_1based"]:
                continue
            source_rows["design"] += 1
            identifier = cells.get(0)
            if identifier is None or identifier == "":
                source_rows["blank_design_id"] += 1
                exclusions.append({"source": "S7", "source_row": number, "reason": "blank_design_id"})
                continue
            if not isinstance(identifier, str) or not identifier or "\n" in identifier:
                raise ValueError("invalid stable S7 ID")
            if identifier in designs:
                raise ValueError("duplicate S7 source ID")
            item = {"id": identifier, "source_row": number, "source_category": cells.get(1), "valid": False}
            designs[identifier] = item
            try:
                full = cells.get(10)
                if not isinstance(full, str) or len(full) != 230 or set(full) - set(BASES):
                    raise ValueError("invalid_230nt_assayed_oligo")
                start, end = _source_integer(cells.get(7)), _source_integer(cells.get(8))
                if start < 0 or end - start != 200:
                    raise ValueError("non200_or_invalid_genomic_interval")
                raw_chromosome = cells.get(6)
                if isinstance(raw_chromosome, ExcelNumber):
                    raw_chromosome = str(_source_integer(raw_chromosome))
                if not isinstance(raw_chromosome, str) or not raw_chromosome:
                    raise ValueError("missing_or_invalid_chromosome")
                candidates = [raw_chromosome] if raw_chromosome in names else ["chr" + raw_chromosome]
                if raw_chromosome == "MT":
                    candidates = ["chrM"]
                if len(candidates) != 1 or candidates[0] not in names:
                    raise ValueError("source_contig_absent_from_hg38")
                chromosome = candidates[0]
                chromosome_aliases[raw_chromosome] = chromosome
                strand = cells.get(9)
                if strand not in {"+", "-"}:
                    raise ValueError("missing_or_invalid_declared_strand")
                core = full[15:215]
                if genome.fetch(chromosome, start, end, strand=strand) != core:
                    raise ValueError("core_mismatch_at_declared_hg38_interval_and_strand")
                orientation_counts[strand] += 1
                item.update({"valid": True, "sequence": core, "chromosome": chromosome,
                             "start": start, "end": end, "strand": strand,
                             "parent_id": identifier,
                             "assayed_oligo_sha256": hashlib.sha256(full.encode()).hexdigest(),
                             "left_adapter_sha256": hashlib.sha256(full[:15].encode()).hexdigest(),
                             "right_adapter_sha256": hashlib.sha256(full[215:].encode()).hexdigest()})
            except (ValueError, KeyError) as error:
                item["reason"] = str(error)
                exclusions.append({"source": "S7", "source_row": number, "id": identifier, "reason": str(error)})
        if not header_seen:
            raise ValueError("S7 header row was absent")
    with IndependentWorkbook(source_paths["folds_workbook"]) as book:
        header_seen = False
        for number, cells in book.rows(score_binding["sheet"]):
            if number == score_binding["header_row_1based"]:
                if [cells.get(i) for i in range(len(expected_score_header))] != expected_score_header:
                    raise ValueError("S9 header differs from frozen binding")
                header_seen = True
            if number < score_binding["first_data_row_1based"]:
                continue
            source_rows["scores"] += 1
            identifier = cells.get(1)
            if identifier is None or identifier == "":
                source_rows["blank_score_id"] += 1
                exclusions.append({"source": "S9", "source_row": number, "reason": "blank_score_id"})
                continue
            if not isinstance(identifier, str) or not identifier or "\n" in identifier:
                raise ValueError("invalid stable S9 ID")
            if identifier in scores:
                raise ValueError("duplicate S9 source ID")
            values, literals, missing = {}, {}, []
            for index, context in enumerate(contexts, start=2):
                cell = cells.get(index)
                if not isinstance(cell, ExcelNumber):
                    missing.append(context)
                    continue
                try:
                    values[context] = cell.finite_float()
                    literals[context] = cell.literal
                except ValueError:
                    missing.append(context)
            scores[identifier] = {"id": identifier, "source_row": number, "observations": values,
                                  "source_literals": literals, "missing_contexts": missing}
            if missing:
                exclusions.append({"source": "S9", "source_row": number, "id": identifier,
                                   "reason": "missing_or_nonfinite_published_ESS", "contexts": missing})
        if not header_seen:
            raise ValueError("S9 header row was absent")
    for identifier, score in scores.items():
        if identifier not in designs:
            exclusions.append({"source": "join", "id": identifier, "reason": "S9_only_ID"})
        elif designs[identifier]["valid"] and not score["missing_contexts"]:
            design = designs[identifier]
            eligible[identifier] = {key: design[key] for key in ("id", "sequence", "chromosome", "start", "end", "parent_id")}
            eligible[identifier].update({"observations": score["observations"], "source_literals": score["source_literals"],
                                         "design_source_row": design["source_row"], "score_source_row": score["source_row"]})
    normalized = {row["id"]: row for row in normalized_records}
    if len(normalized) != len(normalized_records) or set(normalized) != set(eligible):
        raise ValueError("normalized panel eligibility differs from independent complete source reconstruction")
    for identifier, raw in eligible.items():
        for key in ("sequence", "chromosome", "start", "end", "parent_id", "observations"):
            if raw[key] != normalized[identifier][key]:
                raise ValueError(f"normalized panel differs from source at {identifier}/{key}")
    report = {"schema_version": 1, "source_workbooks": commitments, "reference": {"path": str(reference_path), **expected_reference},
              "source_rows": source_rows, "source_design_ids": len(designs), "source_score_ids": len(scores),
              "S7_only_IDs": sorted(set(designs) - set(scores)), "S9_only_IDs": sorted(set(scores) - set(designs)),
              "eligible_elements": len(eligible), "eligible_measurements": len(eligible) * len(contexts),
              "normalized_panel_exact_match": True, "chromosome_name_mapping": dict(sorted(chromosome_aliases.items())),
              "reference_valid_design_orientations": dict(orientation_counts), "exclusions": exclusions,
              "endpoint": ENDPOINT, "source_values_unchanged": True,
              "limits": ["Published ESS is whole-panel-normalized and includes test-element information",
                         "Neither raw activities nor source normalization population/SD are reconstructed",
                         "Reference interval matching does not establish genome-wide uniqueness"]}
    return [eligible[identifier] for identifier in sorted(eligible)], report


def checked_panel_bundle(directory: pathlib.Path, protocol: dict, protocol_sha: str) -> tuple[list[dict], dict]:
    root = directory.resolve()
    manifest_bytes = (root / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    if (manifest.get("protocol_sha256") != protocol_sha or manifest.get("endpoint") != ENDPOINT
            or manifest.get("contexts") != protocol["source"]["contexts"]):
        raise ValueError("panel manifest protocol/endpoint/context mismatch")
    required = {"panel.json", "source-schema-binding-v2.json", "reference-validation.json",
                "v2-ess-inspection-gate.json", "frozen-protocol-v2.json"}
    if not required <= set(manifest["files"]):
        raise ValueError("source manifest omits required semantic gate or data")
    for name, expected in manifest["files"].items():
        path = (root / name).resolve()
        if root not in path.parents or pathlib.Path(name).is_absolute() or ".." in pathlib.Path(name).parts:
            raise ValueError("source manifest path escapes bundle")
        if file_hashes(path)["sha256"] != expected:
            raise ValueError(f"source manifest checksum mismatch: {name}")
    for name in ("source-schema-binding-v2.json", "v2-ess-inspection-gate.json"):
        declaration = json.loads((root / name).read_bytes())
        if declaration.get("protocol_sha256") != protocol_sha:
            raise ValueError("source schema/inspection gate binds a different protocol")
    schema = json.loads((root / "source-schema-binding-v2.json").read_bytes())
    gate = json.loads((root / "v2-ess-inspection-gate.json").read_bytes())
    if not schema.get("semantic_endpoint_acceptance") or schema.get("endpoint") != ENDPOINT:
        raise ValueError("source semantic endpoint has not passed")
    if gate.get("frozen_before_ess_numeric_inspection") is not True:
        raise ValueError("source ESS inspection gate was not frozen first")
    if file_hashes(root / "frozen-protocol-v2.json")["sha256"] != protocol_sha:
        raise ValueError("source protocol copy differs from frozen version")
    panel_bytes = (root / "panel.json").read_bytes()
    panel = json.loads(panel_bytes)
    if panel.get("endpoint") != ENDPOINT or panel.get("contexts") != protocol["source"]["contexts"]:
        raise ValueError("source panel endpoint/context mismatch")
    if manifest.get("record_count") != len(panel["records"]):
        raise ValueError("panel manifest record count mismatch")
    return panel["records"], {"manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
                              "panel_sha256": hashlib.sha256(panel_bytes).hexdigest(),
                              "source_files": manifest["files"], "endpoint": ENDPOINT}


def compare_grouping(independent: dict, declared: dict) -> dict:
    if independent["records"] != declared["records"] or independent["components"] != declared["components"]:
        raise ValueError("independent sequence/identity/genomic components differ from declared groups")
    if len(independent["folds"]) != len(declared["folds"]):
        raise ValueError("independent fold count differs")
    for left, right in zip(independent["folds"], declared["folds"]):
        for key in ("fold", "test_groups", "guard_groups", "train_ids", "test_ids", "guard_ids"):
            if left[key] != right[key]:
                raise ValueError(f"independent fold {left['fold']} differs at {key}")
    if independent["candidate_pairs"] != declared["grouping_report"]["candidate_pairs"]:
        raise ValueError("independent unordered Hamming candidate count differs")
    return {"matched": True, "components": len(independent["components"]),
            "candidate_pairs": independent["candidate_pairs"], "folds": len(independent["folds"]),
            "relation_check": "independent bytewise exact Hamming and merged-expanded-interval guards"}


def write_json(path: pathlib.Path, value) -> None:
    content = (json.dumps(value, indent=2, allow_nan=False) + "\n").encode()
    if path.suffix == ".gz":
        with path.open("xb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream:
                stream.write(content)
    else:
        with path.open("xb") as stream:
            stream.write(content)


def run_replay(protocol_path: pathlib.Path, panel_directory: pathlib.Path, groups_path: pathlib.Path | None,
                reference_path: pathlib.Path, output: pathlib.Path, source_only: bool = False) -> None:
    if output.exists():
        raise FileExistsError("independent replay output already exists; retain each attempt")
    output.mkdir(parents=True)
    shutil.copyfile(__file__, output / "implementation.py")
    reference_reader = pathlib.Path(__file__).with_name("twobit_reference.py")
    shutil.copyfile(reference_reader, output / "twobit_reference.py")
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    try:
        protocol, protocol_sha = checked_protocol(protocol_path)
        normalized, commitments = checked_panel_bundle(panel_directory, protocol, protocol_sha)
        write_json(output / "inputs.json", {"protocol_sha256": protocol_sha, **commitments,
                   "python": platform.python_version(), "numpy": np.__version__,
                   "implementation_sha256": file_hashes(pathlib.Path(__file__))["sha256"], "started_at": started})
        write_json(output / "reference-reader.json", {"file": str(reference_reader),
                   "sha256": file_hashes(reference_reader)["sha256"], "snapshot_file": "twobit_reference.py"})
        print("Verifying every source design, declared reference interval and published ESS cell", flush=True)
        records, source_report = verify_raw_workbooks(panel_directory, protocol, normalized, reference_path)
        write_json(output / "source-verification.json", source_report)
        write_json(output / "source-records.json.gz", records)
        if not source_only:
            if groups_path is None:
                raise ValueError("full replay requires explicit frozen group ledger")
            groups_bytes = groups_path.read_bytes()
            declared = json.loads(groups_bytes)
            if declared.get("protocol_sha256") != protocol_sha or declared.get("panel_sha256") != commitments["panel_sha256"]:
                raise ValueError("group ledger does not bind this protocol and panel")
            write_json(output / "group-input.json", {"file": str(groups_path), "sha256": hashlib.sha256(groups_bytes).hexdigest()})
            print("Reconstructing all source-ID Hamming candidates, components, folds and interval guards", flush=True)
            grouping = protocol["grouping"]
            outer = protocol["outer_folds"]
            reconstructed = reconstruct_grouping(
                records, bin_bp=grouping["bin_bp"], candidate_budget=grouping["max_unordered_candidate_pairs"],
                fold_count=outer["count"], guard_bp=outer["guard_bp"],
                minimum_components=grouping["feasibility"]["minimum_components"],
                largest_component_fraction=grouping["feasibility"]["maximum_largest_component_fraction"],
                minimum_training_fraction=outer["feasibility"]["minimum_training_fraction_of_all_primary_elements"])
            write_json(output / "grouping.json.gz", reconstructed)
            write_json(output / "grouping-comparison.json", compare_grouping(reconstructed, declared))
            group_map = {row["id"]: row["group"] for row in reconstructed["records"]}
            for row in records:
                row["group"] = group_map[row["id"]]
            print("Fitting the frozen five models independently with training-only NumPy scaling", flush=True)
            evaluation = evaluate_complete_panel(records, protocol["source"]["contexts"], reconstructed["folds"],
                                                 protocol["fit"]["alpha"])
            evaluation["protocol_sha256"] = protocol_sha
            evaluation["input_commitments"] = commitments
            write_json(output / "evaluation.json.gz", evaluation)
            aliases = {"additive_context": "additive_context_dna", "independent_context": "independent_context_dna"}
            comparisons = [(aliases.get(row["baseline"], row["baseline"]),
                            aliases.get(row["candidate"], row["candidate"]))
                           for row in protocol["evaluation"]["comparisons"]]
            print("Replaying the common component bootstrap; fixed OOF losses, no refitting", flush=True)
            bootstrap = paired_component_bootstrap(
                [row["group"] for row in records],
                [row["model_losses"]["mse"] for row in evaluation["per_sequence_errors"]],
                [row["model_losses"]["mae"] for row in evaluation["per_sequence_errors"]], comparisons,
                protocol["resampling"]["replicates"], protocol["resampling"]["seed_u64"])
            write_json(output / "resampling.json.gz", bootstrap)
            write_json(output / "summary.json", {key: evaluation[key] for key in (
                "schema_version", "endpoint", "protocol_sha256", "contexts", "model_order", "sequence_count",
                "group_count", "measurement_count", "models", "balanced_design_checks", "fixed_zero_prediction_audit")})
        write_json(output / "status.json", {"status": "source_verified" if source_only else "complete",
                   "protocol_sha256": protocol_sha, "endpoint": ENDPOINT,
                   "completed_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   "source_elements": len(records), "source_measurements": len(records) * len(protocol["source"]["contexts"])})
    except BaseException as error:
        write_json(output / "failure.json", {"status": "failed_retained_attempt", "error_type": type(error).__name__,
                   "error": str(error), "traceback": traceback.format_exc(), "started_at": started,
                   "real_data_models_fitted_only_if_evaluation_exists": (output / "evaluation.json.gz").exists()})
        raise
    finally:
        files = {path.name: file_hashes(path)["sha256"] for path in sorted(output.iterdir())
                 if path.is_file() and path.name != "manifest.json"}
        write_json(output / "manifest.json", {"schema_version": 1, "endpoint": ENDPOINT,
                   "protocol_sha256": FROZEN_V2_SHA256, "files": files})


def read_json_payload(path: pathlib.Path) -> tuple[dict, dict]:
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        with gzip.open(path, "rb") as stream:
            decoded = stream.read(1024 * 1024 * 1024 + 1)
        if len(decoded) > 1024 * 1024 * 1024:
            raise ValueError("JSON payload exceeds 1GiB decompression budget")
    else:
        decoded = raw
    return json.loads(decoded), {"file": str(path), "raw_bytes": len(raw),
                                "raw_sha256": hashlib.sha256(raw).hexdigest(),
                                "decoded_bytes": len(decoded), "decoded_sha256": hashlib.sha256(decoded).hexdigest()}


def verified_output_manifest(root: pathlib.Path, name: str, schema_version: int) -> tuple[dict, dict]:
    root = root.resolve()
    path = root / name
    raw = path.read_bytes()
    manifest = json.loads(raw)
    if manifest["schema_version"] != schema_version:
        raise ValueError("unexpected artifact manifest schema")
    for filename, expected in manifest["files"].items():
        candidate = (root / filename).resolve()
        if root not in candidate.parents or pathlib.Path(filename).is_absolute() or ".." in pathlib.Path(filename).parts:
            raise ValueError("artifact path escapes bundle")
        if file_hashes(candidate)["sha256"] != expected:
            raise ValueError(f"artifact checksum mismatch: {filename}")
    return manifest, {"file": str(path), "sha256": hashlib.sha256(raw).hexdigest()}


class NumericalComparison:
    def __init__(self):
        self.checks = []
        self.failures = []

    def exact(self, category: str, left, right):
        matched = left == right
        self.checks.append({"category": category, "matched": matched, "kind": "exact"})
        if not matched:
            self.failures.append({"category": category, "kind": "exact", "detail": "values, membership or ordering differ"})

    def numeric(self, category: str, left, right, tolerance: float):
        a, b = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
        if a.shape != b.shape:
            self.failures.append({"category": category, "kind": "shape", "left_shape": list(a.shape), "right_shape": list(b.shape)})
            self.checks.append({"category": category, "matched": False, "kind": "numeric"})
            return
        finite = np.isfinite(a) & np.isfinite(b)
        delta = np.abs(a - b)
        mask = ~finite | (delta > tolerance)
        worst = float(np.max(delta[finite])) if np.any(finite) else (0.0 if not a.size else None)
        matched = not np.any(mask)
        self.checks.append({"category": category, "matched": matched, "kind": "numeric",
                            "values_checked": int(a.size), "absolute_tolerance": tolerance,
                            "maximum_absolute_difference": worst})
        if not matched:
            locations = np.argwhere(mask)
            self.failures.append({"category": category, "kind": "numeric", "mismatching_values": int(mask.sum()),
                                  "maximum_absolute_difference": worst, "absolute_tolerance": tolerance,
                                  "first_mismatch_indices": locations[:20].tolist()})

    def report(self):
        return {"matched": not self.failures, "failure_count": len(self.failures),
                "failures": self.failures, "checks": self.checks,
                "rank_metrics": "omitted by frozen protocol; no epsilon or rank acceptance rule"}


def compare_evaluations(independent: dict, rust: dict, comparison: NumericalComparison) -> None:
    p_tol, m_tol = 1e-7, 1e-8
    for key in ("contexts", "model_order", "sequence_count", "group_count", "measurement_count", "alpha"):
        comparison.exact("evaluation_" + key, independent[key], rust[key])
    comparison.exact("fold_count", len(independent["folds"]), len(rust["folds"]))
    for left, right in zip(independent["folds"], rust["folds"]):
        index = left["fold_index_0based"]
        comparison.exact(f"fold_{index}_index", index, right["fold_index_0based"])
        split = right["split"]
        for own_key, rust_key in (("train_ids", "training_sequence_ids"), ("test_ids", "test_sequence_ids"),
                                  ("guard_ids", "excluded_guard_sequence_ids")):
            comparison.exact(f"fold_{index}_{own_key}", left[own_key], split[rust_key])
        fit = right["fit"]
        comparison.exact(f"fold_{index}_feature_names", independent["feature_names"], fit["feature_names"])
        own_shared = left["models"]["shared_dna"]["shared"]
        comparison.numeric(f"fold_{index}_global_mean", own_shared["intercept"], fit["global_mean"], p_tol)
        own_context = left["models"]["context_mean"]["shared"]["context_intercepts"]
        comparison.numeric(f"fold_{index}_context_means", [own_context[c] for c in independent["contexts"]],
                           [fit["training_context_means"][c] for c in independent["contexts"]], p_tol)
        fits = [("shared", own_shared, fit["shared_dna_fit"])]
        for context in independent["contexts"]:
            fits.append((context, left["models"]["independent_context_dna"]["by_context"][context],
                         fit["independent_context_fits"][context]))
        for label, a, b in fits:
            comparison.exact(f"fold_{index}_{label}_training_rows", a["training_rows"], b["training_rows"])
            comparison.numeric(f"fold_{index}_{label}_alpha", a["alpha_summed_sse"], b["alpha"], 0.0)
            for key_a, key_b in (("feature_means", "training_feature_means"),
                                 ("feature_scales", "training_feature_scales"),
                                 ("standardized_coefficients", "standardized_coefficients"), ("intercept", "intercept")):
                comparison.numeric(f"fold_{index}_{label}_{key_a}", a[key_a], b[key_b], p_tol)
    def prediction_key(row):
        return row["sequence_id"], row["context"]
    left_predictions = sorted(independent["predictions"], key=prediction_key)
    right_predictions = sorted(rust["predictions"], key=prediction_key)
    for key in ("sequence_id", "context", "group", "fold_index_0based"):
        comparison.exact("prediction_" + key, [row[key] for row in left_predictions], [row[key] for row in right_predictions])
    comparison.numeric("observed_source_values", [row["observed"] for row in left_predictions],
                       [row["observed"] for row in right_predictions], 0.0)
    for model in MODELS:
        comparison.numeric("predictions_" + model,
                           [row["predictions"][model] for row in left_predictions],
                           [row["predictions"][model] for row in right_predictions], p_tol)
    for field, key in (("per_sequence_errors", "sequence_id"), ("per_group_errors", "group")):
        left = sorted(independent[field], key=lambda row: row[key])
        right = sorted(rust[field], key=lambda row: row[key])
        comparison.exact(field + "_membership", [row[key] for row in left], [row[key] for row in right])
        for loss in ("mse", "mae"):
            comparison.numeric(field + "_" + loss, [row["model_losses"][loss] for row in left],
                               [row["model_losses"][loss] for row in right], m_tol)
    left_models = {row["model"]: row for row in independent["models"]}
    right_models = {row["model"]: row for row in rust["models"]}
    comparison.exact("summary_model_inventory", sorted(left_models), sorted(right_models))
    for key in ("sequence_weighted_mse", "sequence_weighted_mae", "macro_group_mse", "macro_group_mae"):
        comparison.numeric("model_" + key, [left_models[model][key] for model in MODELS],
                           [right_models[model][key] for model in MODELS], m_tol)
    own_predictions = {(row["sequence_id"], row["context"]): row for row in independent["predictions"]}
    sequence_ids = sorted({row["sequence_id"] for row in independent["predictions"]})
    contexts = independent["contexts"]
    group_by_id = {row["sequence_id"]: row["group"] for row in independent["predictions"]}
    full_residual = np.array([
        [[own_predictions[(identifier, context)]["observed"] - own_predictions[(identifier, context)]["predictions"][model]
          for model in MODELS] for context in contexts] for identifier in sequence_ids])
    grouping = defaultdict(list)
    for index, identifier in enumerate(sequence_ids):
        grouping[group_by_id[identifier]].append(index)
    group_order = sorted(grouping)
    context_position = {context: index for index, context in enumerate(contexts)}
    for loss, values in (("mse", full_residual ** 2), ("mae", np.abs(full_residual))):
        right_sequences = sorted(rust["per_sequence_errors"], key=lambda row: row["sequence_id"])
        comparison.numeric("per_sequence_by_context_" + loss, values,
                           [[entry["model_losses"][loss] for entry in sorted(row["by_context"],
                             key=lambda item: context_position[item["context"]])] for row in right_sequences], m_tol)
        group_context = np.array([values[grouping[group]].mean(axis=0) for group in group_order])
        right_groups = sorted(rust["per_group_errors"], key=lambda row: row["group"])
        comparison.numeric("per_group_by_context_" + loss, group_context,
                           [[entry["model_losses"][loss] for entry in sorted(row["by_context"],
                             key=lambda item: context_position[item["context"]])] for row in right_groups], m_tol)
        expected_macro_context = group_context.mean(axis=0)
        comparison.numeric("model_macro_group_by_context_" + loss, expected_macro_context,
                           [[next(row for row in right_models[model]["macro_group_by_context"]
                                  if row["context"] == context)[loss] for model in MODELS] for context in contexts], m_tol)
    pairs = [(contexts[left], contexts[right]) for left in range(len(contexts)) for right in range(left + 1, len(contexts))]
    expected_keys = sorted((identifier, left, right) for identifier in sequence_ids for left, right in pairs)
    comparison.exact("pairwise_prediction_inventory", expected_keys,
                     sorted((row["sequence_id"], row["context_left"], row["context_right"])
                            for row in rust["context_contrast_predictions"]))
    expected_contrasts = []
    for row in rust["context_contrast_predictions"]:
        a = own_predictions[(row["sequence_id"], row["context_left"])]
        b = own_predictions[(row["sequence_id"], row["context_right"])]
        expected_contrasts.append([b["observed"] - a["observed"]]
                                 + [b["predictions"][model] - a["predictions"][model] for model in MODELS])
    comparison.exact("pairwise_contrast_count", len(expected_contrasts),
                     independent["sequence_count"] * len(independent["context_contrasts"]))
    comparison.numeric("pairwise_observations_and_predictions", expected_contrasts,
                       [[row["observed_difference"]] + [row["predicted_difference"][model] for model in MODELS]
                        for row in rust["context_contrast_predictions"]], p_tol)
    contrast_group_losses = {}
    for left, right in pairs:
        a, b = context_position[left], context_position[right]
        pair_residual = full_residual[:, b] - full_residual[:, a]
        for group in group_order:
            values = pair_residual[grouping[group]]
            contrast_group_losses[(group, left, right)] = {"mse": np.mean(values ** 2, axis=0),
                                                          "mae": np.mean(np.abs(values), axis=0)}
    contrast_keys = sorted(contrast_group_losses)
    root_pair_groups = {(row["group"], row["context_left"], row["context_right"]): row
                        for row in rust["context_contrast_group_errors"]}
    comparison.exact("pairwise_group_inventory", contrast_keys, sorted(root_pair_groups))
    if set(contrast_keys) == set(root_pair_groups):
        for loss in ("mse", "mae"):
            comparison.numeric("pairwise_group_" + loss, [contrast_group_losses[key][loss] for key in contrast_keys],
                               [root_pair_groups[key]["model_losses"][loss] for key in contrast_keys], m_tol)
    for key, value in independent["balanced_design_checks"].items():
        comparison.numeric("balanced_identity_" + key, value, rust["balanced_design_checks"][key], m_tol)


def compare_bootstraps(independent: dict, rust: dict, comparison: NumericalComparison,
                       own_models: list[dict]) -> None:
    for key in ("group_order", "model_order", "comparison_order", "replicates", "seed", "interval_probabilities", "resampling_indices"):
        comparison.exact("bootstrap_" + key, independent[key], rust[key])
    comparison.exact("bootstrap_primary_weighting", rust["weighting"], "equal_sequences")
    primary = independent["equal_sequences"]
    for own_key, rust_key in (("sampled_mse_improvements", "sampled_mse_improvements"),
                              ("sampled_mae_improvements", "sampled_mae_improvements"),
                              ("sampled_model_mses", "sampled_model_mse"),
                              ("sampled_model_maes", "sampled_model_mae")):
        comparison.numeric("bootstrap_" + rust_key, primary[own_key], rust[rust_key], 1e-8)
    for index, (left, right) in enumerate(zip(primary["comparisons"], rust["comparisons"])):
        comparison.exact(f"bootstrap_comparison_{index}_identity", left["comparison"], right["comparison"])
        for key in ("observed_mse_improvement", "observed_mae_improvement", "mse_interval", "mae_interval"):
            comparison.numeric(f"bootstrap_comparison_{index}_{key}", left[key], right[key], 1e-8)
    own_model_lookup = {row["model"]: row for row in own_models}
    root_model_lookup = {row["model"]: row for row in rust["models"]}
    comparison.exact("bootstrap_model_inventory", sorted(own_model_lookup), sorted(root_model_lookup))
    for loss, singular, plural in (("mse", "sampled_model_mse", "sampled_model_mses"),
                                   ("mae", "sampled_model_mae", "sampled_model_maes")):
        observed = [own_model_lookup[model]["sequence_weighted_" + loss] for model in MODELS]
        comparison.numeric("bootstrap_observed_model_" + loss, observed, rust["observed_model_" + loss], 1e-8)
        intervals = np.quantile(np.array(primary[plural]), [0.025, 0.975], axis=0, method="linear").T
        comparison.numeric("bootstrap_all_model_" + loss + "_intervals", intervals,
                           [root_model_lookup[model][loss + "_interval"] for model in MODELS], 1e-8)
        comparison.numeric("bootstrap_all_model_observed_" + loss, observed,
                           [root_model_lookup[model]["observed_" + loss] for model in MODELS], 1e-8)


def compare_sequence_features(independent: dict, source_records: list[dict], rust_features: dict,
                              rust_study: dict, comparison: NumericalComparison) -> None:
    comparison.exact("feature_vocabulary", list(FEATURE_NAMES), rust_features["feature_names"])
    comparison.exact("feature_panel_digest", rust_study["panel_metadata_sha256"], rust_features["panel_metadata_sha256"])
    expected = [{"id": row["id"], **integer_sequence_features(row["sequence"])}
                for row in sorted(source_records, key=lambda row: row["id"])]
    actual = sorted(rust_features["records"], key=lambda row: row["id"])
    comparison.exact("feature_ID_inventory", [row["id"] for row in expected], [row["id"] for row in actual])
    for key in ("length", "gc_count", "canonical_trimer_counts"):
        comparison.exact("feature_integer_" + key, [row[key] for row in expected], [row[key] for row in actual])
    own_vectors = [[count / (row["length"] - 2) for count in row["canonical_trimer_counts"]]
                   + [row["gc_count"] / row["length"]] for row in expected]
    comparison.exact("saved_feature_matrix_ID_order", [row["id"] for row in expected], independent["sequence_ids"])
    comparison.numeric("independent_saved_feature_fractions", own_vectors, independent["features"], 0.0)
    comparison.numeric("feature_exact_float_fractions", own_vectors, [row["vector"] for row in actual], 0.0)


def compare_rust_run(independent_root: pathlib.Path, rust_study_path: pathlib.Path, output: pathlib.Path) -> bool:
    if output.exists():
        raise FileExistsError("comparison output already exists; retain each attempt")
    output.mkdir(parents=True)
    shutil.copyfile(__file__, output / "implementation.py")
    try:
        own_manifest, own_commitment = verified_output_manifest(independent_root, "manifest.json", 1)
        root = rust_study_path.resolve().parent
        rust_manifest, rust_commitment = verified_output_manifest(root, "artifacts.json", 4)
        expected_files = {"features.json.gz", "findings.md", "grouping.json", "inputs.json", "ledger.json", "protocol.json",
                          "reference.json", "resampling.json.gz", "source-manifest.json", "study.json.gz", "summary.json"}
        if set(rust_manifest["files"]) != expected_files or rust_study_path.resolve() != (root / "study.json.gz"):
            raise ValueError("Rust core artifact inventory or compared study path differs")
        own_inputs = json.loads((independent_root / "inputs.json").read_bytes())
        rust_inputs = json.loads((root / "inputs.json").read_bytes())
        own_group = json.loads((independent_root / "group-input.json").read_bytes())
        for own_value, rust_value, name in (
            (own_inputs["protocol_sha256"], rust_inputs["protocol_sha256"], "protocol"),
            (own_inputs["manifest_sha256"], rust_inputs["source_manifest_sha256"], "source manifest"),
            (own_inputs["panel_sha256"], rust_inputs["panel_sha256"], "panel"),
            (own_group["sha256"], rust_inputs["grouping_sha256"], "group ledger")):
            if own_value != rust_value:
                raise ValueError(f"independent and Rust input commitments differ: {name}")
        if own_inputs["protocol_sha256"] != FROZEN_V2_SHA256 or rust_inputs["endpoint"] != ENDPOINT:
            raise ValueError("comparison protocol/endpoint mismatch")
        own, own_study_commitment = read_json_payload(independent_root / "evaluation.json.gz")
        rust, rust_study_commitment = read_json_payload(rust_study_path)
        own_boot, own_boot_commitment = read_json_payload(independent_root / "resampling.json.gz")
        rust_boot, rust_boot_commitment = read_json_payload(root / "resampling.json.gz")
        rust_features, feature_commitment = read_json_payload(root / "features.json.gz")
        raw_records, source_commitment = read_json_payload(independent_root / "source-records.json.gz")
        comparison = NumericalComparison()
        compare_sequence_features(own, raw_records, rust_features, rust, comparison)
        compare_evaluations(own, rust, comparison)
        compare_bootstraps(own_boot, rust_boot, comparison, own["models"])
        rust_summary = json.loads((root / "summary.json").read_bytes())
        for left, right in (("sequence_weighted_mse", "mse"), ("sequence_weighted_mae", "mae")):
            comparison.numeric("fixed_zero_" + right, own["fixed_zero_prediction_audit"][left],
                               rust_summary["fixed_zero_prediction_diagnostic"][right], 1e-8)
        result = comparison.report()
        result.update({"schema_version": 1, "endpoint": ENDPOINT, "protocol_sha256": FROZEN_V2_SHA256,
                       "independent_manifest": own_commitment, "rust_artifacts_manifest": rust_commitment,
                       "independent_study": own_study_commitment, "rust_study": rust_study_commitment,
                       "independent_resampling": own_boot_commitment, "rust_resampling": rust_boot_commitment,
                       "rust_features": feature_commitment, "independent_raw_source_records": source_commitment,
                       "implementation_sha256": file_hashes(pathlib.Path(__file__))["sha256"],
                       "prediction_absolute_tolerance": 1e-7, "metric_absolute_tolerance": 1e-8,
                       "limits": ["Numerical implementation agreement does not resolve the source ESS formula discrepancy",
                                  "Source whole-panel target normalization is not prospective or training-only",
                                  "Raw all-row source verification and independent grouping are bound by the independent manifest"]})
        write_json(output / "comparison.json", result)
        print(json.dumps({"matched": result["matched"], "failure_count": result["failure_count"],
                          "comparison": str(output / "comparison.json")}))
        return result["matched"]
    except BaseException as error:
        write_json(output / "failure.json", {"status": "failed_retained_comparison", "error_type": type(error).__name__,
                   "error": str(error), "traceback": traceback.format_exc()})
        raise
    finally:
        write_json(output / "manifest.json", {"schema_version": 1, "protocol_sha256": FROZEN_V2_SHA256,
                   "files": {path.name: file_hashes(path)["sha256"] for path in sorted(output.iterdir())
                             if path.is_file() and path.name != "manifest.json"}})


def published_ess_semantic_audit(records: Sequence[dict], contexts: Sequence[str]) -> dict:
    """Prespecified zero-sum/source-precision audit, never a target correction.

    Source XML numeric literals are parsed as exact rational numbers; no decimal
    context precision, floating summation or selected epsilon can erase them.
    All three fixed protocol thresholds are retained, regardless of result.
    """
    rows = sorted(records, key=lambda row: row["id"])
    if not rows:
        raise ValueError("ESS audit needs nonempty raw verified records")
    thresholds = [Fraction("1e-12"), Fraction("1e-8"), Fraction("1e-6")]
    counts = [0, 0, 0]
    signed = []
    worst = None
    largest_magnitude = Fraction(-1)
    for row in rows:
        literals = row["source_literals"]
        exact = sum((Fraction(literals[context]) for context in contexts), Fraction(0))
        magnitude = abs(exact)
        signed.append(float(exact))
        for index, threshold in enumerate(thresholds):
            counts[index] += magnitude > threshold
        if magnitude > largest_magnitude:
            largest_magnitude = magnitude
            worst = {"id": row["id"], "score_source_row_1based": row["score_source_row"],
                     "source_literals": literals, "sum_exact_numerator": str(exact.numerator),
                     "sum_exact_denominator": str(exact.denominator), "sum_float": float(exact)}
    sums = np.array(signed)
    values = _matrix([[row["observations"][context] for context in contexts] for row in rows], "ESS audit values")
    return {
        "schema_version": 1, "endpoint": ENDPOINT, "elements": len(rows), "contexts": list(contexts),
        "measurement_count": len(rows) * len(contexts),
        "source_cell_mapping": {context: context + ", Element Specificity Score" for context in contexts},
        "source_numeric_arithmetic": "exact fractions parsed independently from raw XLSX XML numeric literals",
        "mean_signed_sum": float(sums.mean()), "mean_absolute_sum": float(np.abs(sums).mean()),
        "root_mean_squared_sum": float(np.sqrt(np.mean(sums ** 2))),
        "maximum_absolute_sum": float(largest_magnitude),
        "count_abs_sum_gt1e-12": counts[0], "count_abs_sum_gt1e-8": counts[1], "count_abs_sum_gt1e-6": counts[2],
        "maximum_source_row": worst,
        "fixed_zero": {"mse": float(np.mean(values ** 2)), "mae": float(np.mean(np.abs(values))),
                       "per_context": [{"context": context, "mse": float(np.mean(values[:, index] ** 2)),
                                        "mae": float(np.mean(np.abs(values[:, index])))}
                                       for index, context in enumerate(contexts)]},
        "simple_element_centering_formula_audit": "discrepancy" if counts[2] else "no_large_residual_at_prespecified_thresholds",
        "policy": "All source values retained unchanged; no recentering, snapping, exclusion, new feature, model or refit",
        "interpretation": "The exact S9 headers identify these numbers as published ESS. Residuals test applicability of the cited simple final-paper formula to this preprint-associated source release; they do not reconstruct raw activity or identify an alternative normalization.",
    }


def run_endpoint_audit(independent_root: pathlib.Path, output: pathlib.Path) -> None:
    if output.exists():
        raise FileExistsError("endpoint-audit output already exists; retain attempts")
    _, commitment = verified_output_manifest(independent_root, "manifest.json", 1)
    records, raw_commitment = read_json_payload(independent_root / "source-records.json.gz")
    inputs = json.loads((independent_root / "inputs.json").read_bytes())
    if inputs["protocol_sha256"] != FROZEN_V2_SHA256:
        raise ValueError("endpoint audit binds another protocol")
    audit = published_ess_semantic_audit(records, ["HepG2", "K562", "WTC11"])
    audit["protocol_sha256"] = FROZEN_V2_SHA256
    audit["input_commitments"] = {"independent_manifest": commitment, "source_records": raw_commitment,
                                  "source_manifest_sha256": inputs["manifest_sha256"], "panel_sha256": inputs["panel_sha256"]}
    output.mkdir(parents=True)
    shutil.copyfile(__file__, output / "implementation.py")
    write_json(output / "endpoint-audit.json", audit)
    write_json(output / "manifest.json", {"schema_version": 1, "protocol_sha256": FROZEN_V2_SHA256,
               "files": {path.name: file_hashes(path)["sha256"] for path in sorted(output.iterdir()) if path.is_file()}})
    print(json.dumps({"endpoint_audit": str(output / "endpoint-audit.json"),
                      "maximum_absolute_sum": audit["maximum_absolute_sum"],
                      "count_abs_sum_gt1e-6": audit["count_abs_sum_gt1e-6"]}))


def describe() -> dict:
    return {
        "status": "frozen_context_specificity_v2_explicit_replay",
        "real_outcome_loading_enabled": "explicit --out only; frozen protocol and source semantic/hash gates enforced",
        "protocol_sha256": FROZEN_V2_SHA256, "endpoint": ENDPOINT,
        "models": list(MODELS), "features": list(FEATURE_NAMES),
        "python": platform.python_version(), "numpy": np.__version__,
        "implementation_sha256": hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
        "ridge_objective": "summed_row_SSE + explicit_alpha * squared_standardized_slopes",
        "context_intercepts": "unpenalized; additive residualizes both X and y within context",
        "source_loading": "independent all-row raw XLSX verification, unchanged published ESS; not raw activities",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--describe", action="store_true", help="print helper status and runtime versions")
    parser.add_argument("--protocol", type=pathlib.Path, default=pathlib.Path("experiments/context-specificity-v2.json"))
    parser.add_argument("--dataset", type=pathlib.Path, default=pathlib.Path("datasets/context-panel"))
    parser.add_argument("--groups", type=pathlib.Path)
    parser.add_argument("--reference", type=pathlib.Path, default=pathlib.Path("/workspace/.tools/reference/hg38.2bit"))
    parser.add_argument("--out", type=pathlib.Path)
    parser.add_argument("--source-check-only", action="store_true")
    parser.add_argument("--compare-rust", type=pathlib.Path)
    parser.add_argument("--independent-run", type=pathlib.Path)
    parser.add_argument("--audit-endpoint", action="store_true")
    args = parser.parse_args()
    if args.audit_endpoint:
        if args.out is None or args.independent_run is None:
            parser.error("--audit-endpoint requires --independent-run and new --out")
        run_endpoint_audit(args.independent_run, args.out)
    elif args.compare_rust is not None:
        if args.out is None or args.independent_run is None:
            parser.error("--compare-rust requires --independent-run and new --out")
        if not compare_rust_run(args.independent_run, args.compare_rust, args.out):
            raise SystemExit(2)
    elif args.out is None or args.describe:
        print(json.dumps(describe(), indent=2, allow_nan=False))
    else:
        run_replay(args.protocol, args.dataset, args.groups, args.reference, args.out, args.source_check_only)


if __name__ == "__main__":
    main()
