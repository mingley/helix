#!/usr/bin/env python3
"""Independent NumPy/SciPy replication of the frozen regulatory-prediction study.

This implementation does not call Rust feature, fold, or fitting code. Real
outcome parsing/fitting is gated on a frozen protocol and verified input hashes.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import pathlib
import re
from dataclasses import dataclass
import itertools
import sys
import platform
import argparse
import datetime
from fractions import Fraction
import gzip
from scipy.linalg import cho_factor, cho_solve
import scipy

import numpy as np
from scipy.stats import rankdata

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def blocked_folds(
    positions: np.ndarray,
    blocks: list[tuple[int, int]],
    guard: int,
) -> list[dict]:
    """Construct explicit inclusive-position blocks and expanded exclusion zones."""
    p = np.asarray(positions, dtype=np.int64)
    if p.ndim != 1 or len(p) == 0 or guard < 0:
        raise ValueError("fold positions must be nonempty and guard nonnegative")
    if any(start > end for start, end in blocks):
        raise ValueError("fold block is reversed")
    coverage = np.zeros(len(p), dtype=np.int64)
    folds = []
    for fold, (start, end) in enumerate(blocks):
        test = (p >= start) & (p <= end)
        if not test.any():
            raise ValueError("every fold must have test observations")
        # The protocol guards actual retained test centers, not unobserved sites.
        distance = np.abs(p[:, None] - np.unique(p[test])[None, :]).min(axis=1)
        train = distance > guard
        excluded = ~(test | train)
        coverage += test.astype(np.int64)
        if not test.any() or not train.any():
            raise ValueError("every fold must have training and test observations")
        folds.append({
            "fold": fold,
            "start": start,
            "end": end,
            "train_indices": np.flatnonzero(train).tolist(),
            "test_indices": np.flatnonzero(test).tolist(),
            "guard_indices": np.flatnonzero(excluded).tolist(),
        })
    if not np.all(coverage == 1):
        raise ValueError("each observation must belong to exactly one test block")
    return folds


def motif_delta(
    sequence: str,
    position_0based: int,
    alternate: str,
    counts: dict[str, list[float]],
) -> dict:
    """A fixed PFM's max score change over identical overlapping window sets."""
    if set(sequence) - set(BASES) or alternate not in BASES:
        raise ValueError("motif scoring requires unambiguous DNA")
    if not 0 <= position_0based < len(sequence):
        raise ValueError("motif variant position is outside target")
    matrix = np.array([counts[base] for base in BASES], dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != 4 or matrix.shape[1] == 0:
        raise ValueError("PFM must have four nonempty aligned base rows")
    if not np.isfinite(matrix).all() or np.any(matrix < 0):
        raise ValueError("PFM counts must be finite and nonnegative")
    width = matrix.shape[1]
    probability = (matrix + 0.25) / (matrix.sum(axis=0) + 1.0)
    weight = np.log2(probability / 0.25)
    starts = range(max(0, position_0based - width + 1),
                   min(position_0based, len(sequence) - width) + 1)
    reference_scores = []
    alternate_scores = []

    def score(window):
        indices = np.fromiter((BASES.index(base) for base in window), dtype=np.int64)
        return float(weight[indices, np.arange(width)].sum())

    for start in starts:
        reference = sequence[start:start + width]
        index = position_0based - start
        changed = reference[:index] + alternate + reference[index + 1:]
        for left, right in [(reference, changed),
                            (reference.translate(COMPLEMENT)[::-1],
                             changed.translate(COMPLEMENT)[::-1])]:
            reference_scores.append(score(left))
            alternate_scores.append(score(right))
    if not reference_scores:
        raise ValueError("no complete overlapping motif window")
    before = max(reference_scores)
    after = max(alternate_scores)
    return {
        "reference_max": before,
        "alternate_max": after,
        "delta": after - before,
        "coordinate_orientation_candidates": len(reference_scores),
    }


def checked_motifs(directory: pathlib.Path, ids: list[str]) -> tuple[list[dict], dict]:
    """Load prior-frozen PFMs, checking both raw and normalized representations."""
    root = directory.resolve()
    raw_manifest = (root / "manifest.json").read_bytes()
    manifest = json.loads(raw_manifest)
    contents = {}
    for name, expected in manifest["files"].items():
        path = (root / name).resolve()
        if pathlib.Path(name).name != name or path.parent != root:
            raise ValueError("motif input path escapes bundle")
        raw = path.read_bytes()
        if sha256(raw) != expected:
            raise ValueError(f"motif checksum mismatch: {name}")
        contents[name] = raw
    parsed = json.loads(contents[manifest["motifs_file"]])
    if parsed["schema_version"] != 1:
        raise ValueError("unsupported motif schema")
    profiles = {profile["matrix_id"]: profile for profile in parsed["motifs"]}
    if len(profiles) != len(parsed["motifs"]) or set(profiles) != set(ids):
        raise ValueError("motif inventory disagrees with frozen protocol")
    result = []
    for matrix_id in ids:
        profile = profiles[matrix_id]
        raw = json.loads(contents[profile["source_file"]])
        if raw["matrix_id"] != matrix_id or raw["pfm"] != profile["pfm"]:
            raise ValueError("normalized motif disagrees with source PFM")
        if any(len(profile["pfm"][base]) != profile["width"] for base in BASES):
            raise ValueError("motif width disagrees with PFM")
        result.append(profile)
    return result, {
        "directory": str(directory),
        "manifest_sha256": sha256(raw_manifest),
        "files": manifest["files"],
        "motif_ids": ids,
    }


def ridge_predict(
    train_x: np.ndarray,
    train_y: np.ndarray,
    test_x: np.ndarray,
    alpha: float,
) -> tuple[np.ndarray, dict]:
    """Ridge with training-only population-SD scaling and a free intercept.

    Objective: sum((y - intercept - X_scaled beta)^2) + alpha * ||beta||².
    There is no division by training sample count in the ridge penalty.
    """
    x = np.asarray(train_x, dtype=np.float64)
    y = np.asarray(train_y, dtype=np.float64)
    q = np.asarray(test_x, dtype=np.float64)
    if x.ndim != 2 or q.ndim != 2 or x.shape[1] != q.shape[1]:
        raise ValueError("training/test features must have matching column dimensions")
    if y.ndim != 1 or len(y) != len(x) or len(y) == 0:
        raise ValueError("training targets must be a nonempty matching vector")
    if not np.isfinite(x).all() or not np.isfinite(y).all() or not np.isfinite(q).all():
        raise ValueError("ridge inputs must be finite")
    if not math.isfinite(alpha) or alpha <= 0:
        raise ValueError("ridge alpha must be finite and positive")
    mean = x.mean(axis=0)
    scale = x.std(axis=0, ddof=0)
    scale[scale <= 1e-12] = 1.0
    z = (x - mean) / scale
    center_y = float(y.mean())
    if x.shape[1]:
        gram = z.T @ z + alpha * np.eye(x.shape[1])
        beta = np.linalg.solve(gram, z.T @ (y - center_y))
        prediction = center_y + ((q - mean) / scale) @ beta
    else:
        beta = np.zeros(0, dtype=np.float64)
        prediction = np.full(len(q), center_y, dtype=np.float64)
    return prediction, {
        "training_feature_mean": mean.tolist(),
        "training_feature_scale": scale.tolist(),
        "training_target_mean": center_y,
        "coefficients_scaled": beta.tolist(),
    }


def correlation(left: np.ndarray, right: np.ndarray) -> float | None:
    a = np.asarray(left, dtype=np.float64)
    b = np.asarray(right, dtype=np.float64)
    if len(a) < 2:
        return None
    if np.all(a == a[0]) or np.all(b == b[0]):
        return None
    a = a - a.mean()
    b = b - b.mean()
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denominator == 0.0:
        return None
    return float(np.dot(a, b) / denominator)


def metrics(actual: np.ndarray, predicted: np.ndarray) -> dict:
    y = np.asarray(actual, dtype=np.float64)
    p = np.asarray(predicted, dtype=np.float64)
    if y.ndim != 1 or p.shape != y.shape or len(y) == 0:
        raise ValueError("metrics require matching nonempty vectors")
    if not np.isfinite(y).all() or not np.isfinite(p).all():
        raise ValueError("metrics require finite values")
    error = y - p
    squared = float(np.dot(error, error))
    target_variance = 0.0 if np.all(y == y[0]) else float(np.dot(y - y.mean(), y - y.mean()))
    return {
        "count": len(y),
        "mse": squared / len(y),
        "rmse": math.sqrt(squared / len(y)),
        "mae": float(np.abs(error).mean()),
        "r_squared": None if target_variance == 0.0 else 1.0 - squared / target_variance,
        "pearson": correlation(y, p),
        "spearman": correlation(rankdata(y, method="average"), rankdata(p, method="average")),
    }


@dataclass(frozen=True)
class Observation:
    assay_id: str
    accession: str
    position: int
    reference: str
    alternate: str
    score: float
    hgvs: str
    genomic_key: str = ""
    score_raw: str = ""

    @property
    def key(self):
        return self.assay_id + "|" + self.hgvs


def parse_rows(raw: bytes, target: str, assay_id: str, expected_rows: int,
               mapping: dict | None = None) -> tuple[list[Observation], dict]:
    """Apply only the frozen representation and input-quality rules."""
    observed = []
    excluded = []
    seen = set()
    accessions = set()
    count = 0
    for row in csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))):
        count += 1
        accession = row["accession"]
        if not accession or accession in accessions:
            raise ValueError("missing or duplicate source accession")
        accessions.add(accession)
        hgvs = row["hgvs_nt"]
        match = re.fullmatch(r"n\.(\d+)([ACGT])>([ACGT])", hgvs)
        if not match:
            equality = re.fullmatch(r"n\.(?:\d+[ACGT]?)?=", hgvs)
            excluded.append({"accession": accession, "hgvs": hgvs,
                             "reasons": ["equality" if equality else "non_snv"]})
            continue
        position, reference, alternate = int(match[1]), match[2], match[3]
        if not 1 <= position <= len(target) or target[position - 1] != reference:
            raise ValueError("SNV reference/coordinate disagrees with target")
        if reference == alternate:
            raise ValueError("SNV requires differing alleles")
        key = (position, reference, alternate)
        if key in seen:
            raise ValueError("duplicate single-nucleotide substitution")
        seen.add(key)

        def number(name):
            try:
                value = float(row[name])
            except (TypeError, ValueError):
                return None
            return value if math.isfinite(value) else None

        score = number("score")
        tags = number("unique_tags")
        dna = number("dna_sequences")
        rna = number("rna_sequences")
        reasons = []
        if score is None:
            reasons.append("nonfinite_or_missing_score")
        if tags is None or tags < 10:
            reasons.append("unique_tags_below_10_or_invalid")
        if dna is None or dna <= 0:
            reasons.append("dna_count_nonpositive_or_invalid")
        if rna is None or rna <= 0:
            reasons.append("rna_count_nonpositive_or_invalid")
        if reasons:
            excluded.append({"accession": accession, "hgvs": hgvs, "reasons": reasons})
            continue
        genomic_key = ""
        if mapping is not None:
            strand = mapping.get("target_reference_alignment_strand",
                                 mapping.get("target_strand"))
            start = mapping.get("start0", mapping.get("start", mapping.get("target_start")))
            end = mapping.get("end0", mapping.get("end", mapping.get("target_end")))
            chromosome = mapping["chromosome"].removeprefix("chr")
            if strand == "+":
                offset = start + position - 1
                ref, alt = reference, alternate
            elif strand == "-":
                offset = end - position
                ref, alt = reference.translate(COMPLEMENT), alternate.translate(COMPLEMENT)
            else:
                raise ValueError("missing reference alignment orientation")
            genomic_key = f"GRCh38:{chromosome}:{offset}:{ref}>{alt}"
        observed.append(Observation(
            assay_id=assay_id, accession=accession, position=position,
            reference=reference, alternate=alternate, score=score, hgvs=hgvs,
            genomic_key=genomic_key,
            score_raw=row["score"],
        ))
    if count != expected_rows:
        raise ValueError("source row count disagrees with pinned metadata")
    observed.sort(key=lambda row: (row.assay_id, row.position, row.alternate))
    allele_counts = {}
    for row in observed:
        allele_counts[row.position] = allele_counts.get(row.position, 0) + 1
    return observed, {
        "source_rows": count, "included_rows": len(observed),
        "excluded_rows": excluded, "position_allele_counts": allele_counts,
        "incomplete_positions": [p for p, count in sorted(allele_counts.items()) if count != 3],
    }


def feature_matrix(target: str, rows: list[Observation], profiles: list[dict],
                   protocol: dict) -> tuple[np.ndarray, dict[str, list[int]], list[str]]:
    """Reconstruct all 82 frozen features directly from reference strings."""
    mutation_order = protocol["features"]["mutation_identity_order"]
    kmer_order = protocol["features"]["delta3mer_order"]
    if kmer_order != ["".join(s) for s in itertools.product(BASES, repeat=3)]:
        raise ValueError("unexpected frozen 3mer vocabulary")
    mutation_index = {value: i for i, value in enumerate(mutation_order)}
    kmer_index = {value: i for i, value in enumerate(kmer_order)}
    columns = (
        ["mutation:" + value for value in mutation_order]
        + ["position:x", "position:x_squared"]
        + ["delta3mer:" + value for value in kmer_order]
        + ["motif_delta:" + profile["matrix_id"] for profile in profiles]
    )
    result = np.zeros((len(rows), len(columns)), dtype=np.float64)
    if len(target) < 2:
        raise ValueError("target must contain at least two bases")
    for i, row in enumerate(rows):
        position = row.position - 1
        if row.reference != target[position]:
            raise ValueError("feature row reference disagrees with target")
        identity = row.reference + ">" + row.alternate
        result[i, mutation_index[identity]] = 1.0
        x = 2.0 * position / (len(target) - 1) - 1.0
        result[i, 12:14] = [x, x * x]
        for start in range(max(0, position - 2), min(position, len(target) - 3) + 1):
            ref = target[start:start + 3]
            local_index = position - start
            alt = ref[:local_index] + row.alternate + ref[local_index + 1:]
            result[i, 14 + kmer_index[ref]] -= 1.0
            result[i, 14 + kmer_index[alt]] += 1.0
        for j, profile in enumerate(profiles):
            result[i, 78 + j] = motif_delta(
                target, position, row.alternate, profile["pfm"])["delta"]
    groups = {
        "mutation_identity": list(range(12)),
        "position": list(range(12, 14)),
        "delta3mer": list(range(14, 78)),
        "motif_deltas": list(range(78, 82)),
    }
    models = {
        model["id"]: [column for group in model["feature_sets"] for column in groups[group]]
        for model in protocol["models"]
    }
    return result, models, columns


def checked_assay(directory: pathlib.Path, protocol_frozen: bool) -> tuple[str, list[Observation], dict]:
    """Validate raw commitments before parsing any real outcome values."""
    if not protocol_frozen:
        raise ValueError("real outcomes cannot be parsed before protocol freeze")
    root = directory.resolve()
    manifest_bytes = (root / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    contents = {}
    for name, expected in manifest["files"].items():
        if pathlib.Path(name).name != name:
            raise ValueError("assay paths must be basenames")
        path = (root / name).resolve()
        if path.parent != root:
            raise ValueError("assay path escapes dataset")
        raw = path.read_bytes()
        if sha256(raw) != expected:
            raise ValueError(f"input checksum mismatch: {name}")
        contents[name] = raw
    metadata = json.loads(contents["metadata.json"])
    target = metadata["targetGenes"][0]["targetSequence"]["sequence"].upper()
    if set(target) - set("ACGT"):
        raise ValueError("target sequence must be unambiguous DNA")
    fasta = contents["target.fasta"].decode("ascii").splitlines()
    if not fasta or not fasta[0].startswith(">") or "".join(fasta[1:]).upper() != target:
        raise ValueError("FASTA target does not match source metadata")
    if manifest["expected_rows"] != metadata["numVariants"]:
        raise ValueError("source row count disagrees with pinned metadata")
    observed, accounting = parse_rows(
        contents["scores.csv"], target, manifest["assay_urn"], manifest["expected_rows"], manifest)
    return target, observed, {
        "directory": str(directory),
        "manifest_sha256": sha256(manifest_bytes),
        "files": manifest["files"],
        "assay_id": manifest["assay_urn"],
        "target_length": len(target),
        "row_accounting": accounting,
    }


def checked_cohort(directory: pathlib.Path, protocol: dict, protocol_hash: str):
    root = directory.resolve()
    manifest_raw = (root / "manifest.json").read_bytes()
    manifest = json.loads(manifest_raw)
    contents = {}
    for name, expected in manifest["files"].items():
        path = (root / name).resolve()
        if pathlib.Path(name).name != name or path.parent != root:
            raise ValueError("cohort path escapes input bundle")
        raw = path.read_bytes()
        if sha256(raw) != expected:
            raise ValueError(f"cohort checksum mismatch: {name}")
        contents[name] = raw
    if manifest["frozen_protocol_sha256"] != protocol_hash:
        raise ValueError("cohort acquisition protocol differs from requested protocol")
    if sha256(contents["frozen-protocol.json"]) != protocol_hash:
        raise ValueError("cohort does not contain exact frozen acquisition protocol")
    gate = json.loads(contents["outcome-acquisition-gate.json"])
    if gate["protocol_sha256"] != protocol_hash or not gate["frozen_before_new_outcome_acquisition"]:
        raise ValueError("cohort acquisition gate is not bound to this protocol")
    roles = {item["assay_id"]: item["role"] for item in protocol["cohort"]["assay_roles"]}
    if gate["assay_roles"] != roles:
        raise ValueError("cohort acquisition roles differ from protocol")
    assays = {}
    accounting = {}
    for assay in manifest["assays"]:
        assay_id = assay["assay_id"]
        if assay_id in assays or roles.get(assay_id) != assay["role"]:
            raise ValueError("cohort assay identity/role mismatch")
        metadata = json.loads(contents[assay["metadata_file"]])
        target = metadata["targetGenes"][0]["targetSequence"]["sequence"].upper()
        if target != assay["target_sequence"] or len(target) != assay["target_length"]:
            raise ValueError("cohort target disagrees with metadata")
        source_reference = json.loads(contents[assay["reference_file"]])
        mapping = assay["mapping"]
        left = mapping["start"] - source_reference["start"]
        right = mapping["end"] - source_reference["start"]
        reference = source_reference["dna"].upper()[left:right]
        strand = mapping["target_reference_alignment_strand"]
        if strand == "-":
            reference = reference.translate(COMPLEMENT)[::-1]
        if strand not in {"+", "-"} or reference != target:
            raise ValueError("cohort target does not map exactly to frozen reference")
        if metadata["numVariants"] != assay["expected_source_rows"]:
            raise ValueError("cohort source row count disagrees with metadata")
        rows, counts = parse_rows(
            contents[assay["scores_file"]], target, assay_id,
            assay["expected_source_rows"], mapping)
        if not rows:
            raise ValueError("cohort assay has no eligible SNVs")
        assays[assay_id] = {"metadata": assay, "target": target, "rows": rows}
        accounting[assay_id] = counts
    if set(assays) != set(roles):
        raise ValueError("cohort inventory differs from frozen roles")
    return assays, {
        "directory": str(directory),
        "manifest_sha256": sha256(manifest_raw),
        "files": manifest["files"],
        "row_accounting": accounting,
    }


class SplitMix64:
    def __init__(self, state):
        self.state = state

    def next(self):
        mask = (1 << 64) - 1
        self.state = (self.state + 0x9E3779B97F4A7C15) & mask
        z = ((self.state ^ (self.state >> 30)) * 0xBF58476D1CE4E5B9) & mask
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & mask
        return (z ^ (z >> 31)) & mask

    def index(self, size):
        return self.next() % size


def grouped_metrics(rows: list[Observation], actual, predicted):
    y = np.asarray(actual, dtype=np.float64)
    p = np.asarray(predicted, dtype=np.float64)
    groups = {}
    for i, row in enumerate(rows):
        groups.setdefault(row.position, []).append(i)
    position_mse = []
    position_mae = []
    position_actual = []
    position_predicted = []
    position_losses = []
    exact_observed = [
        Fraction(row.score_raw) if row.score_raw and float(y[i]) == row.score
        else Fraction(str(float(y[i])))
        for i, row in enumerate(rows)
    ]
    exact_position_actual = []

    def exact_mean(values):
        # Mixed one/two/three-allele groups must not turn a constant predictor
        # into tiny rank differences through floating-point re-averaging.
        return float(values[0]) if np.all(values == values[0]) else float(np.mean(values))

    for position, indices in sorted(groups.items()):
        error = y[indices] - p[indices]
        mse = float(np.mean(error * error))
        mae = float(np.mean(np.abs(error)))
        position_mse.append(mse)
        position_mae.append(mae)
        position_actual.append(exact_mean(y[indices]))
        position_predicted.append(exact_mean(p[indices]))
        exact_position_actual.append(sum((exact_observed[i] for i in indices), Fraction()) / len(indices))
        position_losses.append({"position": position, "mse": mse, "mae": mae})
    return {
        "rows": len(rows),
        "positions": len(groups),
        "mse": float(np.mean(position_mse)),
        "mae": float(np.mean(position_mae)),
        "spearman_rows": correlation(rankdata(y, method="average"), rankdata(p, method="average")),
        "spearman_positions": correlation(
            rankdata(exact_position_actual, method="average"),
            rankdata(position_predicted, method="average")),
        "pearson_rows": correlation(y, p),
        "position_losses": position_losses,
    }


@dataclass
class CachedFit:
    train: np.ndarray
    test: np.ndarray
    train_scaled: np.ndarray
    test_scaled: np.ndarray
    mean: np.ndarray
    scale: np.ndarray
    factor: object
    columns: list[int]
    unique_test_scaled: np.ndarray
    test_inverse: np.ndarray

    @classmethod
    def build(cls, features, train, test, columns, alpha):
        train = np.asarray(train, dtype=np.int64)
        test = np.asarray(test, dtype=np.int64)
        x = features[train][:, columns]
        q = features[test][:, columns]
        mean = x.mean(axis=0)
        scale = x.std(axis=0, ddof=0)
        scale[scale <= 1e-12] = 1.0
        z = (x - mean) / scale
        if len(columns):
            _, unique_indices, inverse = np.unique(q, axis=0, return_index=True, return_inverse=True)
            unique_test = (q[unique_indices] - mean) / scale
        else:
            unique_test = np.empty((0, 0))
            inverse = np.zeros(len(q), dtype=int)
        q = (q - mean) / scale
        factor = None
        if len(columns):
            factor = cho_factor(z.T @ z + alpha * np.eye(len(columns)),
                                lower=True, check_finite=True)
        return cls(train, test, z, q, mean, scale, factor, columns, unique_test, inverse)

    def predict(self, all_y):
        y = np.asarray(all_y, dtype=np.float64)
        one_dimensional = y.ndim == 1
        if one_dimensional:
            y = y[:, None]
        training = y[self.train]
        mean_y = training.mean(axis=0)
        if self.factor is None:
            beta = np.zeros((0, y.shape[1]))
            prediction = np.broadcast_to(mean_y, (len(self.test), y.shape[1])).copy()
        else:
            rhs = self.train_scaled.T @ (training - mean_y)
            beta = cho_solve(self.factor, rhs, check_finite=True)
            # Equal feature vectors are the same mathematical input. Evaluate
            # each unique UNSCALED row once so BLAS tail/alignment kernels cannot
            # split their predictions by one ULP and invent Spearman ranks.
            unique_prediction = mean_y + self.unique_test_scaled @ beta
            prediction = unique_prediction[self.test_inverse]
        if one_dimensional:
            return prediction[:, 0], float(mean_y[0]), beta[:, 0]
        return prediction, mean_y, beta


def rotation_targets(rows, y, protocol):
    """One shared table for HBB plus primary loci, in frozen stream order."""
    seed = protocol["spatial_null"]["seed"]
    replicates = protocol["spatial_null"]["replicates"]
    rng = SplitMix64(seed)
    classes = protocol["features"]["mutation_identity_order"]
    assays = sorted(set(row.assay_id for row in rows))
    strata = []
    unchanged = []
    for assay in assays:
        for identity in classes:
            members = [i for i, row in enumerate(rows)
                       if row.assay_id == assay and row.reference + ">" + row.alternate == identity]
            members.sort(key=lambda i: (rows[i].position, rows[i].accession))
            if len(members) < 2:
                unchanged.append({"assay_id": assay, "identity": identity, "rows": len(members)})
            strata.append(members)
    result = np.tile(np.asarray(y)[:, None], (1, replicates))
    shifts = []
    for repetition in range(replicates):
        offsets = []
        for members in strata:
            if len(members) < 2:
                offsets.append(0)
                continue
            offset = 1 + rng.index(len(members) - 1)
            offsets.append(offset)
            source = np.roll(np.asarray(members), -offset)
            result[members, repetition] = np.asarray(y)[source]
        shifts.append(offsets)
    return result, {"seed": seed, "unchanged_strata": unchanged, "offsets": shifts}


def study_splits(rows, assays, protocol):
    hbb_id = "urn:mavedb:00000018-a-1"
    hbb = [i for i, row in enumerate(rows) if row.assay_id == hbb_id]
    hbb_positions = np.array([rows[i].position for i in hbb])
    length = len(assays[hbb_id]["target"])
    boundaries = [k * length // 5 for k in range(6)]
    blocks = [(boundaries[k] + 1, boundaries[k + 1]) for k in range(5)]
    plans = []
    for fold in blocked_folds(hbb_positions, blocks, 20):
        plans.append({
            "scheme": "hbb_guarded_position_cv", "fold": str(fold["fold"]),
            "train": [hbb[i] for i in fold["train_indices"]],
            "test": [hbb[i] for i in fold["test_indices"]],
            "guard": [hbb[i] for i in fold["guard_indices"]],
            "block_start": fold["start"], "block_end": fold["end"],
        })
    primary = protocol["evaluation"]["primary_leave_element_out"]["test_assays"]
    for assay_id in primary:
        test = [i for i, row in enumerate(rows) if row.assay_id == assay_id]
        train = [i for i, row in enumerate(rows)
                 if row.assay_id in primary and row.assay_id != assay_id]
        plans.append({"scheme": "primary_leave_element_out", "fold": assay_id,
                      "train": train, "test": test, "guard": []})
        plans.append({"scheme": "hbb_to_external_transfer", "fold": assay_id,
                      "train": hbb, "test": test, "guard": []})
    random_fold = {}
    for i in hbb:
        row = rows[i]
        text = f"{row.assay_id}|{row.position}|{row.reference}>{row.alternate}|row-random-v1"
        random_fold[i] = int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little") % 5
    for fold in range(5):
        test = [i for i in hbb if random_fold[i] == fold]
        train = [i for i in hbb if random_fold[i] != fold]
        plans.append({"scheme": "row_random_diagnostic", "fold": str(fold),
                      "train": train, "test": test, "guard": []})
    return plans


def summarize_predictions(rows, y, family_predictions, protocol):
    result = {}
    for scheme, by_model in family_predictions.items():
        by_assay = {}
        for assay in sorted(set(row.assay_id for row in rows)):
            indices = [i for i, row in enumerate(rows) if row.assay_id == assay]
            if not np.isfinite(by_model["mean"][indices]).all():
                continue
            local_rows = [rows[i] for i in indices]
            by_assay[assay] = {
                model: grouped_metrics(local_rows, y[indices], prediction[indices])
                for model, prediction in by_model.items()
            }
            mean_mse = by_assay[assay]["mean"]["mse"]
            for model in by_model:
                by_assay[assay][model]["relative_mse_vs_mean"] = (
                    None if mean_mse == 0 else by_assay[assay][model]["mse"] / mean_mse)
        aggregate = {}
        for model in by_model:
            aggregate[model] = {
                "macro_mse": float(np.mean([values[model]["mse"] for values in by_assay.values()])),
                "macro_mae": float(np.mean([values[model]["mae"] for values in by_assay.values()])),
            }
        comparisons = {}
        for comparison in protocol["comparisons"]:
            candidate, baseline = comparison["candidate"], comparison["baseline"]
            comparisons[comparison["id"]] = {
                "candidate": candidate, "baseline": baseline,
                "observed_improvement": (
                    aggregate[baseline]["macro_mse"] - aggregate[candidate]["macro_mse"]),
            }
        result[scheme] = {"assays": by_assay, "models": aggregate,
                          "comparisons": comparisons}
    return result


def bootstrap_comparisons(summary, protocol):
    repetitions = protocol["resampling"]["replicates"]
    families = {}
    for scheme in ["hbb_guarded_position_cv", "primary_leave_element_out"]:
        source = summary[scheme]
        if scheme == "hbb_guarded_position_cv":
            config = protocol["resampling"]["hbb_position_block_bootstrap"]
            hbb = source["assays"]["urn:mavedb:00000018-a-1"]
            positions = [item["position"] for item in hbb["mean"]["position_losses"]]
            blocks = {}
            for i, position in enumerate(positions):
                blocks.setdefault((position - 1) // config["block_bp"], []).append(i)
            groups = [members for _, members in sorted(blocks.items())]
            deltas = {}
            for comparison in protocol["comparisons"]:
                left = hbb[comparison["baseline"]]["position_losses"]
                right = hbb[comparison["candidate"]]["position_losses"]
                deltas[comparison["id"]] = np.array(
                    [a["mse"] - b["mse"] for a, b in zip(left, right)])
        else:
            config = protocol["resampling"]["external_element_bootstrap"]
            assay_ids = sorted(source["assays"])
            groups = [[i] for i in range(len(assay_ids))]
            deltas = {}
            for comparison in protocol["comparisons"]:
                deltas[comparison["id"]] = np.array([
                    source["assays"][assay][comparison["baseline"]]["mse"]
                    - source["assays"][assay][comparison["candidate"]]["mse"]
                    for assay in assay_ids])
        rng = SplitMix64(config["seed"])
        draws = []
        samples = {key: [] for key in deltas}
        for _ in range(repetitions):
            sampled_groups = [rng.index(len(groups)) for _ in groups]
            draws.append(sampled_groups)
            indices = [index for group in sampled_groups for index in groups[group]]
            for key, values in deltas.items():
                samples[key].append(float(values[indices].mean()))
        families[scheme] = {
            "seed": config["seed"],
            "replicates": repetitions,
            "resampled_group_indices": draws,
            "comparisons": {
                key: {"statistics": values,
                      "interval_type7_025_975": np.quantile(values, [0.025, 0.975], method="linear").tolist()}
                for key, values in samples.items()
            },
        }
    return families


def null_family_metrics(rows, actual, predicted):
    """Vectorized fixed-OOF metrics for all null columns, with position weighting."""
    included = np.flatnonzero(np.isfinite(predicted[:, 0]))
    assays = {}
    for i in included:
        row = rows[i]
        assays.setdefault(row.assay_id, {}).setdefault(row.position, []).append(i)
    per_assay = []
    for groups in assays.values():
        per_position = []
        for members in groups.values():
            error = actual[members] - predicted[members]
            per_position.append(np.mean(error * error, axis=0))
        per_assay.append(np.mean(per_position, axis=0))
    return np.mean(per_assay, axis=0)


def spatial_null(rows, y, caches, summary, protocol):
    rotated, rotation_info = rotation_targets(rows, y, protocol)
    repetitions = rotated.shape[1]
    families = ["hbb_guarded_position_cv", "primary_leave_element_out"]
    predictions = {
        family: {model["id"]: np.full((len(rows), repetitions), np.nan)
                 for model in protocol["models"]}
        for family in families
    }
    for plan, model, fit in caches:
        if plan["scheme"] not in families:
            continue
        prediction, _, _ = fit.predict(rotated)
        predictions[plan["scheme"]][model][fit.test] = prediction
    result = {}
    model_mses = {}
    for family in families:
        losses = {model: null_family_metrics(rows, rotated, values)
                  for model, values in predictions[family].items()}
        model_mses[family] = {model: values.tolist() for model, values in losses.items()}
        comparisons = {}
        for comparison in protocol["comparisons"]:
            key = comparison["id"]
            statistics = losses[comparison["baseline"]] - losses[comparison["candidate"]]
            observed = summary[family]["comparisons"][key]["observed_improvement"]
            count = int(np.sum(statistics >= observed))
            tail = (1 + count) / (1 + repetitions)
            comparisons[key] = {
                "observed_improvement": observed,
                "statistics": statistics.tolist(),
                "null_greater_equal_observed": count,
                "tail_fraction": tail,
                "bonferroni_six": min(1.0, tail * 6.0),
            }
        ordered = sorted(comparisons, key=lambda key: comparisons[key]["tail_fraction"])
        previous = 0.0
        for rank, key in enumerate(ordered):
            adjusted = min(1.0, max(previous, (len(ordered) - rank)
                                    * comparisons[key]["tail_fraction"]))
            comparisons[key]["holm_within_family"] = adjusted
            previous = adjusted
        result[family] = comparisons
    return {"rotation": rotation_info, "families": result, "model_mses": model_mses}


def compare_rust(independent_directory, rust_path, output):
    """Check every target/fold/preprocessing vector before declaring reproduction."""
    root = pathlib.Path(independent_directory)
    manifest = json.loads((root / "artifacts.json").read_text())
    for name, expected in manifest["files"].items():
        if pathlib.Path(name).name != name or sha256((root / name).read_bytes()) != expected:
            raise ValueError("independent artifact integrity check failed")
    independent = json.loads((root / "replication.json").read_text())
    own_predictions = json.loads((root / "predictions.json").read_text())
    own_folds = json.loads((root / "folds.json").read_text())
    own_fits = json.loads((root / "fits.json").read_text())
    own_features = json.loads((root / "features.json").read_text())
    own_resampling = json.loads((root / "resampling.json").read_text())
    inputs = json.loads((root / "inputs.json").read_text())
    raw_rust = pathlib.Path(rust_path).read_bytes()
    decoded_rust = gzip.decompress(raw_rust) if raw_rust.startswith(b"\x1f\x8b") else raw_rust
    rust = json.loads(decoded_rust)
    if rust["protocol_sha256"] != FROZEN_PROTOCOL_SHA256 or independent["protocol_sha256"] != FROZEN_PROTOCOL_SHA256:
        raise ValueError("comparison protocols differ")
    if rust["motif_manifest_sha256"] != inputs["motifs"]["manifest_sha256"]:
        raise ValueError("comparison motif bundles differ")
    prediction_tolerance = 1e-7
    metric_tolerance = 1e-8
    failures = []
    maxima = {}
    checked = {}

    def fail(category, identity, expected, actual):
        failures.append({"category": category, "identity": str(identity),
                         "independent": expected, "rust": actual})

    def numeric(category, identity, expected, actual, tolerance=metric_tolerance):
        checked[category] = checked.get(category, 0) + 1
        if expected is None or actual is None:
            if expected is not actual:
                fail(category, identity, expected, actual)
            return
        a = np.asarray(expected, dtype=float)
        b = np.asarray(actual, dtype=float)
        if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
            fail(category, identity, expected, actual)
            return
        error = float(np.max(np.abs(a - b))) if a.size else 0.0
        maxima[category] = max(maxima.get(category, 0.0), error)
        if error > tolerance:
            fail(category, identity, expected, actual)

    def equal(category, identity, expected, actual):
        checked[category] = checked.get(category, 0) + 1
        if expected != actual:
            fail(category, identity, expected, actual)

    vocabulary = own_features["vocabulary"]
    for model, columns in own_features["models"].items():
        equal("feature_schema", model, [vocabulary[i] for i in columns], rust["feature_schema"][model])
    expected = {}
    key_to_accession = {}
    for item in own_predictions:
        key = (item["scheme"], item["model"], item["assay_id"],
               item["position"], item["reference"], item["alternate"])
        if key in expected:
            raise ValueError("duplicate independent target prediction")
        expected[key] = item
        key_to_accession[item["assay_id"] + "|" + item["hgvs"]] = item["accession"]
    own_fold_map = {(item["scheme"], item["fold"]): item for item in own_folds}
    own_fit_map = {(item["scheme"], item["model"], item["fold"]): item for item in own_fits}
    observed_keys = set()
    equal("evaluation_families", "all", sorted(independent["evaluation"]),
          sorted(item["family"] for item in rust["analyses"]))
    for analysis in rust["analyses"]:
        family = analysis["family"]
        fold_keys = {}
        for fold in analysis["folds"]:
            index = fold["fold_index_0based"]
            key = str(index) if family in {"hbb_guarded_position_cv", "row_random_diagnostic"} else fold["target_assay"]
            fold_keys[index] = key
            own = own_fold_map[(family, key)]
            for ours, theirs in [("train_ids", "train_ids"), ("test_ids", "test_ids"),
                                 ("guard_ids", "excluded_guard_ids")]:
                expected_ids = sorted(key_to_accession[key] for key in own[ours])
                equal("fold_membership", (family, key, ours), expected_ids, sorted(fold[theirs]))
        for model in analysis["models"]:
            name = model["model"]
            own_summary = independent["evaluation"][family]
            numeric("macro_mse", (family, name),
                    own_summary["models"][name]["macro_mse"], model["macro_assay_position_mse"])
            numeric("macro_mae", (family, name),
                    own_summary["models"][name]["macro_mae"], model["macro_assay_position_mae"])
            for fit in model["fits"]:
                fold = fold_keys[fit["fold_index_0based"]]
                own = own_fit_map[(family, name, fold)]
                theirs = fit["fit"]
                equal("fit_dimensions", (family, name, fold), len(own["feature_names"]), theirs["feature_count"])
                equal("fit_training_rows", (family, name, fold), own["training_rows"], theirs["training_rows"])
                numeric("fit_means", (family, name, fold), own["training_feature_mean"], theirs["training_feature_means"])
                numeric("fit_scales", (family, name, fold), own["training_feature_scale"], theirs["training_feature_scales"])
                numeric("fit_intercepts", (family, name, fold), own["training_target_mean"], theirs["intercept"])
                numeric("fit_coefficients", (family, name, fold), own["coefficients_scaled"], theirs["standardized_coefficients"])
            for prediction in model["predictions"]:
                key = (family, name, prediction["assay_id"], prediction["position_1based"],
                       prediction["reference"], prediction["alternate"])
                if key in observed_keys:
                    fail("duplicate_rust_prediction", key, "unique", "duplicate")
                observed_keys.add(key)
                if key not in expected:
                    fail("unexpected_rust_target", key, "absent", "present")
                    continue
                own = expected[key]
                equal("target_accessions", key, own["accession"], prediction["id"])
                equal("target_fold", key, own["fold"], fold_keys[prediction["fold_index_0based"]])
                numeric("observations", key, own["observed"], prediction["observed"], 1e-12)
                numeric("predictions", key, own["predicted"], prediction["predicted"], prediction_tolerance)
            for assay in model["per_assay_metrics"]:
                own = own_summary["assays"][assay["assay_id"]][name]
                for ours, theirs in [
                    ("mse", "position_weighted_mse"), ("mae", "position_weighted_mae"),
                    ("relative_mse_vs_mean", "relative_mse_versus_training_mean")]:
                    numeric("per_assay_" + ours, (family, name, assay["assay_id"]), own[ours], assay[theirs])
                numeric("spearman_rows", (family, name, assay["assay_id"]),
                        own["spearman_rows"], assay["row_metrics"]["spearman"])
                numeric("spearman_positions", (family, name, assay["assay_id"]),
                        own["spearman_positions"], assay["position_mean_metrics"]["spearman"])
    equal("prediction_target_inventory", "all", len(expected), len(observed_keys))
    missing = sorted(set(expected) - observed_keys)
    if missing:
        fail("missing_rust_predictions", "all", missing[:20], [])
    comparisons = ["delta3mer_vs_chemistry", "motif_vs_chemistry", "full_vs_coordinate"]
    for own_family, rust_family in [
        ("hbb_guarded_position_cv", "hbb_position_block"),
        ("primary_leave_element_out", "primary_element")]:
        own = own_resampling["bootstrap"][own_family]
        theirs = rust["resampling"][rust_family]
        equal("bootstrap_draws", own_family, own["resampled_group_indices"], theirs["resampling_indices"])
        for index, key in enumerate(comparisons):
            numeric("bootstrap_statistics", (own_family, key),
                    own["comparisons"][key]["statistics"],
                    [row[index] for row in theirs["sampled_statistics"]])
            summary = next(item for item in theirs["comparisons"] if item["comparison_id"] == key)
            numeric("bootstrap_intervals", (own_family, key),
                    own["comparisons"][key]["interval_type7_025_975"],
                    [summary["percentile_95_lower"], summary["percentile_95_upper"]])
    model_order = ["mean", "chemistry", "coordinate", "delta3mer", "motif", "full"]
    for own_family, loss_field, model_field in [
        ("hbb_guarded_position_cv", "hbb_loss_improvements", "hbb_model_mses"),
        ("primary_leave_element_out", "primary_loss_improvements", "primary_macro_model_mses")]:
        for index, key in enumerate(comparisons):
            numeric("spatial_null_statistics", (own_family, key),
                    own_resampling["spatial_null"]["families"][own_family][key]["statistics"],
                    [row[loss_field][index] for row in rust["spatial_null"]["statistics"]])
        if "model_mses" in own_resampling["spatial_null"]:
            for index, model in enumerate(model_order):
                numeric("spatial_null_model_mses", (own_family, model),
                        own_resampling["spatial_null"]["model_mses"][own_family][model],
                        [row[model_field][index] for row in rust["spatial_null"]["statistics"]])
    result = {
        "schema_version": 1,
        "matched": not failures,
        "primary_predictions_losses_and_controls_matched": not any(
            item["category"] not in {"spearman_rows", "spearman_positions"} for item in failures),
        "secondary_exact_rank_metrics_matched": not any(
            item["category"] in {"spearman_rows", "spearman_positions"} for item in failures),
        "protocol_sha256": FROZEN_PROTOCOL_SHA256,
        "independent_directory": str(root),
        "independent_artifacts_sha256": sha256((root / "artifacts.json").read_bytes()),
        "rust_study_path": str(rust_path), "rust_study_sha256": sha256(raw_rust),
        "rust_decoded_study_sha256": sha256(decoded_rust),
        "prediction_absolute_tolerance": prediction_tolerance,
        "metric_absolute_tolerance": metric_tolerance,
        "checked": checked, "maximum_absolute_errors": maxima,
        "failure_count": len(failures), "failures": failures,
        "scope": "Independent numerical implementation agreement, not independent biological replication",
    }
    output = pathlib.Path(output)
    if output.exists():
        raise ValueError("comparison output already exists")
    output.mkdir(parents=True)
    (output / "comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: result[key] for key in [
        "matched", "primary_predictions_losses_and_controls_matched",
        "secondary_exact_rank_metrics_matched", "checked",
        "maximum_absolute_errors", "failure_count"]}, indent=2))
    return result


def context_pairs(assays, protocol):
    result = []
    for pair in protocol["validation_context_analysis"]["pairs"]:
        left, right = assays[pair["left"]], assays[pair["right"]]
        a = {row.genomic_key: row for row in left["rows"]}
        b = {row.genomic_key: row for row in right["rows"]}
        keys = sorted(set(a) & set(b))
        y = np.array([a[key].score for key in keys])
        p = np.array([b[key].score for key in keys])
        result.append({
            **pair, "matched_variants": len(keys),
            "unmatched_left": len(set(a) - set(b)), "unmatched_right": len(set(b) - set(a)),
            "mse": float(np.mean((y - p) ** 2)) if keys else None,
            "mae": float(np.mean(np.abs(y - p))) if keys else None,
            "signed_mean_left_minus_right": float(np.mean(y - p)) if keys else None,
            "pearson": correlation(y, p) if keys else None,
            "spearman": correlation(rankdata(y), rankdata(p)) if keys else None,
            # Two contexts, variant-specific common predictor: half the pair
            # difference on each context, hence MSE floor = mean(diff²)/4.
            "finite_context_blind_mse_lower_bound": float(np.mean((y - p) ** 2) / 4.0) if keys else None,
        })
    return result


FROZEN_PROTOCOL_SHA256 = "c8e97fabef2afe8ccaa5a01a7ef7bb952f9c9f15857e5c540af112a66e78cb58"


def run_study(protocol_path, output, cohort_directory=None, motif_directory=None):
    implementation_snapshot = pathlib.Path(__file__).read_bytes()
    protocol_raw = pathlib.Path(protocol_path).read_bytes()
    protocol_hash = sha256(protocol_raw)
    protocol = json.loads(protocol_raw)
    if (protocol_hash != FROZEN_PROTOCOL_SHA256 or protocol.get("status") != "frozen"
            or protocol.get("protocol_id") != "regulatory-prediction-v1"):
        raise ValueError("this implementation requires the exact frozen v1 protocol")
    output = pathlib.Path(output)
    if output.exists():
        raise ValueError("output already exists; choose a new independent run directory")
    profiles, motif_inputs = checked_motifs(
        pathlib.Path(motif_directory or protocol["inputs"]["motif_directory"]),
        protocol["features"]["motif_ids"])
    assays, cohort_inputs = checked_cohort(
        pathlib.Path(cohort_directory or protocol["inputs"]["cohort_directory"]),
        protocol, protocol_hash)
    hbb_target, hbb_rows, hbb_inputs = checked_assay(
        pathlib.Path(protocol["inputs"]["hbb_directory"]), protocol_frozen=True)
    hbb_id = "urn:mavedb:00000018-a-1"
    if hbb_target != assays[hbb_id]["target"] or hbb_rows != assays[hbb_id]["rows"]:
        raise ValueError("cohort HBB copy differs from independently pinned original")
    active = set(protocol["evaluation"]["primary_leave_element_out"]["test_assays"]) | {hbb_id}
    rows = []
    matrices = []
    feature_info = {}
    vocabulary = None
    models = None
    for assay_id in sorted(active):
        assay = assays[assay_id]
        print(f"Independent features: {assay_id} ({len(assay['rows'])} eligible rows)", flush=True)
        matrix, model_columns, columns = feature_matrix(
            assay["target"], assay["rows"], profiles, protocol)
        if vocabulary is not None and columns != vocabulary:
            raise ValueError("inconsistent feature vocabularies")
        vocabulary, models = columns, model_columns
        rows.extend(assay["rows"])
        matrices.append(matrix)
        feature_info[assay_id] = {
            "rows": len(matrix), "columns": matrix.shape[1],
            "little_endian_f64_sha256": sha256(matrix.astype("<f8").tobytes()),
        }
    features = np.vstack(matrices)
    y = np.array([row.score for row in rows])
    plans = study_splits(rows, assays, protocol)
    schemes = sorted(set(plan["scheme"] for plan in plans))
    family_predictions = {
        scheme: {model: np.full(len(rows), np.nan) for model in models}
        for scheme in schemes
    }
    predictions = []
    fits = []
    fold_artifacts = []
    caches = []
    for plan in plans:
        print(f"Independent fits: {plan['scheme']} / {plan['fold']}", flush=True)
        train_positions = {}
        for i in plan["train"]:
            train_positions.setdefault(rows[i].assay_id, set()).add(rows[i].position)
        test_positions = {}
        for i in plan["test"]:
            test_positions.setdefault(rows[i].assay_id, set()).add(rows[i].position)
        exact_overlap = 0
        support_overlap = 0
        for assay_id, positions in test_positions.items():
            train = train_positions.get(assay_id, set())
            exact_overlap += sum(position in train for position in positions)
            support_overlap += sum(any(abs(position - other) <= 20 for other in train)
                                   for position in positions)
        fold_artifacts.append({
            **{key: value for key, value in plan.items()
               if key not in {"train", "test", "guard"}},
            "train_ids": [rows[i].key for i in plan["train"]],
            "test_ids": [rows[i].key for i in plan["test"]],
            "guard_ids": [rows[i].key for i in plan["guard"]],
            "test_positions_also_in_training": exact_overlap,
            "test_positions_with_possible_feature_support_overlap": support_overlap,
        })
        for model, columns in models.items():
            fit = CachedFit.build(
                features, plan["train"], plan["test"], columns,
                protocol["fit"]["ridge_alpha"])
            prediction, mean_y, coefficients = fit.predict(y)
            caches.append((plan, model, fit))
            family_predictions[plan["scheme"]][model][fit.test] = prediction
            fits.append({
                "scheme": plan["scheme"], "fold": plan["fold"], "model": model,
                "training_rows": len(fit.train), "test_rows": len(fit.test),
                "feature_names": [vocabulary[column] for column in columns],
                "training_feature_mean": fit.mean.tolist(),
                "training_feature_scale": fit.scale.tolist(),
                "training_target_mean": mean_y,
                "coefficients_scaled": coefficients.tolist(),
            })
            for i, value in zip(fit.test, prediction):
                row = rows[i]
                predictions.append({
                    "scheme": plan["scheme"], "fold": plan["fold"], "model": model,
                    "assay_id": row.assay_id, "accession": row.accession,
                    "hgvs": row.hgvs, "position": row.position,
                    "reference": row.reference, "alternate": row.alternate,
                    "observed": row.score, "predicted": float(value),
                })
    summary = summarize_predictions(rows, y, family_predictions, protocol)
    print("Independent block/element bootstraps", flush=True)
    bootstrap = bootstrap_comparisons(summary, protocol)
    print("Independent shared-table spatial null (99 cached refits)", flush=True)
    null = spatial_null(rows, y, caches, summary, protocol)
    contexts = context_pairs(assays, protocol)
    software = {
        "python": sys.version,
        "implementation": platform.python_implementation(),
        "numpy": np.__version__, "scipy": scipy.__version__,
        "script_sha256": sha256(implementation_snapshot),
    }
    report = {
        "schema_version": 1,
        "purpose": "Independent numerical implementation, not independent experimental replication",
        "protocol_sha256": protocol_hash,
        "software": software,
        "completed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "eligible_active_rows": len(rows),
        "assays": sorted(active),
        "model_dimensions": {model: len(columns) for model, columns in models.items()},
        "evaluation": summary,
        "context_pairs": contexts,
        "limitations": [
            "HBB was previously explored; this is a retrospective within-element diagnostic.",
            "Nine external promoter loci share one study and source processing; they are not nine independent studies.",
            "Motif deltas are dual-strand sequence similarity changes, not binding or expression measurements.",
            "Source target orientation is kept for chemistry/kmer baselines; dual-strand motif deltas are reverse-complement invariant.",
            "The physical HBB insert is native-minus relative to the plus coordinate target; reporter TSS/transcript identity is not established.",
            "Bootstrap intervals and stratified spatial-null tails are conditional robustness diagnostics, not calibrated biological significance.",
        ],
        "rank_arithmetic": {
            "observed_position_means": "Exact Fraction(raw CSV decimal score) sums/counts ranked directly; no float regrouping or epsilon.",
            "predicted_identical_inputs": "Evaluate each identical unscaled feature row once per fit, then broadcast exactly; no score quantization.",
            "predicted_position_means": "Preserve identical values exactly; otherwise NumPy arithmetic mean in ALT-sorted row order.",
            "remaining_near_ties": "No epsilon or quantization. Cross-implementation rank mismatches remain failures if exact float rank inputs differ.",
        },
    }
    artifacts = {
        "replication.json": report,
        "predictions.json": predictions,
        "folds.json": fold_artifacts,
        "fits.json": fits,
        "features.json": {
            "vocabulary": vocabulary, "models": models, "assays": feature_info,
        },
        "inputs.json": {"cohort": cohort_inputs, "motifs": motif_inputs, "hbb": hbb_inputs},
        "resampling.json": {"bootstrap": bootstrap, "spatial_null": null},
    }
    output.mkdir(parents=True)
    commitments = {}
    for name, artifact in artifacts.items():
        content = (json.dumps(artifact, indent=2, allow_nan=False) + "\n").encode()
        (output / name).write_bytes(content)
        commitments[name] = sha256(content)
    (output / "implementation.py").write_bytes(implementation_snapshot)
    commitments["implementation.py"] = sha256(implementation_snapshot)
    (output / "artifacts.json").write_text(json.dumps({
        "schema_version": 1, "files": commitments,
        "protocol_sha256": protocol_hash,
    }, indent=2) + "\n")
    print(f"Independent replication saved to {output}", flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", default="experiments/regulatory-prediction-v1.json")
    parser.add_argument("--out", required=True)
    parser.add_argument("--cohort")
    parser.add_argument("--motifs")
    parser.add_argument("--compare-rust", help="Rust StudyResult JSON; compare without refitting")
    parser.add_argument("--independent-run", help="Completed independent output for --compare-rust")
    args = parser.parse_args()
    if args.compare_rust:
        if not args.independent_run:
            parser.error("--compare-rust requires --independent-run")
        result = compare_rust(args.independent_run, args.compare_rust, args.out)
        if not result["matched"]:
            raise SystemExit(2)
    else:
        run_study(args.protocol, args.out, args.cohort, args.motifs)


if __name__ == "__main__":
    main()
