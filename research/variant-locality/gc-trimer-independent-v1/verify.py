#!/usr/bin/env python3
"""Independently check the panel GC/trimer identity without importing its producer.

Only the existing panel's DNA sequences are used.  Parsing the JSON also
instantiates previously explored observation fields, which are not accessed.
"""
import argparse
from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PANEL = ROOT / "datasets/context-panel/panel.json"
PRIMARY = ROOT / "research/variant-locality/gc-trimer-identity-v1"
PANEL_SHA256 = "fce6b2a1aaf5aa4d530a2e298dccb840bc7c594a05b6831fdb0c6c0f9e2a9be3"
SCOPE = "Existing DNA panel only; previously explored observation fields instantiated but not used"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def canonical(trimer):
    opposite = "".join(COMPLEMENT[base] for base in reversed(trimer))
    return min(trimer, opposite)


def inspect_sequence(sequence):
    """Count raw bases in sliding windows, separately from canonical features."""
    if not isinstance(sequence, str) or len(sequence) < 4 or set(sequence) - set(COMPLEMENT):
        raise ValueError("expected at least four uppercase A/C/G/T bases")
    gc = sequence.count("G") + sequence.count("C")
    weighted_windows = sum(
        sum(base == "G" or base == "C" for base in sequence[i:i + 3])
        for i in range(len(sequence) - 2)
    )
    # The first/last positions appear in one window; their neighbors in two.
    edge = sum((3 - min(i + 1, 3, len(sequence) - i)) * (base == "G" or base == "C")
               for i, base in enumerate(sequence))
    if 3 * gc != weighted_windows + edge:
        raise ValueError("raw window identity failed")
    return gc, weighted_windows, edge


def canonical_counts(sequence):
    inspect_sequence(sequence)
    return Counter(canonical(sequence[i:i + 3]) for i in range(len(sequence) - 2))


def reparameterized_scores(sequence, gc_weight=Fraction(5, 7)):
    """Check a heterogeneous exact-rational score before and after reparameterization."""
    counts = canonical_counts(sequence)
    length = len(sequence)
    gc, _, edge = inspect_sequence(sequence)
    # Independent nonuniform rational weights on each canonical trimer class.
    def weight(word):
        code = sum((i + 1) * ("ACGT".index(base) + 1) for i, base in enumerate(word))
        return Fraction(code - 12, 11)
    original = sum((weight(word) * Fraction(count, length - 2)
                    for word, count in counts.items()), Fraction())
    original += gc_weight * Fraction(gc, length)
    converted = sum(((weight(word) + gc_weight * Fraction(
        (length - 2) * sum(base in "GC" for base in word), 3 * length)) *
        Fraction(count, length - 2) for word, count in counts.items()), Fraction())
    converted += gc_weight * Fraction(edge, 3 * length)
    return original, converted


def verify_coefficients(length):
    """Verify each possible count's coefficient, independently of its weight."""
    if length < 4:
        raise ValueError("length must be at least four")
    for letters in itertools.product("ACGT", repeat=3):
        word = "".join(letters)
        k = word.count("G") + word.count("C")
        reversed_k = sum(COMPLEMENT[base] in "GC" for base in word)
        if k != reversed_k:
            raise ValueError("reverse complement does not preserve GC")
        # A unit GC coefficient shifts the frequency coefficient by the left
        # factor. Multiplying it by count/(length-2) must yield k/(3*length).
        if Fraction((length - 2) * k, 3 * length) / (length - 2) != Fraction(k, 3 * length):
            raise ValueError("frequency coefficient reparameterization failed")


def compare_summary(observed, expected):
    if observed != expected:
        keys = sorted(key for key in set(observed) | set(expected)
                      if observed.get(key) != expected.get(key))
        raise ValueError(f"primary summary differs: {keys}")


def reconstruct(panel, source_sha):
    records = panel["records"]
    if len(records) != 55276:
        raise ValueError("unexpected panel record count")
    counts = Counter()
    samples = {}
    sg = st = sg2 = st2 = sgt = 0
    for row in records:
        sequence = row["sequence"]
        if not isinstance(sequence, str) or len(sequence) != 200:
            raise ValueError("unexpected panel sequence length")
        gc, windows, edge = inspect_sequence(sequence)
        counts[edge] += 1
        samples.setdefault(edge, sequence)
        sg += gc
        st += windows
        sg2 += gc * gc
        st2 += windows * windows
        sgt += gc * windows
    if set(counts) != set(range(7)):
        raise ValueError("expected to observe every edge count from zero through six")
    n = len(records)
    covariance = n * sgt - sg * st
    denominator = (n * sg2 - sg * sg) * (n * st2 - st * st)
    if covariance <= 0 or denominator <= 0:
        raise ValueError("panel correlation is not defined and positive")
    exact_squared = Fraction(covariance * covariance, denominator)
    # Compute sqrt of the exact reduced fraction at high precision, then round
    # independently to the published 30 significant decimal digits.
    with localcontext() as context:
        context.prec = 80
        high = (Decimal(exact_squared.numerator) / Decimal(exact_squared.denominator)).sqrt()
        context.prec = 30
        correlation = str(+high)
    first = "".join("AAC"[i % 3] for i in range(200))
    second = "".join("AAC"[(i + 1) % 3] for i in range(200))
    first_counts, second_counts = canonical_counts(first), canonical_counts(second)
    first_gc, _, first_edge = inspect_sequence(first)
    second_gc, _, second_edge = inspect_sequence(second)
    if first_counts != second_counts or first_gc == second_gc:
        raise ValueError("constructed equal-trimer/different-GC witness failed")
    for length in (4, 5, 200):
        verify_coefficients(length)
    for sequence in [first, second, *samples.values()]:
        original, converted = reparameterized_scores(sequence)
        if original != converted:
            raise ValueError("exact score reparameterization failed")
    return {
        "scope": SCOPE,
        "panel_sha256": PANEL_SHA256,
        "implementation_sha256": source_sha,
        "records": n,
        "sequence_length": 200,
        "identity": "3 * GC_count = GC_weighted_canonical_trimer_count + boundary_GC",
        "boundary_GC": "2*GC(first) + GC(second) + GC(penultimate) + 2*GC(last)",
        "boundary_count_histogram": {str(k): counts[k] for k in range(7)},
        "GC_count_vs_weighted_trimer_count_pearson_decimal": correlation,
        "correlation_squared_exact": {"numerator": exact_squared.numerator,
                                      "denominator": exact_squared.denominator},
        "correlation_sufficient_statistics": {"sum_G": sg, "sum_T": st,
                                               "sum_G2": sg2, "sum_T2": st2,
                                               "sum_GT": sgt},
        "same_trimers_different_GC_example": {
            "first_sequence": first, "second_sequence": second,
            "canonical_trimer_counts": dict(sorted(first_counts.items())),
            "first_GC_count": first_gc, "second_GC_count": second_gc,
            "first_boundary_GC": first_edge, "second_boundary_GC": second_edge,
        },
        "biological_effects_established": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, default=PANEL)
    parser.add_argument("--summary", type=Path, default=PRIMARY / "summary.json")
    parser.add_argument("--receipt", type=Path, required=True, help="new receipt file")
    args = parser.parse_args()
    if args.receipt.exists() or args.receipt.is_symlink():
        raise FileExistsError("receipt path already exists")
    if sha256(args.panel) != PANEL_SHA256:
        raise ValueError("unexpected panel SHA-256")
    panel = json.loads(args.panel.read_text())
    expected = reconstruct(panel, sha256(PRIMARY / "analyze.py"))
    observed = json.loads(args.summary.read_text())
    compare_summary(observed, expected)
    receipt = {
        "method": "Raw sliding-trimer characters and positional multiplicity; separate exact rational and Decimal calculations",
        "panel_sha256": sha256(args.panel),
        "primary_summary_sha256": sha256(args.summary),
        "primary_implementation_sha256": sha256(PRIMARY / "analyze.py"),
        "implementation_sha256": sha256(Path(__file__).resolve()),
        "records_checked": expected["records"],
        "panel_length": expected["sequence_length"],
        "boundary_histogram": expected["boundary_count_histogram"],
        "correlation_squared_exact": expected["correlation_squared_exact"],
        "correlation_decimal": expected["GC_count_vs_weighted_trimer_count_pearson_decimal"],
        "all_summary_fields_compared": True,
        "trimer_gc_identity_checked_per_panel_sequence": True,
        "rational_score_checks": "two constructed witnesses and one sequence from each of seven edge-count classes; coefficient algebra for all 64 trimers at lengths 4, 5 and 200",
        "old_panel_observation_fields_instantiated": True,
        "old_panel_observation_fields_used": False,
        "biological_effects_established": False,
        "all_checks_passed": True,
    }
    with args.receipt.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
