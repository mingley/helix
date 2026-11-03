#!/usr/bin/env python3
"""Check an exact GC/trimer identity on the frozen DNA panel.

Only sequences are analyzed. The panel JSON parser also instantiates previously
explored observation fields; those fields are not accessed by this analysis.
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
PANEL_SHA256 = "fce6b2a1aaf5aa4d530a2e298dccb840bc7c594a05b6831fdb0c6c0f9e2a9be3"
PUBLISHED = Path(__file__).with_name("summary.json")
BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
WORDS = tuple("".join(word) for word in itertools.product(BASES, repeat=3))
CANON = {word: min(word, word.translate(COMPLEMENT)[::-1]) for word in WORDS}
CLASSES = tuple(sorted(set(CANON.values())))
GC = {word: sum(base in "GC" for base in word) for word in CLASSES}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def features(sequence):
    """Return exact canonical trimer counts, GC count and boundary statistic."""
    if not isinstance(sequence, str) or len(sequence) < 4 or set(sequence) - set(BASES):
        raise ValueError("expected at least four uppercase A/C/G/T bases")
    counts = Counter(CANON[sequence[i:i + 3]] for i in range(len(sequence) - 2))
    g = sum(base in "GC" for base in sequence)
    boundary = (2 * (sequence[0] in "GC") + (sequence[1] in "GC") +
                (sequence[-2] in "GC") + 2 * (sequence[-1] in "GC"))
    weighted_trimers = sum(GC[word] * count for word, count in counts.items())
    if weighted_trimers + boundary != 3 * g:
        raise ValueError("GC/trimer identity failed")
    return counts, g, weighted_trimers, boundary


def exact_scores(sequence, weights, gc_weight):
    """Evaluate the original and boundary-reparameterized scores as rationals."""
    counts, g, _, boundary = features(sequence)
    length = len(sequence)
    old = sum((weights[word] * Fraction(n, length - 2) for word, n in counts.items()), Fraction())
    old += gc_weight * Fraction(g, length)
    new = sum(((weights[word] + gc_weight * Fraction((length - 2) * GC[word], 3 * length)) *
               Fraction(n, length - 2) for word, n in counts.items()), Fraction())
    new += gc_weight * Fraction(boundary, 3 * length)
    return old, new


def analyze():
    if sha256(PANEL) != PANEL_SHA256:
        raise ValueError("frozen sequence-panel SHA-256 mismatch")
    panel = json.loads(PANEL.read_text())
    records = panel["records"]
    if len(records) != 55276:
        raise ValueError("unexpected sequence-panel size")
    histogram = Counter()
    sum_g = sum_t = sum_g2 = sum_t2 = sum_gt = 0
    for row in records:
        sequence = row["sequence"]
        if len(sequence) != 200:
            raise ValueError("panel contains a non-200-base sequence")
        _, g, t, boundary = features(sequence)
        histogram[boundary] += 1
        sum_g += g
        sum_t += t
        sum_g2 += g * g
        sum_t2 += t * t
        sum_gt += g * t
    n = len(records)
    covariance = n * sum_gt - sum_g * sum_t
    var_g = n * sum_g2 - sum_g * sum_g
    var_t = n * sum_t2 - sum_t * sum_t
    if covariance <= 0 or var_g <= 0 or var_t <= 0:
        raise ValueError("unexpected degenerate or negative panel correlation")
    correlation_squared = Fraction(covariance * covariance, var_g * var_t)
    with localcontext() as ctx:
        ctx.prec = 30
        correlation = str(Decimal(covariance) / (Decimal(var_g) * Decimal(var_t)).sqrt())

    # A cyclic word represented with two different starting positions has the
    # same multiset of 3-mer windows but different terminal GC composition.
    first, second = "AAC" * 66 + "AA", "ACA" * 66 + "AC"
    c1, g1, _, b1 = features(first)
    c2, g2, _, b2 = features(second)
    if len(first) != 200 or len(second) != 200 or c1 != c2 or g1 == g2:
        raise ValueError("constructive example did not have the required features")
    return {
        "scope": "Existing DNA panel only; previously explored observation fields instantiated but not used",
        "panel_sha256": PANEL_SHA256,
        "implementation_sha256": sha256(Path(__file__)),
        "records": n,
        "sequence_length": 200,
        "identity": "3 * GC_count = GC_weighted_canonical_trimer_count + boundary_GC",
        "boundary_GC": "2*GC(first) + GC(second) + GC(penultimate) + 2*GC(last)",
        "boundary_count_histogram": {str(k): histogram[k] for k in range(7)},
        "GC_count_vs_weighted_trimer_count_pearson_decimal": correlation,
        # These exact integers exceed JavaScript's safe integer range. Keep
        # them as decimal strings so JSON consumers can parse them losslessly.
        "correlation_squared_exact": {"numerator": str(correlation_squared.numerator),
                                      "denominator": str(correlation_squared.denominator)},
        "correlation_sufficient_statistics": {"sum_G": sum_g, "sum_T": sum_t,
                                               "sum_G2": sum_g2, "sum_T2": sum_t2,
                                               "sum_GT": sum_gt},
        "same_trimers_different_GC_example": {
            "first_sequence": first, "second_sequence": second,
            "canonical_trimer_counts": dict(sorted(c1.items())),
            "first_GC_count": g1, "second_GC_count": g2,
            "first_boundary_GC": b1, "second_boundary_GC": b2,
        },
        "biological_effects_established": False,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--out", type=Path, help="new output directory")
    mode.add_argument("--verify", action="store_true", help="compare with committed summary")
    args = p.parse_args()
    result = analyze()
    if args.verify:
        if result != json.loads(PUBLISHED.read_text()):
            raise ValueError("published summary differs from panel reanalysis")
    else:
        args.out.mkdir(parents=True, exist_ok=False)
        (args.out / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"records": result["records"],
                      "boundary_count_histogram": result["boundary_count_histogram"],
                      "pearson": result["GC_count_vs_weighted_trimer_count_pearson_decimal"]}, sort_keys=True))


if __name__ == "__main__":
    main()
