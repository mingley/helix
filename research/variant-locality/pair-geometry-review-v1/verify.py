#!/usr/bin/env python3
"""Independent full-sequence check of every published nearby-pair class.

Run from the repository root; no measurements are read. Generates a JSON receipt.
"""
import collections
import csv
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path

ROOT = Path('research/variant-locality')
TABLE = ROOT / 'pair-geometry-v1/all-pairs.csv'
SUMMARY = ROOT / 'pair-geometry-v1/summary.json'
FITS = ROOT / 'model-envelope-v1/saved-fits.json'
OUT = ROOT / 'pair-geometry-review-v1/receipt.json'
BASES = 'ACGT'
COMP = str.maketrans('ACGT', 'TGCA')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def feature_counts(seq, index):
    result = [0] * 33
    for j in range(198):
        word = seq[j:j+3]
        reverse = word.translate(COMP)[::-1]
        result[index[min(word, reverse)]] += 1
    result[32] = sum(b in 'CG' for b in seq)
    return result


def replace(seq, pos, base):
    return seq[:pos] + base + seq[pos+1:]


def main():
    fits = json.loads(FITS.read_text())
    summary = json.loads(SUMMARY.read_text())
    words = sorted({min(w, w.translate(COMP)[::-1]) for w in
                    map(''.join, itertools.product(BASES, repeat=3))})
    assert fits['feature_names'] == ['rc3mer_' + w for w in words] + ['gc_fraction']
    index = {w: i for i, w in enumerate(words)}
    models = fits['models']
    assert len(fits['raw_slopes']) == 5
    # Work in exact rational arithmetic on the five saved binary64 slopes.
    slopes = {m: [sum((Fraction.from_float(f[m][i]) for f in fits['raw_slopes']),
                      Fraction()) / 5 for i in range(33)] for m in models}
    expected = set()
    for d, n, p, q in ((1, 4, 1, 2), (2, 3, 0, 2)):
        for ref in map(''.join, itertools.product(BASES, repeat=n)):
            expected.update((d, ref, a, b) for a in BASES for b in BASES
                            if a != ref[p] and b != ref[q])
    seen = set()
    distributions = {1: collections.Counter(), 2: collections.Counter()}
    distinct = {1: set(), 2: set()}
    max_exact = {(d, m): Fraction() for d in (1, 2) for m in models}
    worst = 0.0
    values_checked = 0
    with TABLE.open(newline='') as handle:
        for row in csv.DictReader(handle):
            d = int(row['distance_nt'])
            ref = row['reference_minimal_context']
            p, q = int(row['first_index']), int(row['second_index'])
            a, b = row['alt_first'], row['alt_second']
            key = (d, ref, a, b)
            assert key in expected and key not in seen and (p, q) == ((1, 2) if d == 1 else (0, 2))
            seen.add(key)
            assert row['double_context'] == replace(replace(ref, p, a), q, b)
            # Full 200nt counts and four separate score vectors, with a flank
            # unrelated to the generator's A...C flank.
            left = 'ACGT' * 24
            original = left + ref + ('TGCA' * 26)[:104-len(ref)]
            assert len(original) == 200
            one = replace(original, len(left)+p, a)
            two = replace(original, len(left)+q, b)
            both = replace(one, len(left)+q, b)
            vectors = [feature_counts(x, index) for x in (both, original, one, two)]
            delta = tuple(x+y-z-w for x, y, z, w in zip(*vectors))
            assert delta == tuple(map(int, row['interaction_count_vector'].split(',')))
            l1 = sum(map(abs, delta))
            assert l1 == int(row['count_l1']) and delta[-1] == 0 and any(delta)
            distributions[d][l1] += 1
            distinct[d].add(delta)
            for model in models:
                exact = sum((slopes[model][i] * delta[i] / (198 if i < 32 else 200)
                             for i in range(33)), Fraction())
                observed = float(row[model])
                assert math.isfinite(observed)
                err = abs(observed - float(exact))
                worst = max(worst, err)
                assert err < 1e-15
                values_checked += 1
                max_exact[d, model] = max(max_exact[d, model], abs(exact))
    assert seen == expected and len(seen) == 2880
    for d in (1, 2):
        item = summary['pair_summary'][str(d)]
        assert item['distinct_interaction_vectors'] == len(distinct[d])
        assert item['count_l1_histogram'] == {str(k): v for k, v in sorted(distributions[d].items())}
        for m in models:
            result = summary['prediction_extrema'][str(d)][m]
            assert max_exact[d, m] == Fraction(int(result['exact_max_numerator']),
                                               int(result['exact_max_denominator']))
    receipt = {'method': 'Four complete 200nt counts per class with independent flank; exact rational scores from five saved binary64 coefficient vectors',
               'input_sha256': {str(p): digest(p) for p in (TABLE, SUMMARY, FITS)},
               'classes_checked': len(seen), 'feature_counts_checked': len(seen)*33,
               'predictions_checked': values_checked, 'maximum_binary64_error': worst,
               'exact_summary_extrema_checked': len(max_exact),
               'all_checks_passed': True, 'outcomes_read': False}
    OUT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
