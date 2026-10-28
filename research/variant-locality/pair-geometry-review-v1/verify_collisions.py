#!/usr/bin/env python3
"""Rebuild collision grouping using full 200nt vectors, without outcomes."""
import collections
import csv
import hashlib
import itertools
import json
from pathlib import Path
from verify import BASES, COMP, feature_counts, replace

ROOT = Path('research/variant-locality')
SOURCE = ROOT / 'panel-scan-v1/opportunity-classes.csv'
SUMMARY = ROOT / 'collision-v1/summary.json'
OUT = ROOT / 'pair-geometry-review-v1/collision-receipt.json'


def main():
    words = sorted({min(w, w.translate(COMP)[::-1]) for w in
                    map(''.join, itertools.product(BASES, repeat=3))})
    index = {w: i for i, w in enumerate(words)}
    groups = collections.defaultdict(list)
    keys = set()
    with SOURCE.open(newline='') as handle:
        for row in csv.DictReader(handle):
            word, alt = row['reference_5mer'], row['alternate_base']
            key = (word, alt)
            assert len(word) == 5 and word[2] != alt and key not in keys
            keys.add(key)
            ref = 'TGCA' * 24 + word + ('ACGT' * 25)[:99]
            assert len(ref) == 200
            mutant = replace(ref, 98, alt)
            before, after = feature_counts(ref, index), feature_counts(mutant, index)
            delta = tuple(b-a for a, b in zip(before, after))
            assert delta == tuple(int(row[w+'_count_delta']) for w in words) + (int(row['GC_numerator_delta']),)
            orbit = min(key, (word.translate(COMP)[::-1], alt.translate(COMP)))
            groups[delta].append((orbit, int(row['panel_slot_count'])))
    assert len(keys) == 3072
    colliding = [group for group in groups.values() if len({x[0] for x in group}) > 1]
    nonzero = [group for vector, group in groups.items()
               if any(vector) and len({x[0] for x in group}) > 1]
    summary = json.loads(SUMMARY.read_text())
    assert len(groups) == summary['exact_delta_vectors']
    for label, group in [('beyond_rc_collision', colliding),
                         ('nonzero_beyond_rc_collision', nonzero)]:
        assert len(group) == summary[label+'_vector_count']
        assert sum(map(len, group)) == summary[label+'_directed_classes']
        assert sum(count for g in group for _, count in g) == summary[label+'_panel_slots']
    assert sum(count for g in groups.values() for _, count in g) == summary['panel_slot_denominator']
    receipt = {'method': 'Recompute counts from two full 200nt sequences per edit class and independently group by exact integer delta and RC orbit; source panel slot counts inherited',
               'input_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (SOURCE, SUMMARY)},
               'classes_checked': len(keys), 'exact_vectors': len(groups),
               'collision_vectors': len(colliding), 'nonzero_collision_vectors': len(nonzero),
               'all_checks_passed': True, 'outcomes_read': False}
    OUT.write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    print(json.dumps(receipt, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
