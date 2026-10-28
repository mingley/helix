#!/usr/bin/env python3
"""Exhaustively group panel edit classes by exact feature difference, outcome-free."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = ROOT / 'research/variant-locality/panel-scan-v1/opportunity-classes.csv'
EXPECTED_INPUT = 'd7e10b3ae852796de1511b133a591ba698abb21c8bcbe180a2f30c2bbda91cea'
BASES = 'ACGT'
COMPLEMENT = str.maketrans('ACGT', 'TGCA')

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]

def features(sequence, kmers):
    counts = collections.Counter(min(sequence[i:i+3], rc(sequence[i:i+3])) for i in range(3))
    return tuple(counts[k] for k in kmers) + (sum(c in 'GC' for c in sequence),)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--out', type=Path, required=True, help='new or empty output directory')
    args = parser.parse_args()
    if digest(args.input) != EXPECTED_INPUT:
        raise ValueError('input hash mismatch')
    args.out.mkdir(parents=True, exist_ok=True)
    if any(args.out.iterdir()):
        raise ValueError('output directory must be empty')
    with args.input.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    kmers = sorted({min(a+b+c, rc(a+b+c)) for a in BASES for b in BASES for c in BASES})
    assert len(kmers) == 32 and len(rows) == 3072
    groups = collections.defaultdict(list)
    lookup = {}
    for row in rows:
        seq, alt = row['reference_5mer'], row['alternate_base']
        assert len(seq) == 5 and set(seq) <= set(BASES) and alt in BASES and seq[2] != alt
        assert row['reference_base'] == seq[2]
        mutant = seq[:2] + alt + seq[3:]
        delta = tuple(b-a for a, b in zip(features(seq, kmers), features(mutant, kmers)))
        stored = tuple(int(row[k+'_count_delta']) for k in kmers) + (int(row['GC_numerator_delta']),)
        assert delta == stored
        identity = (seq, alt)
        reverse = (rc(seq), rc(alt))
        swapped = (mutant, seq[2])
        assert (row['RC_reference_5mer'], row['RC_alternate_base']) == reverse
        assert (row['swapped_reference_5mer'], row['swapped_alternate_base']) == swapped
        assert identity not in lookup
        entry = {'reference': seq, 'alternate': alt, 'mutant': mutant,
                 'rc_orbit': min(identity, reverse), 'slot_count': int(row['panel_slot_count']),
                 'delta': delta}
        assert entry['slot_count'] > 0
        lookup[identity] = entry
        groups[delta].append(entry)
    for (seq, alt), entry in lookup.items():
        assert lookup[(rc(seq), rc(alt))]['delta'] == entry['delta']
        assert lookup[(entry['mutant'], seq[2])]['delta'] == tuple(-x for x in entry['delta'])
    assert len(groups) == 1405 and sum(e['slot_count'] for e in lookup.values()) == 32502288
    summaries = []
    for delta, entries in sorted(groups.items()):
        orbits = sorted(set(e['rc_orbit'] for e in entries))
        summaries.append({'delta': list(delta), 'is_zero': not any(delta),
                          'directed_classes': len(entries), 'rc_orbits': len(orbits),
                          'panel_slots': sum(e['slot_count'] for e in entries),
                          'collides_beyond_rc': len(orbits) > 1,
                          'representative_edits': ['%s>%s' % orbit for orbit in orbits]})
    collisions = [s for s in summaries if s['collides_beyond_rc']]
    nonzero = [s for s in collisions if not s['is_zero']]
    histogram = collections.Counter(s['rc_orbits'] for s in summaries)
    summary = {
      'input_sha256': EXPECTED_INPUT,
      'representation': '32 canonical overlapping trimer integer count deltas plus GC count delta; equal 200nt denominators 198 and 200 preserve equality',
      'directed_edit': 'ordered reference 5mer and alternate center base; mutant 5mer differs only at center',
      'rc_orbit': 'identify each directed edit with its reverse complement; do not merge an edit with its allele swap unless their delta vectors happen to coincide',
      'allele_swap': 'reverse the direction; the integer delta negates, so a swap is not generally equal and is never separately quotiented',
      'source_records': 55276, 'interior_positions_per_record': 196, 'alternates_per_position': 3,
      'panel_slot_denominator': sum(e['slot_count'] for e in lookup.values()),
      'directed_classes': len(lookup), 'exact_delta_vectors': len(groups),
      'rc_orbits': len({e['rc_orbit'] for e in lookup.values()}),
      'all_vector_orbit_count_histogram': dict(sorted(histogram.items())),
      'beyond_rc_collision_vector_count': len(collisions),
      'beyond_rc_collision_directed_classes': sum(s['directed_classes'] for s in collisions),
      'beyond_rc_collision_panel_slots': sum(s['panel_slots'] for s in collisions),
      'nonzero_beyond_rc_collision_vector_count': len(nonzero),
      'nonzero_beyond_rc_collision_directed_classes': sum(s['directed_classes'] for s in nonzero),
      'nonzero_beyond_rc_collision_panel_slots': sum(s['panel_slots'] for s in nonzero),
      'maximum_rc_orbits_per_vector': max(histogram),
      'limitations': ['A slot is an edit opportunity in an assay record, not a unique genomic variant, population frequency, mutation probability, or measured effect.',
                      'Equality holds for any fixed linear predictor on these features, not arbitrary nonlinear models or other feature sets.',
                      'These calculations read no assay outcomes and make no claim of biological equivalence.']}
    (args.out/'summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n')
    (args.out/'vectors.json').write_text(json.dumps(summaries, indent=2, sort_keys=True)+'\n')
    with (args.out/'collision-classes.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['reference_5mer','alternate_base','mutant_5mer','rc_orbit_reference','rc_orbit_alternate','panel_slot_count','zero_delta','rc_orbits_in_vector',*['delta_'+k for k in kmers],'delta_GC'])
        for delta, entries in sorted(groups.items()):
            n_orbits = len({e['rc_orbit'] for e in entries})
            if n_orbits <= 1:
                continue
            for e in sorted(entries, key=lambda x:(x['reference'],x['alternate'])):
                writer.writerow([e['reference'],e['alternate'],e['mutant'],*e['rc_orbit'],e['slot_count'],int(not any(delta)),n_orbits,*delta])
    manifest = {'input': {'sha256': digest(args.input), 'repository_path': 'research/variant-locality/panel-scan-v1/opportunity-classes.csv'},
                'implementation': {'sha256': digest(Path(__file__)), 'repository_path': 'research/variant-locality/collision-v1/analyze.py'},
                'outputs': {name: {'sha256': digest(args.out/name), 'bytes': (args.out/name).stat().st_size} for name in ['summary.json','vectors.json','collision-classes.csv']}}
    (args.out/'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if 'collision_' in k},sort_keys=True))

if __name__ == '__main__':
    main()
