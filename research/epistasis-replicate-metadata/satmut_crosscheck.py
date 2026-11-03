#!/usr/bin/env python3
"""Compare selected satmut background DNA with S18, without response values.

Long identifiers and sequences are compared directly. Short lookup identifiers
are generated hypotheses, not a verified original-design ID crosswalk.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import inventory

ROOT = Path(__file__).resolve().parents[2]
DESIGN_READER = ROOT / 'research/epistasis-reference-audit/audit.py'
READER_SHA = '66413f09afe3ad4ecb71ff8667cf8a95403ea9e4aa2a308585b3f73888ac8c4f'
S18_SHA = '92929b4cef82797237c4c2be770feb744daa97255ac6dbfe3d53d5a34c996fbb'
SATMUT_SHA = 'dd0080c8d451710e766e6b2dadb35e952da680183d47eb5775adc7fbaf91dea0'
HEADERS = ['sequence', 'mut_pos', 'mut_base', 'final_id', 'final_lsid']


def long_name(group, label):
    first = group['variant1'].removeprefix('chr').split(':')
    second = group['variant2'].removeprefix('chr').split(':')
    if len(first) != 4 or len(second) != 4 or label not in ('RR', 'AR', 'RA', 'AA'):
        raise ValueError('unexpected variant or genotype')
    if group['center'] == 'var1':
        center, other = (*first, label[0]), (*second, label[1])
    elif group['center'] == 'var2':
        center, other = (*second, label[1]), (*first, label[0])
    else:
        raise ValueError('unexpected center')
    windows = {'left': 'L', 'middle': 'C', 'right': 'R'}
    if group['window'] not in windows:
        raise ValueError('unexpected window')
    return ':'.join(center) + ':w' + windows[group['window']] + ':Alt_' + ':'.join(other) + ':m0'


def base_differences(a, b):
    if len(a) != 200 or len(b) != 200 or (set(a) | set(b)) - set('ACGT'):
        raise ValueError('expected two literal 200-base sequences')
    return [{'offset_zero_based': i, 's18': x, 'satmut': y}
            for i, (x, y) in enumerate(zip(a, b)) if x != y]


def design_reader():
    if inventory.sha256(DESIGN_READER) != READER_SHA:
        raise ValueError('S18 reader SHA-256 mismatch')
    spec = importlib.util.spec_from_file_location('s18_ref_reader', DESIGN_READER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def collect_lookup(lookup, names):
    if inventory.sha256(lookup) != inventory.LOOKUP_SHA:
        raise ValueError('lookup SHA-256 mismatch')
    any_seen, flag0 = set(), set()
    with gzip.open(lookup, 'rb') as source:
        for line in source:
            row_names, flag = inventory.lookup_fields(line)
            for name in row_names:
                if name in names:
                    any_seen.add(name)
                    if flag == 0 and len(row_names) == 1:
                        flag0.add(name)
    return any_seen, flag0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--satmut', type=Path, required=True, help='design-sequence member only')
    parser.add_argument('--lookup', type=Path, required=True, help='hash-bound OL27 mapping gzip')
    parser.add_argument('--out', type=Path, required=True, help='new output directory')
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('output already exists')
    reader = design_reader()
    bindings = [(args.satmut, SATMUT_SHA), (reader.S18, S18_SHA),
                (inventory.DESIGN, inventory.DESIGN_SHA), (inventory.FEATURES, inventory.FEATURES_SHA)]
    if any(inventory.sha256(path) != expected for path, expected in bindings):
        raise ValueError('input hash mismatch')
    groups = json.loads(inventory.DESIGN.read_text())['eligible_groups']
    classifications = inventory.feature_classes(inventory.FEATURES, groups)
    sequences = reader.source_sequences(groups, reader.S18)
    targets = defaultdict(list)
    for group, seqs in zip(groups, sequences):
        names = inventory.candidate_names(group)
        for label in ('RR', 'AR', 'RA', 'AA'):
            targets[long_name(group, label)].append((group, label, seqs[label], names[label]))

    rows, long_counts, matched_groups = [], Counter(), {}
    total_rows, background_rows = 0, 0
    with gzip.open(args.satmut, 'rt', encoding='utf-8', newline='') as source:
        table = csv.DictReader(source, delimiter='\t')
        if table.fieldnames != HEADERS:
            raise ValueError('unexpected satmut sequence header')
        for row_number, row in enumerate(table, 2):
            total_rows += 1
            name = row['final_id']
            if not name.endswith(':m0') or ':Alt_' not in name:
                continue
            background_rows += 1
            long_counts[name] += 1
            literal_200 = len(row['sequence']) == 200 and not (set(row['sequence']) - set('ACGT'))
            matches = []
            for group, label, sequence, short_name in targets.get(name, []):
                key = inventory.key_for(group)
                matched_groups[key] = group
                diffs = base_differences(sequence, row['sequence']) if literal_200 else None
                matches.append({'key': list(key), 'diplotype': label, 's18_row': group['rows'][label],
                                's18_sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
                                'hypothesized_short_name': short_name,
                                'feature_classification': classifications[key], 'differences': diffs,
                                'sequence_comparison': 'compared' if literal_200 else 'not_literal_200nt_acgt'})
            rows.append({'satmut_row': row_number, 'long_id': name,
                         'satmut_sequence_sha256': hashlib.sha256(row['sequence'].encode()).hexdigest(),
                         'satmut_sequence_length': len(row['sequence']),
                         'satmut_sequence_literal_200nt_acgt': literal_200,
                         'matches': matches})

    needed = {name for group in matched_groups.values() for name in inventory.candidate_names(group).values()}
    any_seen, flag0_seen = collect_lookup(args.lookup, needed)
    for row in rows:
        for match in row['matches']:
            key = tuple(match['key'])
            group = matched_groups[key]
            candidates = set(inventory.candidate_names(group).values())
            match['hypothesized_short_name_seen_any_flag'] = match['hypothesized_short_name'] in any_seen
            match['hypothesized_short_name_seen_flag0_single'] = match['hypothesized_short_name'] in flag0_seen
            match['all_four_hypothesized_names_seen_flag0_single'] = candidates <= flag0_seen

    match_rows = [row for row in rows if len(row['matches']) == 1]
    unmatched = sum(not row['matches'] for row in rows)
    ambiguous = sum(len(row['matches']) > 1 for row in rows)
    mismatched = sum(match['differences'] is not None and bool(match['differences'])
                     for row in rows for match in row['matches'])
    invalid = sum(match['differences'] is None for row in rows for match in row['matches'])
    complete_collision = {tuple(match['key']) for row in rows for match in row['matches']
                          if match['feature_classification'] != 'no_collision' and
                          match['all_four_hypothesized_names_seen_flag0_single']}
    summary = {
        'scope': 'selected/explored satmut design DNA only; no biological response or RNA count values read',
        'source_satmut_gzip_sha256': SATMUT_SHA,
        's18_sha256': S18_SHA, 'design_sha256': inventory.DESIGN_SHA,
        'feature_ledger_sha256': inventory.FEATURES_SHA, 'lookup_sha256': inventory.LOOKUP_SHA,
        's18_reader_sha256': READER_SHA, 'inventory_implementation_sha256': inventory.sha256(Path(inventory.__file__)),
        'implementation_sha256': inventory.sha256(Path(__file__)),
        'satmut_design_rows': total_rows, 'satmut_pair_m0_background_rows': background_rows,
        'duplicate_satmut_pair_m0_long_identifiers': sum(n > 1 for n in long_counts.values()),
        'pair_background_rows_matching_one_s18_long_candidate': len(match_rows),
        'pair_background_rows_with_no_s18_long_candidate': unmatched,
        'pair_background_rows_with_multiple_s18_long_candidates': ambiguous,
        'long_candidate_sequence_mismatches': mismatched,
        'long_candidate_sequences_not_literal_200nt_acgt': invalid,
        'pair_background_rows_not_literal_200nt_acgt': sum(not row['satmut_sequence_literal_200nt_acgt'] for row in rows),
        'distinct_matched_s18_windows': len(matched_groups),
        'matching_rows_whose_hypothesized_short_name_is_seen_any_flag':
            sum(match['hypothesized_short_name_seen_any_flag'] for row in rows for match in row['matches']),
        'matching_rows_whose_hypothesized_short_name_is_seen_flag0_single':
            sum(match['hypothesized_short_name_seen_flag0_single'] for row in rows for match in row['matches']),
        'collision_windows_with_satmut_match_and_all_four_hypothesized_short_names_seen_flag0_single': len(complete_collision),
        'short_identifier_to_sequence_contract_verified': False,
        'prior_outcome_based_selection_of_satmut_design': True,
    }
    args.out.mkdir(parents=True, exist_ok=False)
    ledger = args.out / 'pair-backgrounds.jsonl'
    with ledger.open('x', encoding='utf-8', newline='\n') as out:
        for row in rows:
            out.write(json.dumps(row, sort_keys=True, separators=(',', ':')) + '\n')
    summary['ledger_sha256'] = inventory.sha256(ledger)
    with (args.out / 'summary.json').open('x', encoding='utf-8', newline='\n') as out:
        out.write(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
