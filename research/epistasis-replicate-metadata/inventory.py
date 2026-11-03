#!/usr/bin/env python3
"""Exploratory OL27 lookup-name inventory; does not read reporter outcomes.

The generated names are hypotheses, not verified S18-to-lookup identities.
See README.md before using any result as an analysis cohort.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / 'research/epistasis-source-audit/s18-design-audit.json'
DESIGN_SHA = '70dc4a6bd1813ea879b9699565cf75e8e9c21bff9ea1af5a46c25c1ad4686b0f'
FEATURES = ROOT / 'research/epistasis-feature-collisions-v1/results/ledger.jsonl'
FEATURES_SHA = '53c2751cd093da20a2c1bad0cd602e8ed173b7f164456c60f9938a7e1efa9ce7'
LOOKUP_SHA = '58f0c48c6d0a21b4a558794f31dcf38ccd6e8c766d56450512a28a27de059a15'
LOOKUP_MD5 = 'c18ba6775c64fbef1da73cede0d29dea'
LOOKUP_UNCOMPRESSED_BYTES = 1039687724


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def candidate_names(group):
    """Guess strings from a visible naming pattern; do not assert identity."""
    first = group['variant1'].removeprefix('chr').split(':')
    second = group['variant2'].removeprefix('chr').split(':')
    if len(first) != 4 or len(second) != 4:
        raise ValueError('unexpected variant identifier')
    if group['center'] not in ('var1', 'var2'):
        raise ValueError('unexpected center label')
    windows = {'left': 'L', 'middle': 'C', 'right': 'R'}
    if group['window'] not in windows:
        raise ValueError('unexpected window label')
    result = {}
    for label in ('RR', 'AR', 'RA', 'AA'):
        if group['center'] == 'var1':
            centered, other = (*first, label[0]), (second[0], second[1], label[1])
        else:
            centered, other = (*second, label[1]), (first[0], first[1], label[0])
        result[label] = ':'.join(centered) + ':w' + windows[group['window']] + ':Alt_' + '-'.join(other)
    return result


def lookup_fields(line):
    """Read only the oligo-name field and overall mapping flag."""
    fields = line.rstrip(b'\r\n').split(b'\t')
    if len(fields) != 10:
        raise ValueError('lookup row is not ten-field TSV')
    if fields[4] not in (b'0', b'1', b'2'):
        raise ValueError('unexpected overall mapping flag')
    names = fields[1].decode('utf-8').split(',')
    if not names or any(not name for name in names):
        raise ValueError('empty oligo name')
    return names, int(fields[4])


def key_for(group):
    return tuple(group[x] for x in ('variant1', 'variant2', 'center', 'window'))


def feature_classes(path, groups):
    classes = {}
    allowed = {'no_collision', 'first_blind_only', 'second_blind_only',
               'same_nonzero', 'opposite_nonzero', 'both_blind'}
    with path.open(encoding='utf-8') as source:
        for line in source:
            entry = json.loads(line)
            key = key_for(entry)
            if key in classes or entry['classification'] not in allowed:
                raise ValueError('duplicate key or unexpected feature classification')
            classes[key] = entry['classification']
    if set(classes) != {key_for(x) for x in groups}:
        raise ValueError('feature ledger and design groups differ')
    return classes


def overlap(groups, seen, classifications):
    by_count = Counter()
    collision_by_count, collision_full_classes = Counter(), Counter()
    complete_pairs, collision_full_pairs = set(), set()
    for group in groups:
        n = sum(name in seen for name in candidate_names(group).values())
        by_count[n] += 1
        if n == 4:
            complete_pairs.add((group['variant1'], group['variant2']))
        classification = classifications[key_for(group)]
        if classification != 'no_collision':
            collision_by_count[n] += 1
            if n == 4:
                collision_full_pairs.add((group['variant1'], group['variant2']))
                collision_full_classes[classification] += 1
    return {'windows_by_number_of_four_candidate_strings_present':
            {str(n): by_count[n] for n in range(5)},
            'pair_ids_with_at_least_one_four_string_window': len(complete_pairs),
            'collision_windows_by_number_of_four_candidate_strings_present':
            {str(n): collision_by_count[n] for n in range(5)},
            'collision_pair_ids_with_four_string_window': len(collision_full_pairs),
            'collision_classifications_for_four_string_windows': dict(sorted(collision_full_classes.items()))}


def write_report(path, result):
    # The check before a potentially long scan is only a convenience; this
    # exclusive create also protects a file created while the scan runs.
    with path.open('x', encoding='utf-8', newline='\n') as out:
        out.write(json.dumps(result, indent=2, sort_keys=True) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lookup', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True, help='new report path')
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('output already exists')
    if sha256(DESIGN) != DESIGN_SHA or sha256(FEATURES) != FEATURES_SHA or sha256(args.lookup) != LOOKUP_SHA:
        raise ValueError('design, feature ledger, or compressed lookup hash mismatch')
    groups = json.loads(DESIGN.read_text())['eligible_groups']
    classifications = feature_classes(FEATURES, groups)
    candidates = defaultdict(set)
    for group in groups:
        key = key_for(group)
        for label, name in candidate_names(group).items():
            candidates[name].add((key, label))
    all_seen, flag0_single_seen = set(), set()
    counts = Counter()
    digest = hashlib.md5()
    size = 0
    with gzip.open(args.lookup, 'rb') as stream:
        for line in stream:
            digest.update(line)
            size += len(line)
            names, flag = lookup_fields(line)
            counts['rows'] += 1
            counts[f'flag_{flag}_rows'] += 1
            if len(names) > 1:
                counts['multiple_name_rows'] += 1
            for name in names:
                if name in candidates:
                    all_seen.add(name)
                    if flag == 0 and len(names) == 1:
                        flag0_single_seen.add(name)
    if digest.hexdigest() != LOOKUP_MD5 or size != LOOKUP_UNCOMPRESSED_BYTES:
        raise ValueError('lookup uncompressed MD5 or size mismatch')
    result = {
        'schema_version': 1,
        'scope': 'exploratory string overlap only; not a verified S18-to-lookup crosswalk or outcome cohort',
        'design_sha256': DESIGN_SHA,
        'feature_ledger_sha256': FEATURES_SHA,
        'lookup_gzip_sha256': LOOKUP_SHA,
        'lookup_uncompressed_md5': LOOKUP_MD5,
        'lookup_uncompressed_bytes': size,
        'implementation_sha256': sha256(Path(__file__)),
        'design_windows': len(groups),
        'distinct_candidate_strings': len(candidates),
        'candidate_strings_reused_by_multiple_window_labels': sum(len(v) > 1 for v in candidates.values()),
        'lookup_row_counts': dict(sorted(counts.items())),
        'seen_under_any_flag': overlap(groups, all_seen, classifications),
        'seen_on_a_flag_0_single_name_row': overlap(groups, flag0_single_seen, classifications),
        'exposure': 'All lookup rows and QC/barcode fields instantiated as bytes; only oligo names and overall flags interpreted. No RNA or reporter-outcome table read.',
    }
    write_report(args.out, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
