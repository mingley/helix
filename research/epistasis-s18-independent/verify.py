#!/usr/bin/env python3
"""Independently verify S18 pair-design audit using XLSX XML, not openpyxl.

This reads only the official S18 design workbook, never any outcome workbook.
Run from the repository root with the source path as the only argument.
"""
import collections
import hashlib
import itertools
import json
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
EXPECTED = '92929b4cef82797237c4c2be770feb744daa97255ac6dbfe3d53d5a34c996fbb'
REPORT = Path('research/epistasis-source-audit/s18-design-audit.json')
RECEIPT = Path('research/epistasis-s18-independent/receipt.json')
VARIANT = re.compile(r'chr(?:[1-9]|1\d|2[0-2]|X|Y):[0-9]+:[ACGT]:[ACGT]')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def extract(source):
    assert sha(source) == EXPECTED
    with zipfile.ZipFile(source) as z:
        workbook = ET.fromstring(z.read('xl/workbook.xml'))
        sheets = workbook.findall('.//' + NS + 'sheet')
        assert len(sheets) == 1 and sheets[0].attrib['name'] == 'S18.Variant Pair Oligos'
        strings = []
        with z.open('xl/sharedStrings.xml') as data:
            for _, elem in ET.iterparse(data, events=('end',)):
                if elem.tag == NS + 'si':
                    strings.append(''.join(t.text or '' for t in elem.iter(NS + 't')))
                    elem.clear()
        with z.open('xl/worksheets/sheet1.xml') as data:
            for _, elem in ET.iterparse(data, events=('end',)):
                if elem.tag != NS + 'row':
                    continue
                number = int(elem.attrib['r'])
                cells = {}
                for cell in elem.findall(NS + 'c'):
                    match = re.fullmatch(r'([A-G])([0-9]+)', cell.attrib['r'])
                    assert match and int(match.group(2)) == number
                    assert cell.find(NS + 'f') is None
                    value = cell.find(NS + 'v')
                    assert value is not None and value.text is not None
                    typ = cell.attrib.get('t', 'n')
                    if typ == 's':
                        obj = strings[int(value.text)]
                    else:
                        assert typ == 'n' and re.fullmatch(r'\d+', value.text)
                        obj = int(value.text)
                    assert match.group(1) not in cells
                    cells[match.group(1)] = obj
                yield number, cells
                elem.clear()


def eligible(key, items):
    v1, v2, center, window = key
    if not all(isinstance(v, str) and VARIANT.fullmatch(v) for v in (v1, v2)):
        return None
    c1, n1, r1, a1 = v1.split(':')
    c2, n2, r2, a2 = v2.split(':')
    n1, n2 = int(n1), int(n2)
    d = abs(n1-n2)
    if c1 != c2 or d < 3 or r1 == a1 or r2 == a2:
        return None
    if center not in ('var1', 'var2') or window not in ('left', 'middle', 'right'):
        return None
    if len(items) != 4 or {x[2] for x in items} != {'RR', 'AR', 'RA', 'AA'}:
        return None
    if any(x[1] != d or not isinstance(x[3], str) or len(x[3]) != 200 or
           set(x[3])-set('ACGT') for x in items):
        return None
    seq = {item[2]: item[3] for item in items}
    # Enumerate the four different bases at each position and require exactly
    # two varying sites whose genotype patterns are AR and RA respectively.
    polymorphic = [i for i in range(200) if len({s[i] for s in seq.values()}) != 1]
    if len(polymorphic) != 2:
        return None
    patterns = {i: tuple(seq[label][i] for label in ('RR', 'AR', 'RA', 'AA'))
                for i in polymorphic}
    first = [i for i in polymorphic if patterns[i] == (r1, a1, r1, a1)]
    second = [i for i in polymorphic if patterns[i] == (r2, r2, a2, a2)]
    if len(first) != 1 or len(second) != 1 or first[0] == second[0]:
        return None
    if abs(first[0]-second[0]) != d:
        return None
    return {'variant1': v1, 'variant2': v2, 'center': center, 'window': window,
            'distance': d, 'rows': {item[2]: item[0] for item in items},
            'position1_zero_based': first[0], 'position2_zero_based': second[0]}


def main():
    assert len(sys.argv) == 2, 'provide S18 workbook path'
    source = Path(sys.argv[1])
    expected_headers = ('Variant 1', 'Variant 2', 'Distance', 'Center Variant',
                        'Window', 'Diplotype', 'Oligo')
    groups = collections.defaultdict(list)
    rowcount = 0
    labelcounts = collections.Counter()
    for rownum, cells in extract(source):
        if rownum == 1:
            assert cells['A'] == 'Supplementary Table 18 | Oligo sequences for variant pairs'
            continue
        if rownum == 2:
            assert tuple(cells[x] for x in 'ABCDEFG') == expected_headers
            continue
        assert set(cells) == set('ABCDEFG')
        key = tuple(cells[x] for x in ('A', 'B', 'D', 'E'))
        groups[key].append((rownum, cells['C'], cells['F'], cells['G']))
        labelcounts[cells['F']] += 1
        rowcount += 1
    accepted, rejected = {}, {}
    for key, items in groups.items():
        result = eligible(key, items)
        if result is None:
            rejected[key] = sorted(item[0] for item in items)
        else:
            accepted[key] = result
    report = json.loads(REPORT.read_text())
    expected_accepted = {(x['variant1'], x['variant2'], x['center'], x['window']): x
                         for x in report['eligible_groups']}
    expected_rejected = {tuple(x['key']): sorted(x['rows']) for x in report['ineligible_groups']}
    assert len(expected_accepted) == len(report['eligible_groups'])
    assert len(expected_rejected) == len(report['ineligible_groups'])
    assert accepted == expected_accepted
    assert rejected == expected_rejected
    assert rowcount == report['data_rows']
    assert len(groups) == report['design_groups']
    assert dict(labelcounts) == report['diplotype_rows']
    assert len({(k[0], k[1]) for k in groups}) == report['distinct_pair_identifiers']
    assert len({(k[0], k[1]) for k in accepted}) == report['eligible_distinct_pairs']
    distances = collections.Counter(str(x['distance']) for x in accepted.values())
    assert dict(distances) == report['eligible_distance_distribution']
    receipt = {'method': 'Separate stdlib streaming XLSX XML parser; all four sequences compared by polymorphic-column genotype patterns',
               'source_sha256': sha(source), 'primary_report_sha256': sha(REPORT),
               'rows_checked': rowcount, 'groups_checked': len(groups),
               'accepted_groups_exactly_matched': len(accepted),
               'rejected_group_keys_and_rows_exactly_matched': len(rejected),
               'accepted_distinct_pair_identifiers': report['eligible_distinct_pairs'],
               'accepted_group_details_checked': ['row numbers', 'sequence edit positions', 'distance', 'labels', 'pair and window identifiers'],
               'diplotype_row_counts_and_distance_distribution_match': True,
               'reject_reason_categories_independently_checked': False,
               'outcome_workbooks_read': False, 'all_checks_passed': True}
    RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
