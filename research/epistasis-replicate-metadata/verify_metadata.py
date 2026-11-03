#!/usr/bin/env python3
"""Verify retained GEO metadata and optionally replay a bounded count header.

No count-table data row is returned or interpreted by this program. The optional
range request buffers compressed bytes that may encode rows; see README.md.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import urllib.request
import zlib

HERE = Path(__file__).resolve().parent
MANIFEST_SHA = 'ecf807627442f4dec3d2000e758388d4393357a86dd0a8ee91402594d7d3c421'
RNA = [f'GSM825448{x}' for x in range(4, 9)]
DNA = [f'GSM82547{x}' for x in range(37, 42)]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def values(text, key):
    prefix = key + ' = '
    return [line[len(prefix):] for line in text.splitlines() if line.startswith(prefix)]


def decode_first_line(compressed, limit=4096):
    """Return only the first gzip line and number of compressed bytes fed."""
    dec = zlib.decompressobj(31)
    result = bytearray()
    for index, byte in enumerate(compressed):
        pending = bytes([byte])
        while pending:
            returned = dec.decompress(pending, 1)
            pending = dec.unconsumed_tail
            if returned:
                result.extend(returned)
                if returned == b'\n':
                    return bytes(result), index + 1
                if len(result) > limit:
                    raise ValueError('first line exceeds limit')
    raise ValueError('gzip prefix does not contain a complete first line')


def verify_sources(directory):
    data = (directory / 'manifest.json').read_bytes()
    if sha(data) != MANIFEST_SHA:
        raise ValueError('metadata manifest hash mismatch')
    manifest = json.loads(data)
    seen = set()
    for entry in manifest['records']:
        name = entry['file']
        if Path(name).name != name or name in seen:
            raise ValueError('invalid or duplicate metadata filename')
        seen.add(name)
        contents = (directory / name).read_bytes()
        if len(contents) != entry['saved_bytes'] or sha(contents) != entry['saved_sha256']:
            raise ValueError('metadata source hash or size mismatch: ' + name)
    if len(seen) != 16:
        raise ValueError('unexpected metadata record count')
    return manifest


def verify_sample_metadata(directory):
    for accession, expected in [('GSE266781', RNA), ('GSE266843', DNA)]:
        text = (directory / (accession + '.soft')).read_text()
        if sorted(values(text, '!Series_sample_id')) != sorted(expected):
            raise ValueError('GEO series sample membership mismatch')
    groups = {}
    for label, samples in [('rna', RNA), ('plasmid', DNA)]:
        reps, sra = [], set()
        for accession in samples:
            text = (directory / (accession + '.soft')).read_text()
            if values(text, '!Sample_geo_accession') != [accession]:
                raise ValueError('GEO sample accession mismatch')
            bio = values(text, '!Sample_description')
            biological = [x.removeprefix('biological replicate number: ') for x in bio
                          if x.startswith('biological replicate number: ')]
            technical = [x.removeprefix('technical replicate number: ') for x in bio
                         if x.startswith('technical replicate number: ')]
            relation = [x for x in values(text, '!Sample_relation') if x.startswith('SRA: ')]
            if len(biological) != 1 or technical != ['1'] or len(relation) != 1:
                raise ValueError('unexpected replicate or SRA metadata')
            reps.append(int(biological[0]))
            sra.add(relation[0])
        if sorted(reps) != [1, 2, 3, 4, 5] or len(sra) != 5:
            raise ValueError('unexpected replicate labels or duplicated SRA records')
        groups[label] = {'samples': samples, 'biological_replicate_labels': sorted(reps),
                         'technical_replicate_labels': [1], 'distinct_sra_records': len(sra)}
    return groups


def verify_workbook(path, evidence):
    if sha(path.read_bytes()) != evidence['s29_workbook_sha256']:
        raise ValueError('publisher workbook SHA-256 mismatch')
    import openpyxl
    book = openpyxl.load_workbook(path, read_only=True, data_only=False)
    try:
        sheet = book[evidence['s29_sheet']]
        row = next(sheet.iter_rows(min_row=evidence['s29_row_1based'],
                                   max_row=evidence['s29_row_1based'], max_col=6,
                                   values_only=True))
        if list(row) != evidence['s29_row']:
            raise ValueError('S29 metadata row mismatch')
    finally:
        book.close()


def verify_live_header(record):
    req = urllib.request.Request(record['url'], headers={
        'Range': f"bytes={record['range_start']}-{record['range_end_inclusive']}",
        'User-Agent': 'helix-research-metadata/1.0',
    })
    with urllib.request.urlopen(req, timeout=60) as response:
        expected = f"bytes 0-65535/{record['gzip_bytes']}"
        if response.status != 206 or response.headers.get('Content-Range') != expected:
            raise ValueError('server did not return the expected bounded range')
        raw = response.read(65537)
    if len(raw) != 65536 or sha(raw) != record['range_sha256']:
        raise ValueError('compressed prefix hash or size mismatch')
    header, used = decode_first_line(raw)
    if header != (record['header'] + '\n').encode() or used != record['compressed_bytes_fed_to_inflater']:
        raise ValueError('count table header mismatch')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', type=Path, help='optional hash-bound publisher S30 bundle')
    parser.add_argument('--live-header', action='store_true', help='fetch bounded count-table prefix')
    args = parser.parse_args()
    evidence = json.loads((HERE / 'evidence.json').read_text())
    manifest = verify_sources(HERE / 'sources')
    groups = verify_sample_metadata(HERE / 'sources')
    rna_series = (HERE / 'sources/GSE266781.soft').read_text()
    rna_readme = (HERE / 'sources/GSM8254484_README.txt').read_text()
    lookup_readme = (HERE / 'sources/GSM8254614_README.txt').read_text()
    if ('ENCFF970OLE' not in rna_series or 'ENCFF823YWM' not in rna_readme or
            evidence['count_table']['uncompressed_md5_per_geo_readme'] not in rna_readme or
            evidence['lookup_table']['uncompressed_md5_per_geo_readme'] not in lookup_readme):
        raise ValueError('source-file accession or checksum metadata mismatch')
    if args.workbook:
        verify_workbook(args.workbook, evidence)
    if args.live_header:
        verify_live_header(evidence['count_table'])
    print(json.dumps({'metadata_sources_verified': len(manifest['records']), 'groups': groups,
                      'workbook_row_checked': bool(args.workbook),
                      'live_count_header_checked': args.live_header,
                      'rna_outcome_rows_interpreted': False}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
