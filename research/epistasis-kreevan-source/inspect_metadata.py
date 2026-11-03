#!/usr/bin/env python3
"""Inspect Kreevan 2026 design metadata without reading assay count rows.

Inputs for the offline inspection are the separately downloaded article XML,
GEO series and GEO sample metadata, plus a header-only capture. No S2/S4/S8
data rows are parsed. Output is a new JSON file and never overwrites evidence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import urllib.request
import xml.etree.ElementTree as ET
import zlib

ARTICLE_SHA = "71c19fda7b746be06ee7e9959be3eb37f7e9516627a702b37ff5be20e6f29b66"
SERIES_SHA = "8891217610f90ab81392de85e7df7e3bab46883aabc620202f9e793a9ef21dd6"
SAMPLES_SHA = "4f7e3facc4c7354ce72db508f910768016c675a719bf764efe2dca04c5bfbeb3"
HEADER_SHA = "fbdc8b7202e1870700232a1c47192fcffad0bc6d7e4d185c193d12f4991a7502"
FIRST_OLIGO_SHA = "daf20bfcec572b9f2a6faaf07de855026b746c77643d23ab7a345d322232fd07"
COUNT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE349nnn/GSE349490/suppl/GSE349490_TableS2_MPRA_BarcodeCount.txt.gz"
DATA_PREFIX = ("Barcode", "Oligo", "Error", "CIGAR", "MD", "Aln_Start:Stop")


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def save_new(path, data):
    with path.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write('\n')


def capture_header(out):
    """Decompress only the gzip header line from a bounded HTTPS range.

    Socket/HTTP layers may buffer compressed bytes from the requested range;
    only the prefix through the first newline is fed to the decompressor. A
    decompressed byte after the header is treated as an error and never saved.
    """
    request = urllib.request.Request(COUNT_URL, headers={'Range': 'bytes=0-8191'})
    with urllib.request.urlopen(request, timeout=30) as response:
        content_range = response.headers.get('Content-Range', '')
        require(response.status == 206 and content_range.startswith('bytes 0-8191/'),
                'server did not honor the bounded range')
        decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
        decoded, compressed, pending = bytearray(), bytearray(), b''
        while not decoded.endswith(b'\n'):
            if pending:
                chunk = decoder.decompress(pending, 1)
                pending = decoder.unconsumed_tail
            else:
                chunk = decoder.decompress(b'', 1)
                if not chunk:
                    pending = response.read(1)
                    require(bool(pending) and len(compressed) < 8192, 'header missing in bounded range')
                    compressed.extend(pending)
                    continue
            decoded.extend(chunk)
            require(len(decoded) <= 4096, 'unexpectedly long header')
        fields = bytes(decoded).decode('ascii').rstrip('\r\n').split('\t')
        save_new(out, {'source_url': COUNT_URL, 'http_status': response.status,
                       'content_range': content_range, 'compressed_bytes_fed_to_decoder': len(compressed),
                       'compressed_prefix_sha256': hashlib.sha256(compressed).hexdigest(),
                       'decoded_bytes_after_header': 0, 'columns': fields})


def decode_first_oligo(stream):
    """Decode only the first record's barcode and Oligo fields, not counts.

    A one-byte decompressed-output limit makes the second tab a hard stop.
    The gzip decoder may hold compressed input internally, and the HTTP layer
    may buffer the requested range; no bytes beyond that tab are decompressed.
    """
    decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
    output, compressed, pending = bytearray(), bytearray(), b''
    header = None
    while True:
        if pending:
            chunk = decoder.decompress(pending, 1)
            pending = decoder.unconsumed_tail
        else:
            chunk = decoder.decompress(b'', 1)
            if not chunk:
                pending = stream.read(1)
                require(bool(pending) and len(compressed) < 8192, 'first oligo missing in range')
                compressed.extend(pending)
                continue
        if not chunk:
            continue
        output.extend(chunk)
        require(len(output) <= 4096, 'unexpectedly long header or first fields')
        if header is None and chunk == b'\n':
            header = bytes(output).decode('ascii').rstrip('\r\n').split('\t')
            require(tuple(header[:6]) == DATA_PREFIX and len(header) == 34, 'unexpected count header')
            output.clear()
        elif header is not None and output.count(b'\t') == 2:
            barcode, oligo, remainder = bytes(output).split(b'\t', 2)
            require(not remainder and barcode and oligo, 'malformed first data fields')
            text = oligo.decode('ascii')
            return {'oligo_field': text, 'oligo_length': len(oligo),
                    'oligo_is_literal_acgt': re.fullmatch(r'[ACGT]+', text) is not None,
                    'oligo_sha256': hashlib.sha256(oligo).hexdigest(),
                    'leading_barcode_length': len(barcode),
                    'compressed_bytes_fed_to_decoder': len(compressed),
                    'compressed_prefix_sha256': hashlib.sha256(compressed).hexdigest(),
                    'decoded_bytes_after_second_tab': 0}


def capture_first_oligo(out):
    request = urllib.request.Request(COUNT_URL, headers={'Range': 'bytes=0-8191'})
    with urllib.request.urlopen(request, timeout=30) as response:
        content_range = response.headers.get('Content-Range', '')
        require(response.status == 206 and content_range.startswith('bytes 0-8191/'),
                'server did not honor the bounded range')
        record = decode_first_oligo(response)
        record.update({'source_url': COUNT_URL, 'http_status': response.status,
                       'content_range': content_range})
        save_new(out, record)


def parse_samples(text):
    """Read only metadata fields identifying each replicate/material sample."""
    samples = []
    current = None
    for line in text.splitlines():
        if line.startswith('^SAMPLE = '):
            if current is not None:
                samples.append(current)
            current = {'accession': line.split(' = ', 1)[1]}
        elif current is not None and line.startswith('!Sample_title = '):
            require('title' not in current, 'duplicate sample title')
            current['title'] = line.split(' = ', 1)[1]
        elif current is not None and line.startswith('!Sample_description = Column name in TableS2_MPRA_BarcodeCount.txt: '):
            require('column' not in current, 'duplicate sample column')
            current['column'] = line.rsplit(': ', 1)[1]
        elif current is not None and line.startswith('!Sample_relation = SRA: '):
            require('sra' not in current, 'duplicate SRA relation')
            current['sra'] = line.split(' = SRA: ', 1)[1]
        elif current is not None and line == '!Sample_data_processing = Assembly: hg19':
            require('assembly' not in current, 'duplicate assembly metadata')
            current['assembly'] = 'hg19'
    if current is not None:
        samples.append(current)
    seen = set()
    for sample in samples:
        require(set(sample) == {'accession', 'title', 'column', 'sra', 'assembly'}, 'incomplete sample metadata')
        match = re.fullmatch(r'(K562|GM12878) replicate ([1-7]) (cDNA|pDNA) MPRA', sample['title'])
        require(match is not None, 'unexpected sample title')
        cell, replicate, material = match.groups()
        expected = ('K' if cell == 'K562' else 'G') + replicate + '_' + material
        require(sample['column'] == expected, 'sample title and column disagree')
        require(sample['column'] not in seen, 'duplicate sample column')
        seen.add(sample['column'])
        sample.update({'cell': cell, 'replicate': int(replicate), 'material': material})
    expected = {prefix + str(i) + '_' + material for prefix in ('K', 'G')
                for i in range(1, 8) for material in ('cDNA', 'pDNA')}
    require(len(samples) == 28 and seen == expected, 'missing or extra replicate/material samples')
    require(len({x['accession'] for x in samples}) == 28 and len({x['sra'] for x in samples}) == 28,
            'sample or SRA accession reused')
    return sorted(samples, key=lambda x: x['column'])


def inspect_methods(article_path):
    root = ET.parse(article_path).getroot()
    checks = {
        's4a1': ('up to 75 bp', 'ancestral', 'derived', 'Table S1'),
        's4a2': ('170 bp of genomic sequence', '15 bp adapter sequences on each end'),
        's4b1': ('90% of the cells were used for RNA extraction', 'two aliquots of 5%', '7 replicates'),
        's4c': ('Table S2', 'median ratios method', 'mean raw count of >400'),
    }
    for ident, phrases in checks.items():
        matches = [s for s in root.iter('sec') if s.get('id') == ident]
        require(len(matches) == 1, f'missing or duplicate methods section {ident}')
        text = ' '.join(' '.join(''.join(p.itertext()).split()) for p in matches[0].findall('p'))
        for phrase in phrases:
            require(phrase in text, f'methods binding changed in {ident}')
    license_node = root.find('.//article-meta/permissions/license')
    require(license_node is not None and license_node.get('{http://www.w3.org/1999/xlink}href') ==
            'http://creativecommons.org/licenses/by/4.0/', 'article license metadata differs')


def inspect(series_path, sample_path, article_path, header_path, first_oligo_path):
    for path, digest, name in ((series_path, SERIES_SHA, 'GEO series'),
                               (sample_path, SAMPLES_SHA, 'GEO samples'),
                               (article_path, ARTICLE_SHA, 'article XML'),
                               (header_path, HEADER_SHA, 'count header receipt'),
                               (first_oligo_path, FIRST_OLIGO_SHA, 'first oligo receipt')):
        require(sha256(path) == digest, f'{name} hash mismatch')
    inspect_methods(article_path)
    series_lines = series_path.read_text(encoding='utf-8').splitlines()
    require('^SERIES = GSE349490' in series_lines, 'unexpected GEO series')
    supplied = [s.split(' = ', 1)[1] for s in series_lines
                if s.startswith('!Series_supplementary_file = ')]
    expected_suffix = ('GSE349490_TableS2_MPRA_BarcodeCount.txt.gz',
                       'GSE349490_TableS4_k562_results.tsv.gz',
                       'GSE349490_TableS8_gm12878_results.txt.gz')
    require(len(supplied) == 3 and tuple(x.rsplit('/', 1)[-1] for x in supplied) == expected_suffix,
            'GEO supplementary file inventory differs')
    samples = parse_samples(sample_path.read_text(encoding='utf-8'))
    header = json.loads(header_path.read_text(encoding='utf-8'))
    require(header.get('source_url') == COUNT_URL and header.get('http_status') == 206 and
            header.get('content_range') == 'bytes 0-8191/188701018' and
            header.get('decoded_bytes_after_header') == 0, 'unexpected S2 header capture')
    columns = header.get('columns')
    require(isinstance(columns, list) and tuple(columns[:6]) == DATA_PREFIX and
            len(columns) == 34 and len(set(columns)) == 34 and
            set(columns[6:]) == {x['column'] for x in samples}, 'count header and sample metadata differ')
    first = json.loads(first_oligo_path.read_text(encoding='utf-8'))
    require(first.get('source_url') == COUNT_URL and first.get('http_status') == 206 and
            first.get('content_range') == 'bytes 0-8191/188701018' and
            first.get('decoded_bytes_after_second_tab') == 0 and
            first.get('oligo_is_literal_acgt') is False, 'unexpected first-oligo capture')
    return {'source': {'article_doi': '10.64898/2026.09.10.750584',
                       'article_sha256': ARTICLE_SHA, 'article_license': 'CC BY 4.0',
                       'geo_series': 'GSE349490', 'series_sha256': SERIES_SHA,
                       'samples_sha256': SAMPLES_SHA,
                       'count_header_capture_sha256': HEADER_SHA,
                       'first_oligo_capture_sha256': FIRST_OLIGO_SHA,
                       'inspection_implementation_sha256': sha256(Path(__file__).resolve())},
            'design_metadata': {'genomic_insert_length': 170, 'adapter_length_each_side': 15,
                                'maximum_pair_distance': 75, 'assembly': 'hg19 per GEO sample metadata',
                                'haplotype_labels': ['ancestral/ancestral', 'ancestral/derived',
                                                     'derived/ancestral', 'derived/derived'],
                                'design_sequence_table_cited': 'Table S1'},
            'geo_supplementary_files': supplied, 'count_header_columns': columns,
            'samples': samples, 'first_oligo_field_is_literal_dna': False,
            'count_or_response_values_decoded': False,
            'sequence_identity_established': False,
            'independent_biological_replicates_established': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    capture = commands.add_parser('capture-header')
    capture.add_argument('--out', required=True, type=Path)
    first = commands.add_parser('capture-first-oligo')
    first.add_argument('--out', required=True, type=Path)
    audit = commands.add_parser('inspect')
    for name in ('series', 'samples', 'article', 'header', 'first-oligo', 'out'):
        audit.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    if args.command == 'capture-header':
        capture_header(args.out)
    elif args.command == 'capture-first-oligo':
        capture_first_oligo(args.out)
    else:
        save_new(args.out, inspect(args.series, args.samples, args.article, args.header, args.first_oligo))


if __name__ == '__main__':
    main()
