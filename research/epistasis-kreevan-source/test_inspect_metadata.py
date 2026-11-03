import gzip
import io
import random
import unittest

from inspect_metadata import DATA_PREFIX, decode_first_oligo, parse_samples


def sample_metadata():
    rows = []
    number = 1
    for cell, prefix in (('K562', 'K'), ('GM12878', 'G')):
        for replicate in range(1, 8):
            for material in ('cDNA', 'pDNA'):
                rows.extend((f'^SAMPLE = GSM{number}',
                             f'!Sample_title = {cell} replicate {replicate} {material} MPRA',
                             '!Sample_description = Column name in TableS2_MPRA_BarcodeCount.txt: '
                             f'{prefix}{replicate}_{material}',
                             f'!Sample_relation = SRA: https://example.org/SRX{number}',
                             '!Sample_data_processing = Assembly: hg19'))
                number += 1
    return '\n'.join(rows) + '\n'


class SampleControls(unittest.TestCase):
    def test_all_sample_pairs_and_assembly(self):
        samples = parse_samples(sample_metadata())
        self.assertEqual(len(samples), 28)
        self.assertEqual({x['column'] for x in samples if x['cell'] == 'K562'},
                         {f'K{i}_{m}' for i in range(1, 8) for m in ('cDNA', 'pDNA')})
        self.assertEqual({x['assembly'] for x in samples}, {'hg19'})

    def test_missing_mispairing_and_duplicate_sources_fail(self):
        original = sample_metadata()
        with self.assertRaises(ValueError):
            parse_samples(original.replace('K3_cDNA', 'K3_pDNA'))
        with self.assertRaises(ValueError):
            parse_samples(original.replace('Assembly: hg19', 'Assembly: hg38', 1))
        with self.assertRaises(ValueError):
            parse_samples(original.replace('SRX2\n', 'SRX1\n', 1))


class PrefixControls(unittest.TestCase):
    def make_input(self, oligo='opaque_1_A/T_anc/der', bad_header=False):
        columns = list(DATA_PREFIX) + [f'{p}{i}_{m}' for p in ('G', 'K')
                                       for i in range(1, 8) for m in ('cDNA', 'pDNA')]
        if bad_header:
            columns[1] = 'sequence'
        # Data after Oligo represents fields that the scanner must not decode.
        rng = random.Random(7)
        sentinel = ''.join(str(rng.randrange(10)) for _ in range(12000))
        plain = ('\t'.join(columns) + '\n' + 'ACGTACGTACGTACG\t' + oligo + '\t' + sentinel + '\n')
        return gzip.compress(plain.encode('ascii'), mtime=0), sentinel

    def test_first_oligo_stops_before_count_fields(self):
        compressed, sentinel = self.make_input()
        stream = io.BytesIO(compressed)
        result = decode_first_oligo(stream)
        self.assertEqual(result['oligo_field'], 'opaque_1_A/T_anc/der')
        self.assertFalse(result['oligo_is_literal_acgt'])
        self.assertEqual(result['decoded_bytes_after_second_tab'], 0)
        self.assertLess(stream.tell(), len(compressed))
        self.assertNotIn(sentinel, str(result))

    def test_literal_dna_field_and_bad_header(self):
        compressed, _ = self.make_input('ACGTACGT')
        self.assertTrue(decode_first_oligo(io.BytesIO(compressed))['oligo_is_literal_acgt'])
        bad, _ = self.make_input(bad_header=True)
        with self.assertRaisesRegex(ValueError, 'unexpected count header'):
            decode_first_oligo(io.BytesIO(bad))

    def test_truncated_prefix_fails_closed(self):
        compressed, _ = self.make_input()
        with self.assertRaises(ValueError):
            decode_first_oligo(io.BytesIO(compressed[:20]))


if __name__ == '__main__':
    unittest.main()
