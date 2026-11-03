import gzip
from pathlib import Path
import tempfile
import unittest

import inventory
import verify_metadata


class CandidateNames(unittest.TestCase):
    def setUp(self):
        self.group = {'variant1': 'chr1:100:A:C', 'variant2': 'chr1:103:G:T',
                      'center': 'var1', 'window': 'left'}

    def test_center_orientation_and_allele_labels(self):
        names = inventory.candidate_names(self.group)
        self.assertEqual(names['AR'], '1:100:A:C:A:wL:Alt_1-103-R')
        self.group['center'] = 'var2'
        self.group['window'] = 'middle'
        names = inventory.candidate_names(self.group)
        self.assertEqual(names['AR'], '1:103:G:T:R:wC:Alt_1-100-A')

    def test_reject_unexpected_window_or_malformed_lookup_row(self):
        self.group['window'] = 'unknown'
        with self.assertRaises(ValueError):
            inventory.candidate_names(self.group)
        with self.assertRaises(ValueError):
            inventory.lookup_fields(b'not\tenough\tfields\n')
        with self.assertRaises(ValueError):
            inventory.lookup_fields(b'a\tb\tc\td\t3\tf\tg\th\ti\tj\n')
        self.assertEqual(inventory.lookup_fields(b'barcode\tid1,id2\tx\tx\t1\tx\tx\tx\tx\tx\n'),
                         (['id1', 'id2'], 1))

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'report.json'
            path.write_text('existing report\n')
            with self.assertRaises(FileExistsError):
                inventory.write_report(path, {'replacement': True})
            self.assertEqual(path.read_text(), 'existing report\n')


class BoundedHeader(unittest.TestCase):
    def test_decoder_returns_only_first_line(self):
        compressed = gzip.compress(b'header\tcol\n12\t34\n')
        first, used = verify_metadata.decode_first_line(compressed)
        self.assertEqual(first, b'header\tcol\n')
        self.assertLessEqual(used, len(compressed))

    def test_incomplete_and_overlong_header_are_errors(self):
        with self.assertRaises(ValueError):
            verify_metadata.decode_first_line(gzip.compress(b'no newline'))
        with self.assertRaises(ValueError):
            verify_metadata.decode_first_line(gzip.compress(b'abcd\n'), limit=2)

    def test_retained_metadata_hashes_and_replicates(self):
        sources = Path(__file__).with_name('sources')
        self.assertEqual(len(verify_metadata.verify_sources(sources)['records']), 16)
        groups = verify_metadata.verify_sample_metadata(sources)
        self.assertEqual(groups['rna']['biological_replicate_labels'], [1, 2, 3, 4, 5])
        self.assertEqual(groups['plasmid']['biological_replicate_labels'], [1, 2, 3, 4, 5])


if __name__ == '__main__':
    unittest.main()
