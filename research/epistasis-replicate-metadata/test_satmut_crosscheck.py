from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import fetch_satmut_member as fetch
import satmut_crosscheck as audit


class SatmutSequenceControls(unittest.TestCase):
    def test_long_names_follow_center_and_allele_order(self):
        group = {'variant1': 'chr1:100:A:C', 'variant2': 'chr1:103:G:T',
                 'center': 'var1', 'window': 'left'}
        self.assertEqual(audit.long_name(group, 'AR'),
                         '1:100:A:C:A:wL:Alt_1:103:G:T:R:m0')
        group.update(center='var2', window='right')
        self.assertEqual(audit.long_name(group, 'AR'),
                         '1:103:G:T:R:wR:Alt_1:100:A:C:A:m0')

    def test_full_sequence_comparison_retains_mismatches(self):
        a = 'A' * 200
        self.assertEqual(audit.base_differences(a, a), [])
        b = a[:20] + 'C' + a[21:]
        self.assertEqual(audit.base_differences(a, b),
                         [{'offset_zero_based': 20, 's18': 'A', 'satmut': 'C'}])
        with self.assertRaises(ValueError):
            audit.base_differences(a, a[:-1])
        with self.assertRaises(ValueError):
            audit.base_differences(a, 'N' + a[1:])

    def test_range_bounds_and_exclusive_output(self):
        reader = fetch.RangeReader()
        self.assertEqual(reader.seek(-1, 2), fetch.ARCHIVE_BYTES - 1)
        self.assertEqual(reader.read(0), b'')
        with self.assertRaises(ValueError):
            reader.seek(-1)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'member.gz'
            path.write_bytes(b'existing')
            with self.assertRaises(FileExistsError):
                fetch.write_new(path, b'new')
            self.assertEqual(path.read_bytes(), b'existing')

    def test_range_response_must_match_exact_request(self):
        class Response:
            def __init__(self, status, content_range, data):
                self.status = status
                self.headers = {'Content-Range': content_range}
                self.data = data
            def __enter__(self):
                return self
            def __exit__(self, *_args):
                pass
            def read(self, _count):
                return self.data

        good = f'bytes 0-2/{fetch.ARCHIVE_BYTES}'
        with patch.object(fetch.urllib.request, 'urlopen', return_value=Response(206, good, b'abc')):
            self.assertEqual(fetch.RangeReader().read(3), b'abc')
        for response in (Response(200, good, b'abc'),
                         Response(206, f'bytes 1-3/{fetch.ARCHIVE_BYTES}', b'abc'),
                         Response(206, good, b'ab')):
            with patch.object(fetch.urllib.request, 'urlopen', return_value=response):
                with self.assertRaises(ValueError):
                    fetch.RangeReader().read(3)


if __name__ == '__main__':
    unittest.main()
