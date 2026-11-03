import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from verify import eligible


def example():
    key = ('chr1:100:A:C', 'chr1:103:A:G', 'var1', 'middle')
    ref = 'A' * 200
    one = ref[:20] + 'C' + ref[21:]
    two = ref[:23] + 'G' + ref[24:]
    both = one[:23] + 'G' + one[24:]
    items = [(3, 3, 'RR', ref), (4, 3, 'AR', one),
             (5, 3, 'RA', two), (6, 3, 'AA', both)]
    return key, items


class SyntheticSquares(unittest.TestCase):
    def test_valid_and_label_permutation(self):
        key, items = example()
        result = eligible(key, list(reversed(items)))
        self.assertEqual((20, 23), (result['position1_zero_based'], result['position2_zero_based']))
        self.assertEqual({'RR': 3, 'AR': 4, 'RA': 5, 'AA': 6}, result['rows'])

    def test_missing_or_repeated_labels(self):
        key, items = example()
        self.assertIsNone(eligible(key, items[:3]))
        items[3] = (6, 3, 'RA', items[3][3])
        self.assertIsNone(eligible(key, items))

    def test_double_must_equal_combined_edits(self):
        key, items = example()
        number, distance, label, seq = items[3]
        items[3] = (number, distance, label, seq[:100] + 'T' + seq[101:])
        self.assertIsNone(eligible(key, items))
        items[3] = (number, distance, label, items[1][3])
        self.assertIsNone(eligible(key, items))

    def test_source_distance_and_literal_alleles_must_agree(self):
        key, items = example()
        self.assertIsNone(eligible(('chr1:100:A:G', *key[1:]), items))
        self.assertIsNone(eligible((key[0], 'chr1:102:A:G', *key[2:]), items))
        self.assertIsNone(eligible((key[0], 'chr2:103:A:G', *key[2:]), items))
        self.assertIsNone(eligible((key[0], 'chr1:103:AA:G', *key[2:]), items))


class OptimizedExecution(unittest.TestCase):
    """The executable must not publish a success receipt for tampered inputs."""

    def test_source_and_report_gates_survive_python_optimization(self):
        here = Path(__file__).resolve().parent
        root = here.parent.parent
        source = root / 'research/epistasis-source-audit/siraj2026-S18.xlsx'
        report = root / 'research/epistasis-source-audit/s18-design-audit.json'
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            audit_dir = base / 'research/epistasis-source-audit'
            receipt_dir = base / 'research/epistasis-s18-independent'
            audit_dir.mkdir(parents=True)
            receipt_dir.mkdir(parents=True)
            receipt = receipt_dir / 'receipt.json'
            receipt.write_text('previous receipt\n')
            local_report = audit_dir / 's18-design-audit.json'
            local_report.write_bytes(report.read_bytes())

            # A ZIP with trailing bytes still parses, but no longer matches the
            # source commitment. A stripped hash assertion used to accept it.
            altered_source = base / 'altered.xlsx'
            altered_source.write_bytes(source.read_bytes() + b'tampered')
            result = subprocess.run([sys.executable, '-O', str(here / 'verify.py'),
                                     str(altered_source)], cwd=base, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn('unexpected S18 source hash', result.stderr)
            self.assertEqual(receipt.read_text(), 'previous receipt\n')

            altered_report = json.loads(report.read_text())
            altered_report['data_rows'] += 1
            local_report.write_text(json.dumps(altered_report))
            result = subprocess.run([sys.executable, '-O', str(here / 'verify.py'),
                                     str(source)], cwd=base, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn('data row count differs from report', result.stderr)
            self.assertEqual(receipt.read_text(), 'previous receipt\n')


if __name__ == '__main__':
    unittest.main()
