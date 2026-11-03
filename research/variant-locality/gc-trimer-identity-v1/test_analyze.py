from fractions import Fraction
from itertools import product
import json
import shutil
import subprocess
import unittest

import analyze


class GCTrimerControls(unittest.TestCase):
    def test_all_short_sequences_and_reverse_complements(self):
        for length in (4, 5, 6):
            for bases in product("ACGT", repeat=length):
                sequence = "".join(bases)
                counts, g, t, boundary = analyze.features(sequence)
                self.assertEqual(3 * g, t + boundary)
                reverse = sequence.translate(analyze.COMPLEMENT)[::-1]
                self.assertEqual(analyze.features(reverse), (counts, g, t, boundary))

    def test_gc_is_redundant_for_interior_edits_but_not_at_edges(self):
        baseline = "A" * 200
        original = analyze.features(baseline)
        for index in (2, 100, 197):
            changed = baseline[:index] + "C" + baseline[index + 1:]
            _, g, t, boundary = analyze.features(changed)
            self.assertEqual(boundary, original[3])
            self.assertEqual(3 * (g - original[1]), t - original[2])
        for index in (0, 1, 198, 199):
            changed = baseline[:index] + "C" + baseline[index + 1:]
            _, g, t, boundary = analyze.features(changed)
            self.assertNotEqual(boundary, original[3])
            self.assertEqual(3 * (g - original[1]) - (t - original[2]), boundary - original[3])

    def test_same_trimers_can_have_distinct_gc_if_boundaries_change(self):
        first, second = "AAC" * 66 + "AA", "ACA" * 66 + "AC"
        c1, g1, t1, b1 = analyze.features(first)
        c2, g2, t2, b2 = analyze.features(second)
        self.assertEqual(c1, c2)
        self.assertEqual(t1, t2)
        self.assertEqual((g1, g2, b1, b2), (66, 67, 0, 3))

    def test_linear_scores_reparameterize_exactly(self):
        weights = {word: Fraction(i - 10, 17) for i, word in enumerate(analyze.CLASSES)}
        gc_weight = Fraction(-9, 11)
        for sequence in ("ACGT", "ACGTTAA", "AAC" * 66 + "AA", "ACA" * 66 + "AC"):
            left, right = analyze.exact_scores(sequence, weights, gc_weight)
            self.assertEqual(left, right)

    def test_rejects_invalid_sequences(self):
        for sequence in ("", "ACG", "ACGN", "acgt"):
            with self.assertRaises(ValueError):
                analyze.features(sequence)

    def test_published_fraction_is_exact_decimal_strings(self):
        summary = json.loads(analyze.PUBLISHED.read_text())
        stats = summary["correlation_sufficient_statistics"]
        n = summary["records"]
        covariance = n * stats["sum_GT"] - stats["sum_G"] * stats["sum_T"]
        var_g = n * stats["sum_G2"] - stats["sum_G"] ** 2
        var_t = n * stats["sum_T2"] - stats["sum_T"] ** 2
        expected = Fraction(covariance ** 2, var_g * var_t)
        fraction = summary["correlation_squared_exact"]
        self.assertEqual(fraction, {"numerator": str(expected.numerator),
                                    "denominator": str(expected.denominator)})
        self.assertGreater(expected.numerator, 2 ** 53)
        self.assertGreater(expected.denominator, 2 ** 53)

    @unittest.skipUnless(shutil.which("node"), "Node.js is not installed")
    def test_fraction_survives_javascript_json_and_bigint(self):
        receipt = analyze.PUBLISHED.parent.parent / "gc-trimer-independent-v1/receipt.json"
        script = r"""
const [summary, receipt] = process.argv.slice(1).map(JSON.parse);
const f = summary.correlation_squared_exact;
const checked = [f, receipt.correlation_squared_exact];
for (const item of checked) {
  if (typeof item.numerator !== 'string' || typeof item.denominator !== 'string' ||
      !/^[0-9]+$/.test(item.numerator) || !/^[0-9]+$/.test(item.denominator)) {
    throw new Error('exact fraction must use decimal strings');
  }
  if (BigInt(item.numerator) !== BigInt(f.numerator) ||
      BigInt(item.denominator) !== BigInt(f.denominator)) {
    throw new Error('receipt and summary fractions differ');
  }
}
const s = summary.correlation_sufficient_statistics;
const values = [summary.records, ...Object.values(s)];
if (!values.every(Number.isSafeInteger)) throw new Error('unsafe input statistic');
const n = BigInt(summary.records);
const g = BigInt(s.sum_G), t = BigInt(s.sum_T);
const cov = n * BigInt(s.sum_GT) - g * t;
const varG = n * BigInt(s.sum_G2) - g * g;
const varT = n * BigInt(s.sum_T2) - t * t;
if (BigInt(f.numerator) * varG * varT !== BigInt(f.denominator) * cov * cov) {
  throw new Error('fraction is not exact after JavaScript JSON parsing');
}
"""
        result = subprocess.run([shutil.which("node"), "-e", script,
                                 analyze.PUBLISHED.read_text(), receipt.read_text()],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
