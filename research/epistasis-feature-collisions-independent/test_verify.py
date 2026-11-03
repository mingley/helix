import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from verify import (PAIRS, classify_equalities, expected_features, feature,
                    reverse_complement, square)


def change(seq, index, base):
    return seq[:index] + base + seq[index + 1:]


def make_square(ref, p, first, q, second):
    ar = change(ref, p, first)
    ra = change(ref, q, second)
    return {"RR": ref, "AR": ar, "RA": ra, "AA": change(ar, q, second)}


class FullSequenceControls(unittest.TestCase):
    def test_reverse_complement_and_fixed_length(self):
        dna = ("ACACTTAGCGGCTA" * 15)[:200]
        self.assertEqual(len(dna), 200)
        self.assertEqual(feature(dna), feature(reverse_complement(dna)))
        self.assertEqual(sum(feature(dna)[:32]), 198)
        self.assertEqual(feature(dna)[-1], dna.count("C") + dna.count("G"))
        with self.assertRaises(ValueError):
            feature(dna[:-1])
        with self.assertRaises(ValueError):
            feature("N" + dna[1:])

    def test_all_collision_categories_from_complete_sequences(self):
        plain = "A" * 200
        blind_at_first = plain[:28] + "AAATT" + plain[33:]
        blind_at_both = blind_at_first[:78] + "AAATT" + blind_at_first[83:]
        opposite_reference = change(plain, 80, "C")
        cases = [
            (make_square(blind_at_both, 30, "T", 80, "T"), "both_blind", 6),
            (make_square(blind_at_first, 30, "T", 80, "C"), "first_blind_only", 2),
            (make_square(blind_at_first, 80, "C", 30, "T"), "second_blind_only", 2),
            (make_square(plain, 30, "C", 80, "C"), "same_nonzero", 1),
            (make_square(opposite_reference, 30, "C", 80, "A"), "opposite_nonzero", 1),
            (make_square(plain, 30, "C", 80, "T"), "no_collision", 0),
        ]
        for seqs, category, pairs in cases:
            with self.subTest(category=category):
                result = expected_features(seqs)
                self.assertEqual(result["classification"], category)
                self.assertEqual(len(result["feature_equal_pairs"]), pairs)
                self.assertTrue(all(not p["reverse_complement"] for p in result["feature_equal_pairs"]))

    def test_feature_additivity_and_valid_equality_pattern_are_enforced(self):
        with self.assertRaisesRegex(ValueError, "additivity"):
            expected_features(make_square("A" * 200, 30, "C", 31, "G"))
        with self.assertRaisesRegex(ValueError, "unsupported feature equality"):
            classify_equalities([PAIRS[0]])

    def test_source_square_checks_all_four_labels_and_edits(self):
        seqs = make_square("A" * 200, 30, "C", 80, "G")
        key = ("chr1:1000:A:C", "chr1:1050:A:G", "var1", "left")
        rows = [(n, 50, label, seqs[label]) for n, label in enumerate(("RR", "AR", "RA", "AA"), 3)]
        prior, recovered = square(key, rows)
        self.assertEqual((prior["position1_zero_based"], prior["position2_zero_based"]), (30, 80))
        self.assertEqual(recovered, seqs)
        damaged = list(rows)
        n, d, label, seq = damaged[-1]
        damaged[-1] = (n, d, label, change(seq, 150, "T"))
        self.assertIsNone(square(key, damaged))
        self.assertIsNone(square(key, rows[:-1]))


class OptimizedIntegrityControls(unittest.TestCase):
    def test_optimized_rejects_altered_source_and_ledger_without_writing_receipt(self):
        here = Path(__file__).resolve().parent
        primary = here.parent / "epistasis-feature-collisions-v1/results/ledger.jsonl"
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            receipt = tmp / "receipt.json"
            receipt.write_text("prior receipt\n")
            bad_source = tmp / "source.xlsx"
            bad_source.write_bytes(b"not the pinned S18 source")
            run = subprocess.run([sys.executable, "-O", str(here / "verify.py"),
                                  "--source", str(bad_source), "--receipt", str(receipt)],
                                 capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("source input SHA-256 mismatch", run.stderr)
            self.assertEqual(receipt.read_text(), "prior receipt\n")

            bad_ledger = tmp / "ledger.jsonl"
            lines = primary.read_text().splitlines()
            first = json.loads(lines[0])
            first["first_feature_delta"][0] += 1
            lines[0] = json.dumps(first, sort_keys=True, separators=(",", ":"))
            bad_ledger.write_text("\n".join(lines) + "\n")
            run = subprocess.run([sys.executable, "-O", str(here / "verify.py"),
                                  "--ledger", str(bad_ledger), "--receipt", str(receipt)],
                                 capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("feature ledger line 1 differs", run.stderr)
            self.assertEqual(receipt.read_text(), "prior receipt\n")


if __name__ == "__main__":
    unittest.main()
