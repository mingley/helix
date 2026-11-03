import unittest

from audit import classify, features, rc, square_features


def embedded(first, second, first_alt, second_alt):
    """Construct two separated five-base neighborhoods and their SNVs."""
    dna = list("A" * 200)
    dna[48:53] = first
    dna[98:103] = second
    rr = "".join(dna)
    ar, ra, aa = list(rr), list(rr), list(rr)
    ar[50] = first_alt
    ra[100] = second_alt
    aa[50], aa[100] = first_alt, second_alt
    return {"RR": rr, "AR": "".join(ar), "RA": "".join(ra), "AA": "".join(aa)}


class FeatureChecks(unittest.TestCase):
    def test_integer_counts_and_global_reverse_complement(self):
        a = features("A" * 200)
        self.assertEqual(sum(a[:-1]), 198)
        self.assertEqual(a[0], 198)
        self.assertEqual(a[-1], 0)
        self.assertEqual(features("C" * 200)[-1], 200)
        sequence = ("ACGTTACGAGTC" * 17)[:200]
        self.assertEqual(features(sequence), features(rc(sequence)))
        with self.assertRaises(ValueError):
            features("A" * 199 + "N")

    def test_one_or_both_blind_edits_collapse_the_correct_constructs(self):
        for seqs, category, expected in [
            (embedded("AAATT", "AAAAA", "T", "C"), "first_blind_only",
             {("RR", "AR"), ("RA", "AA")}),
            (embedded("AAAAA", "AAATT", "C", "T"), "second_blind_only",
             {("RR", "RA"), ("AR", "AA")}),
            (embedded("AAATT", "AAATT", "T", "T"), "both_blind",
             {("RR", "AR"), ("RR", "RA"), ("RR", "AA"),
              ("AR", "RA"), ("AR", "AA"), ("RA", "AA")}),
        ]:
            first, second, actual, pairs = square_features(seqs)
            self.assertEqual(actual, category)
            self.assertEqual(classify(first, second), category)
            self.assertEqual({(p["construct1"], p["construct2"]) for p in pairs}, expected)
            self.assertTrue(all(not p["reverse_complement"] for p in pairs))

    def test_equal_and_opposite_nonzero_deltas_collapse_different_constructs(self):
        cases = [
            (embedded("AAAAC", "CAAAA", "C", "C"), "same_nonzero", ("AR", "RA")),
            (embedded("AAAAC", "AACAC", "C", "A"), "opposite_nonzero", ("RR", "AA")),
        ]
        for seqs, kind, expected in cases:
            d1, d2, actual, pairs = square_features(seqs)
            self.assertEqual(actual, kind)
            self.assertNotEqual(d1, (0,) * 33)
            self.assertNotEqual(d2, (0,) * 33)
            self.assertEqual([(p["construct1"], p["construct2"]) for p in pairs], [expected])

    def test_noncolliding_square_and_nearby_interaction(self):
        seqs = embedded("AAAAA", "AAAAA", "C", "G")
        self.assertEqual(square_features(seqs)[2:], ("no_collision", []))
        rr = "A" * 200
        ar, ra, aa = list(rr), list(rr), list(rr)
        ar[100] = "C"
        ra[102] = "C"
        aa[100], aa[102] = "C", "C"
        with self.assertRaisesRegex(ValueError, "parallelogram"):
            square_features({"RR": rr, "AR": "".join(ar), "RA": "".join(ra), "AA": "".join(aa)})


if __name__ == "__main__":
    unittest.main()
