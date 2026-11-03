import unittest

from test_verify import example_square
from verify import square
from verify_hg19 import expected_record


class Reference:
    def __init__(self, sequence):
        self.sequence = sequence
        self.requests = []

    def fetch(self, *args):
        self.requests.append(args)
        return self.sequence


class Hg19Controls(unittest.TestCase):
    def test_hg19_exact_match_and_coordinate_conversion(self):
        key, rows = example_square()
        item = square(key, rows)
        reference = Reference("A" * 200)
        result = expected_record(item, reference)
        self.assertEqual(reference.requests, [("chr1", 79, 279, "+")])
        self.assertEqual((result["hg19_start"], result["hg19_end"]), (79, 279))
        self.assertEqual(result["source_rr_row_1based"], 3)
        self.assertEqual(result["differences"], [])
        self.assertTrue(result["exact_match"])

    def test_differences_and_site_labels(self):
        key, rows = example_square()
        item = square(key, rows)
        genomic = list("A" * 200)
        for i in (3, 20, 23):
            genomic[i] = "T"
        result = expected_record(item, Reference("".join(genomic)))
        self.assertFalse(result["exact_match"])
        self.assertEqual(result["differences"], [
            {"offset_zero_based": 3, "position_one_based": 83, "source_base": "A", "hg19_base": "T", "site": "background"},
            {"offset_zero_based": 20, "position_one_based": 100, "source_base": "A", "hg19_base": "T", "site": "variant1"},
            {"offset_zero_based": 23, "position_one_based": 103, "source_base": "A", "hg19_base": "T", "site": "variant2"},
        ])

    def test_signed_start_disagreement_is_an_error(self):
        key, rows = example_square()
        item = square(key, rows)
        item["position2_zero_based"] += 1
        with self.assertRaises(ValueError):
            expected_record(item, Reference("A" * 200))


if __name__ == "__main__":
    unittest.main()
