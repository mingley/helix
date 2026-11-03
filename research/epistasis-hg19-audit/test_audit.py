import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location("hg19_audit", Path(__file__).with_name("audit.py"))
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def example():
    group = {"variant1": "chr1:105:A:G", "variant2": "chr1:110:C:T",
             "center": "var1", "window": "left", "distance": 5,
             "position1_zero_based": 4, "position2_zero_based": 9}
    seq = list("A" * 200)
    seq[9] = "C"
    row = (group["variant1"], group["variant2"], 5, "var1", "left", "RR", "".join(seq))
    return group, row


class CoordinateAndRowTests(unittest.TestCase):
    def test_start_requires_both_coordinates_and_forward_orientation(self):
        group, row = example()
        self.assertEqual(audit.interval_for_group(group), ("chr1", 100, 300, 4, 9))
        self.assertEqual(audit.check_rr_row(group, row)[-1], row[-1])
        group["position2_zero_based"] = 8
        with self.assertRaisesRegex(ValueError, "one valid forward insert"):
            audit.interval_for_group(group)

    def test_source_identity_and_alleles_are_checked(self):
        group, row = example()
        with self.assertRaisesRegex(ValueError, "row identity"):
            audit.check_rr_row(group, (*row[:-2], "AR", row[-1]))
        mutated = row[-1][:9] + "T" + row[-1][10:]
        with self.assertRaisesRegex(ValueError, "alleles"):
            audit.check_rr_row(group, (*row[:-1], mutated))

    def test_rejects_same_allele_or_distance_inconsistency(self):
        group, _ = example()
        group["variant1"] = "chr1:105:A:A"
        with self.assertRaisesRegex(ValueError, "not a substitution"):
            audit.interval_for_group(group)
        group, _ = example()
        group["distance"] = 4
        with self.assertRaisesRegex(ValueError, "distance"):
            audit.interval_for_group(group)


class ComparisonTests(unittest.TestCase):
    def test_exact_and_discrepant_sequences_are_preserved(self):
        _, row = example()
        source = row[-1]
        self.assertEqual(audit.differences(source, source, 100, 4, 9), [])
        reference = list(source)
        reference[4], reference[9], reference[50] = "G", "T", "C"
        self.assertEqual(audit.differences(source, "".join(reference), 100, 4, 9), [
            {"offset_zero_based": 4, "position_one_based": 105,
             "source_base": "A", "hg19_base": "G", "site": "variant1"},
            {"offset_zero_based": 9, "position_one_based": 110,
             "source_base": "C", "hg19_base": "T", "site": "variant2"},
            {"offset_zero_based": 50, "position_one_based": 151,
             "source_base": "A", "hg19_base": "C", "site": "background"},
        ])

    def test_length_mismatch_is_not_silently_truncated(self):
        with self.assertRaisesRegex(ValueError, "exactly 200"):
            audit.differences("A" * 200, "A" * 199, 100, 4, 9)


if __name__ == "__main__":
    unittest.main()
