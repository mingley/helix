import io
import unittest

from audit_design_crosswalk import (
    DESIGN_FIELDS, evaluate_alias, fasta_aliases, fasta_records, select_aliases, tsv_records,
)


class DesignReferenceTests(unittest.TestCase):
    def test_literal_and_compound_aliases_keep_conflicting_sequences(self):
        stream = io.StringIO(">a descriptive text\nAC\nGT\n>(a;b)\nTTGA\n>c\nTCGA\n")
        index = select_aliases(fasta_records(stream), {"a", "b", "absent"}, fasta_aliases)
        self.assertEqual(evaluate_alias(index["a"], "ACGT")["status"], "conflicting_sequences")
        self.assertEqual(evaluate_alias(index["b"], "TTGA")["status"], "exact_sequence")
        self.assertEqual(evaluate_alias(index["b"], "TCAA")["status"], "sequence_mismatch")
        self.assertEqual(evaluate_alias(index.get("absent", {}), "ACGT")["status"], "missing")

    def test_duplicate_aliases_with_same_dna_keep_both_provenance_types(self):
        stream = io.StringIO(">a\nacgt\n>(b;a)\nACGT\n")
        index = select_aliases(fasta_records(stream), {"a"}, fasta_aliases)
        result = evaluate_alias(index["a"], "ACGT")
        self.assertEqual(result["status"], "exact_sequence")
        self.assertEqual(result["sequence_hashes"][0]["records_by_alias_type"],
                         {"literal_fasta_header": 1, "compound_fasta_member": 1})

    def test_sequence_design_header_is_required_for_tsv(self):
        with self.assertRaisesRegex(ValueError, "design TSV header differs"):
            list(tsv_records(io.StringIO("id\tRNA\nexample\t12\n")))
        line = "\t".join(["alias", "acgt"] + [""] * (len(DESIGN_FIELDS) - 2))
        stream = io.StringIO("\t".join(DESIGN_FIELDS) + "\n" + line + "\n")
        self.assertEqual(list(tsv_records(stream)), [("alias", "ACGT")])


if __name__ == "__main__":
    unittest.main()
