"""Failure-sensitive controls for metadata-only sealing before outcome exposure."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("seal_variant_transfer", ROOT / "scripts/seal_variant_transfer.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def component(ids, retained=True):
    return {"component": M.component_identity(ids), "variant_ids": sorted(ids), "size": len(ids),
            "retained_after_prior_exposure": retained, "exclusion_reasons": [] if retained else ["old_exposure"]}


def pair(identity, group):
    return {"id": identity, "component": group["component"], "locus": "source_parent_" + identity,
            "sequence_ref": "A" * 200, "sequence_alt": "A" * 99 + "T" + "A" * 100,
            "gc_preserving": True}


class SelectionTests(unittest.TestCase):
    def test_one_control_excludes_entire_component(self):
        a, b = component(["a", "b"]), component(["c", "d"])
        pairs = [pair(i, g) for g in [a, b] for i in g["variant_ids"]]
        cohort, assignments, ledger = M.select_cohort([a, b], pairs, {"a", "b", "c", "d"}, {"a"})
        self.assertEqual([r["id"] for r in cohort], ["c", "d"])
        self.assertIn("entire_component_contains_reused_S30_control", assignments["b"]["component_exclusion_reasons"])
        self.assertEqual(sum(r["final_cohort_member_count"] for r in ledger), 2)

    def test_missing_S4_membership_preserves_original_component_and_locus(self):
        group = component(["c", "d"])
        cohort, _, ledger = M.select_cohort([group], [pair("c", group), pair("d", group)], {"d"}, set())
        self.assertEqual(len(cohort), 1)
        self.assertEqual(cohort[0]["id"], "d")
        self.assertEqual(cohort[0]["component"], group["component"])
        self.assertNotEqual(cohort[0]["component"], M.component_identity(["d"]))
        self.assertEqual(cohort[0]["locus"], "source_parent_d")
        self.assertEqual(ledger[0]["original_size"], 2)

    def test_prior_excluded_component_is_never_restored(self):
        bad, good = component(["old"], False), component(["new"])
        cohort, _, _ = M.select_cohort([bad, good], [pair("new", good)], {"old", "new"}, set())
        self.assertEqual([r["id"] for r in cohort], ["new"])

    def test_incomplete_retained_DNA_is_rejected_before_membership_filter(self):
        group = component(["c", "d"])
        with self.assertRaisesRegex(ValueError, "complete retained components"):
            M.select_cohort([group], [pair("d", group)], {"d"}, set())

    def test_regrouped_label_is_rejected(self):
        group = component(["c", "d"])
        rows = [pair("c", group), pair("d", group)]
        rows[0]["component"] = M.component_identity(["c"])
        with self.assertRaisesRegex(ValueError, "changes original component"):
            M.select_cohort([group], rows, {"c", "d"}, set())

    def test_unknown_technical_control_is_rejected(self):
        group = component(["a"])
        with self.assertRaisesRegex(ValueError, "outside complete external population"):
            M.select_cohort([group], [pair("a", group)], {"a"}, {"missing_control"})

    def test_source_group_failure_and_integer_true_are_rejected(self):
        for value in [False, 1, "true", None]:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "failed or incomplete"):
                M.check_group_gate({"grouping_and_prior_exposure_gates_passed": value,
                                    "phase": "DNA_components_and_prior_exposure_complete_no_outcomes"})

    def test_unfinished_primary_group_phase_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "failed or incomplete"):
            M.check_group_gate({"grouping_and_prior_exposure_gates_passed": True, "phase": "prepared_not_executed"})

    def test_independent_comparison_must_have_actual_passes(self):
        comparison = {"matched": True, "failed_checks": [], "new_numeric_outcome_access": False,
                      "checks": [{"passed": True}], "checks_passed": 1,
                      "bindings": {"primary_manifest_sha256": "a", "independent_manifest_sha256": "b"}}
        M.check_comparison(comparison, "a", "b")
        comparison["checks"][0]["passed"] = False
        with self.assertRaisesRegex(ValueError, "comparison failed"):
            M.check_comparison(comparison, "a", "b")

    def test_independent_comparison_must_bind_both_manifests(self):
        comparison = {"matched": True, "failed_checks": [], "new_numeric_outcome_access": False,
                      "checks": [{"passed": True}], "checks_passed": 1,
                      "bindings": {"primary_manifest_sha256": "a", "independent_manifest_sha256": "b"}}
        with self.assertRaisesRegex(ValueError, "manifest bindings differ"):
            M.check_comparison(comparison, "changed", "b")


class InventoryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / "payload").write_bytes(b"frozen metadata\n")
        self.receipt = {"schema_version": 1, "files": {"payload": {
            "bytes": (self.root / "payload").stat().st_size, "sha256": M.sha256(self.root / "payload")}}}
        self.write_manifest()

    def tearDown(self):
        self.temp.cleanup()

    def write_manifest(self):
        (self.root / "manifest.json").write_text(json.dumps(self.receipt))
        self.digest = M.sha256(self.root / "manifest.json")

    def test_exact_inventory_and_tampered_payload(self):
        M.verify_inventory(self.root, "manifest.json", self.digest, exact=True)
        (self.root / "payload").write_bytes(b"changed metadata\n")
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            M.verify_inventory(self.root, "manifest.json", self.digest, exact=True)

    def test_wrong_manifest_pin(self):
        with self.assertRaisesRegex(ValueError, "manifest SHA256 mismatch"):
            M.verify_inventory(self.root, "manifest.json", "0" * 64)

    def test_unlisted_extra_file_rejected(self):
        (self.root / "unlisted").write_text("new")
        with self.assertRaisesRegex(ValueError, "extra or missing"):
            M.verify_inventory(self.root, "manifest.json", self.digest, exact=True)

    def test_wrong_byte_count_rejected(self):
        self.receipt["files"]["payload"]["bytes"] += 1
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, "byte count"):
            M.verify_inventory(self.root, "manifest.json", self.digest)

    def test_traversal_and_symlink_rejected(self):
        with self.assertRaisesRegex(ValueError, "unsafe manifest"):
            M.member(self.root, "../payload")
        (self.root / "link").symlink_to(self.root / "payload")
        with self.assertRaisesRegex(ValueError, "symlink"):
            M.member(self.root, "link")


if __name__ == "__main__":
    unittest.main()
