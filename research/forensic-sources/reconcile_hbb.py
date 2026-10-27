#!/usr/bin/env python3
"""Reconcile the frozen MaveDB HBB rows against the author's original OSF table.

This is source-identity forensics, not predictive evaluation or a positional
score-pattern analysis. The raw dataset files are never modified.
"""
import argparse
import collections
import csv
import hashlib
import json
import pathlib
import re
import unittest
import xml.etree.ElementTree as ET
from decimal import Decimal

ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCES = pathlib.Path(__file__).resolve().parent
FROZEN = ROOT / "datasets/hbb-promoter-assay"
EXPECTED_ORIGINAL = "705694a6d45cd61159883250393b165a5813b40cc06e1dcdbe3270faad58af07"
EXPECTED_MAVEDB = "ec8b189f70c740f0ab2a11863a0a464fb421d93fefb9da1240cd3831872f57c6"
GENOMIC_START = 5227021
PATTERN = re.compile(r"n\.(\d+)(?:([ACGT])>([ACGT])|=)")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_rows(path, delimiter=","):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def reconcile(mave_rows, original_rows, target):
    require(len(target) == 187, "unexpected target length")
    original = {}
    for line_number, row in enumerate(original_rows, 2):
        position, alternate = int(row["Position"]), row["Alt"]
        require(row["Chromosome"] == "11", "unexpected original chromosome")
        n = position - GENOMIC_START
        require(1 <= n <= len(target), "original position outside target")
        require(row["Ref"] == target[n - 1], "original reference allele mismatch")
        require(alternate == "-" or (alternate in "ACGT" and alternate != row["Ref"]),
                "unsupported original allele")
        key = position, alternate
        require(key not in original, "duplicate original position/alternate key")
        original[key] = line_number, row
    seen = set()
    accessions = set()
    deletion_mappings = []
    max_pvalue_delta = Decimal(0)
    pvalue_zero_mavedb = 0
    pvalue_zero_original = 0
    snv_count = 0
    equivalents = collections.defaultdict(list)
    for row in mave_rows:
        require(row["accession"] not in accessions, "duplicate MaveDB accession")
        accessions.add(row["accession"])
        match = PATTERN.fullmatch(row["hgvs_nt"])
        require(match is not None, "unsupported MaveDB HGVS representation")
        n = int(match[1])
        require(1 <= n <= len(target), "MaveDB position outside target")
        alternate = match[3] or "-"
        if match[2]:
            require(match[2] == target[n - 1], "MaveDB reference allele mismatch")
        key = GENOMIC_START + n, alternate
        require(key in original, "MaveDB record has no original position/allele")
        require(key not in seen, "MaveDB records do not match original bijectively")
        seen.add(key)
        line_number, source = original[key]
        for mave, original_column in (
            ("score", "Value"), ("unique_tags", "Tags"),
            ("dna_sequences", "DNA"), ("rna_sequences", "RNA")
        ):
            require(Decimal(row[mave]) == Decimal(source[original_column]),
                    f"numeric mismatch: {row['accession']} {mave}")
        original_pvalue = Decimal(source["P-Value"])
        mave_pvalue = Decimal(row["p-value"])
        delta = abs(mave_pvalue - original_pvalue)
        require(delta <= Decimal("0.000005"), "p-value differs beyond five-decimal rounding")
        max_pvalue_delta = max(max_pvalue_delta, delta)
        pvalue_zero_mavedb += mave_pvalue == 0
        pvalue_zero_original += original_pvalue == 0
        if alternate == "-":
            altered_target = target[:n - 1] + target[n:]
            altered_hash = digest(altered_target.encode())
            equivalents[altered_hash].append(n)
            deletion_mappings.append({
                "mavedb_accession": row["accession"],
                "mavedb_literal_hgvs": row["hgvs_nt"],
                "original_table_line": line_number,
                "original_chromosome": source["Chromosome"],
                "original_position_one_based": int(source["Position"]),
                "original_reference": source["Ref"],
                "original_alternate": "-",
                "source_anchored_deletion_description": f"n.{n}del",
                "normalization_status": "source-anchored, not full-reference HGVS 3-prime normalized",
                "altered_target_sha256": altered_hash,
                "score": source["Value"],
                "pvalue_original": source["P-Value"],
                "pvalue_mavedb": row["p-value"],
                "tags": int(source["Tags"]),
                "dna_count": int(source["DNA"]),
                "rna_count": int(source["RNA"]),
                "matched_fields": [
                    "genomic_position", "original_allele", "score",
                    "unique_tags", "dna_sequences", "rna_sequences", "rounded_pvalue"
                ],
            })
        else:
            snv_count += 1
    require(len(seen) == len(original) == len(mave_rows), "unmatched original/MaveDB rows")
    duplicate_groups = sorted(
        [sorted(positions) for positions in equivalents.values() if len(positions) > 1]
    )
    by_altered_sequence = collections.defaultdict(list)
    for mapping in deletion_mappings:
        by_altered_sequence[mapping["altered_target_sha256"]].append(mapping)
    coefficient_groups = []
    total_sse = Decimal(0)
    for altered_hash, mappings in by_altered_sequence.items():
        if len(mappings) < 2:
            continue
        mappings.sort(key=lambda row: row["original_position_one_based"])
        scores = [Decimal(row["score"]) for row in mappings]
        mean = sum(scores) / len(scores)
        sse = sum((score - mean) ** 2 for score in scores)
        total_sse += sse
        coefficient_groups.append({
            "altered_target_sha256": altered_hash,
            "target_positions_one_based": [
                row["original_position_one_based"] - GENOMIC_START for row in mappings
            ],
            "mavedb_accessions": [row["mavedb_accession"] for row in mappings],
            "original_coefficients": [row["score"] for row in mappings],
            "original_tags": [row["tags"] for row in mappings],
            "group_mean_coefficient": str(mean),
            "sum_squared_deviations_from_group_mean": str(sse),
            "coefficient_range": str(max(scores) - min(scores)),
            "distinct_coefficient_count": len(set(scores)),
        })
    coefficient_groups.sort(key=lambda group: group["target_positions_one_based"])
    duplicate_rows = sum(len(group["target_positions_one_based"]) for group in coefficient_groups)
    qualified_by_identity = collections.defaultdict(list)
    for row in deletion_mappings:
        if row["tags"] >= 10:
            qualified_by_identity[row["altered_target_sha256"]].append(row)
    qualified_sse = Decimal(0)
    for group in qualified_by_identity.values():
        scores = [Decimal(row["score"]) for row in group]
        mean = sum(scores) / len(scores)
        qualified_sse += sum((score - mean) ** 2 for score in scores)
    qualified_rows = sum(map(len, qualified_by_identity.values()))
    source_estimate_inconsistency = {
        "scope": "Original deletion source rows, not SNV evaluation, and not a biological noise floor.",
        "identity": "Exact equality of the complete altered187-nt reporter target sequence.",
        "groups": coefficient_groups,
        "groups_with_distinct_coefficients": sum(
            group["distinct_coefficient_count"] > 1 for group in coefficient_groups
        ),
        "total_within_identity_sse": str(total_sse),
        "lower_bound_mse_over_all_deletion_rows": str(total_sse / len(deletion_mappings)) if deletion_mappings else None,
        "lower_bound_mse_over_repeated_identity_rows": str(total_sse / duplicate_rows) if duplicate_rows else None,
        "mean_group_sum_squared_deviation": str(total_sse / len(coefficient_groups)) if coefficient_groups else None,
        "denominators": {
            "all_source_deletion_rows": len(deletion_mappings),
            "rows_in_repeated_identity_groups": duplicate_rows,
            "repeated_identity_groups": len(coefficient_groups),
        },
        "derivation": "For each identical input x with observed source coefficients y1..yn, sum_i(yi-f(x))^2 = sum_i(yi-ybar)^2 + n*(f(x)-ybar)^2. Summing nonnegative residual terms gives the finite observed-label lower bound.",
        "interpretation": "Different additive-regression coefficients attached to identical isolated edited target sequences are source-estimate inconsistency. Their difference can reflect representation, multi-mutant background assignment, estimation or correlated measurement; it does not demonstrate that identical DNA has different intrinsic biology.",
        "uncertainty": "No replicate coefficient covariance or biological error distribution is available here; these algebraic bounds are neither confidence intervals nor estimates of latent experimental noise.",
        "predefined_minimum_ten_tag_qc": {
            "retained_deletion_rows": qualified_rows,
            "distinct_altered_target_sequences": len(qualified_by_identity),
            "repeated_identity_groups": sum(len(group) > 1 for group in qualified_by_identity.values()),
            "total_within_identity_sse": str(qualified_sse),
            "lower_bound_mse_from_repeated_identity": str(qualified_sse / qualified_rows) if qualified_rows else None,
            "raw_groups_with_exactly_one_low_support_row": sum(
                sum(tags < 10 for tags in group["original_tags"]) == 1
                for group in coefficient_groups
            ),
            "interpretation": "Apply the study's preexisting minimum10-tag QC before interpreting source-label inconsistency. In this HBB snapshot every raw equivalent pair contains one low-support row, so this duplicate-identity bound vanishes after QC; that does not establish noise-free measurements.",
        },
    }
    return {
        "source_rows": len(original_rows),
        "database_rows": len(mave_rows),
        "one_to_one_matches": len(seen),
        "snv_matches": snv_count,
        "equality_to_original_deletion_matches": len(deletion_mappings),
        "exact_numeric_fields": ["score", "unique_tags", "dna_sequences", "rna_sequences"],
        "pvalue_absolute_tolerance": "0.000005",
        "maximum_pvalue_absolute_difference": str(max_pvalue_delta),
        "mavedb_literal_zero_pvalues": pvalue_zero_mavedb,
        "original_literal_zero_pvalues": pvalue_zero_original,
        "deletion_rows_below_ten_tags": sum(int(row["Tags"]) < 10 and row["Alt"] == "-"
                                           for row in original_rows),
        "snv_rows_below_ten_tags": sum(int(row["Tags"]) < 10 and row["Alt"] != "-"
                                      for row in original_rows),
        "deletion_rows_at_least_ten_tags": sum(int(row["Tags"]) >= 10 and row["Alt"] == "-"
                                              for row in original_rows),
        "distinct_deletion_altered_target_sequences": len(equivalents),
        "equivalent_deletion_position_groups": duplicate_groups,
        "rows_in_equivalent_deletion_groups": sum(map(len, duplicate_groups)),
        "source_estimate_inconsistency": source_estimate_inconsistency,
        "deletion_mapping": deletion_mappings,
    }


def clinical_ancestry():
    root = ET.parse(SOURCES / "clinvar-overlap-variationid.xml")
    records = []
    for archive in root.iter("VariationArchive"):
        pmids = sorted({element.text for element in archive.iter("ID")
                        if element.attrib.get("Source") == "PubMed"})
        records.append({
            "accession": archive.attrib["Accession"] + "." + archive.attrib["Version"],
            "variation_name": archive.attrib["VariationName"],
            "pubmed_ids_explicitly_cited": pmids,
            "kircher_2019_pmid_explicitly_cited": "31395865" in pmids,
        })
    require(len(records) == 2, "unexpected overlapping clinical record count")
    return {
        "scope": "the two selected ClinVar SNVs overlapping the reporter fragment",
        "retrieval_state": "live XML record versions match the frozen ESummary versions",
        "records": records,
        "conclusion": "No explicit Kircher2019 citation was found in these XML records. Absence of a citation does not establish independent data, cohorts, organizations, or uncited evidence ancestry.",
    }


def compute_report():
    retained_sources = json.loads((SOURCES / "retrieval.json").read_text())
    for item in retained_sources:
        source_path = SOURCES / item["file"]
        require(source_path.parent == SOURCES, "retained source path escapes source directory")
        source_bytes = source_path.read_bytes()
        require(len(source_bytes) == item["size_bytes"], "retained source byte count changed")
        require(digest(source_bytes) == item["sha256"], "retained source commitment changed")
    original_path = SOURCES / "GRCh38_HBB.tsv"
    mave_path = FROZEN / "scores.csv"
    require(digest(original_path.read_bytes()) == EXPECTED_ORIGINAL, "original OSF table hash changed")
    require(digest(mave_path.read_bytes()) == EXPECTED_MAVEDB, "frozen MaveDB scores hash changed")
    listing = json.loads((SOURCES / "osf-grch38-files.json").read_text())
    metadata = next(item["attributes"] for item in listing["data"]
                    if item["attributes"]["name"] == "GRCh38_HBB.tsv")
    require(metadata["extra"]["hashes"]["sha256"] == EXPECTED_ORIGINAL,
            "download hash disagrees with OSF service")
    require(metadata["current_version"] == 1, "original OSF file version changed")
    dataset_metadata = json.loads((FROZEN / "metadata.json").read_text())
    target = dataset_metadata["targetGenes"][0]["targetSequence"]["sequence"]
    reconciliation = reconcile(read_rows(mave_path), read_rows(original_path, "\t"), target)
    require(reconciliation["one_to_one_matches"] == 623, "unexpected match count")
    require(reconciliation["equality_to_original_deletion_matches"] == 62, "unexpected deletion count")
    literature = json.loads((SOURCES / "novelty-literature.json").read_text())
    source_tree = json.loads((SOURCES / "github-tree.json").read_text())
    source_commit = json.loads((SOURCES / "github-commit.json").read_text())
    require(source_commit["sha"] == source_tree["sha"], "author source ref changed between retained metadata requests")
    protocol = ROOT / "experiments/regulatory-prediction-v1.json"
    return {
        "schema_version": 1,
        "analysis_type": "source-identity reconciliation and adversarial measurement audit",
        "score_exposure_disclosure": "All623 source/database rows were read for identity matching, including coefficients and counts. No positional score-pattern analysis, predictor fitting or performance-based model selection was conducted in this script.",
        "conclusion": {
            "status": "confirmed_original_deletions_represented_as_database_equalities",
            "supported": "All62 literal MaveDB n.<position>= records match original source single-base deletion rows by position, coefficient and all three count columns; p-values match five-decimal rounding.",
            "not_established": [
                "The exact import, normalization or migration code path that introduced the equality notation",
                "The date when the representation discrepancy first appeared",
                "That62 source rows represent62 distinct physical alleles after repeat normalization",
                "Global novelty of a legacy MaveDB deletion/equality discrepancy",
                "Independent biological replication, endogenous expression effects or clinical pathogenicity",
            ],
        },
        "input_commitments": {
            "mavedb_scores_sha256": EXPECTED_MAVEDB,
            "original_osf_table_sha256": EXPECTED_ORIGINAL,
            "original_osf_service_hash_matches": True,
            "original_osf_file_version": metadata["current_version"],
            "original_osf_file_created_at": metadata["date_created"],
            "original_osf_materialized_path": metadata["materialized_path"],
            "target_sequence_sha256": digest(target.encode()),
            "source_documentation_git_commit": source_commit["sha"],
            "source_documentation_commit_tree": source_commit["commit"]["tree"]["sha"],
        },
        "reconciliation": reconciliation,
        "measurement_cautions": [
            "The source library used error-prone-PCR multi-mutant constructs. Coefficients are fitted additive effects, not every allele isolated in one clone.",
            "Three transfections share a constructed library; barcode counts are not independent biological replications.",
            "The primary paper explicitly warns that coefficient p-values are inflated by tag-independence assumptions, correlated counts and unmodeled epistasis.",
            "The original raw table retains27 deletion coefficients with fewer than10 tags. The author's visualization filters minimum tag support at display/analysis time; retention is not proof the published downstream filter was violated.",
            "All561 substitution rows have at least10 tags, so the declared substitution-only evaluation is unaffected by the low-tag deletion rows.",
            "Database literal p=0 values are rounded, not true mathematical zeros or infinite evidence.",
            "The62 source deletions yield53 distinct altered187-nt target sequences. Equivalent rows require grouping and provenance-aware resolution before future deletion modeling.",
            "Counts and sequence opportunity may be predictors of estimation noise; held-out sequence positions still share library-wide regression and measurement construction.",
        ],
        "clinical_ancestry": clinical_ancestry(),
        "evaluation_protocol_audit": {
            "protocol": "experiments/regulatory-prediction-v1.json",
            "protocol_sha256": digest(protocol.read_bytes()),
            "role": "Adversarial interpretation review after protocol freeze, without changing the fixed modeling plan.",
            "checks": [
                {
                    "topic": "feature-support leakage",
                    "assessment": "Radius10 sequence supports are disjoint at center distance>20. Grouping alleles by position and guarded contiguous folds addresses direct feature overlap, not source-estimator covariance or biological dependence.",
                },
                {
                    "topic": "jointly estimated assay labels",
                    "assessment": "HBB coefficients were jointly fitted from multi-mutant constructs before downstream cross-validation. Disjoint sequence windows do not create independent coefficient-estimation data; unknown covariance remains.",
                },
                {
                    "topic": "locus and context transport",
                    "assessment": "Nine leave-element-out targets are stronger than random variants but still one selected study. Cell-line and construct differences are real transport conditions, not interchangeable independent human promoters.",
                },
                {
                    "topic": "measurement support and source formatting",
                    "assessment": "Minimum10-tag QC is prespecified, not selected to improve performance. All561HBB substitutions meet it. Do not derive significance labels from rounded MaveDB p=0 or interpret low-support raw deletion conflicts as biological noise floors.",
                },
                {
                    "topic": "context-blind bounds",
                    "assessment": "Finite equal-input observed-label bounds are algebraic statements about the measured source coefficients. They mix technical variation with possible biology and do not identify either cause or latent error distribution.",
                },
                {
                    "topic": "resampling and rotations",
                    "assessment": "Nine-element and position-block resampling conditions on overlapping fitted OOF predictions. Chemistry-stratified rotations break some within-position and library dependence. Retain their declared descriptive/falsification status, not calibrated biological p-values or population confidence intervals.",
                },
                {
                    "topic": "outcome exposure",
                    "assessment": "Forensic matching necessarily read the already exposed HBB coefficients. No new-promoter outcome table was inspected by this analysis, and no model, threshold, motif or split was selected from forensic matches.",
                },
            ],
        },
        "novelty_check": {
            "status": "related_legacy_discrepancy_already_publicly_acknowledged",
            "prior_public_issue": {
                "url": "https://github.com/Open-Athena/VEP-bench/issues/32",
                "created_at": "2026-09-02T21:24:16Z",
                "finding": "The issue explicitly warns that a legacy MaveDB F9 source deletion appears as a target-equality expression and recommends source-crosswalk validation.",
            },
            "additional_related_issue": "https://github.com/Open-Athena/VEP-bench/issues/62",
            "bounded_literature_query": literature["request"]["queryString"],
            "bounded_literature_hit_count": literature["hitCount"],
            "interpretation": "This is a reproduced source-quality finding and HBB-specific reconciliation. A bounded search is not proof of novelty; no globally novel discrepancy is claimed.",
        },
        "reuse": {
            "mavedb_score_set": "CC0 1.0 as retained in the frozen score-set metadata",
            "primary_article": "CC BY4.0; Kircher et al.2019, DOI10.1038/s41467-019-11526-w, PMID31395865",
            "original_osf_table": "Public author deposit with the producer's historical data-use statement. No separate OSF Creative Commons grant was established; do not relabel it CC0.",
            "retained_author_source_code_and_documentation": "GPL-3.0 repository snapshot; license retained in author-GPL-3.0-LICENSE. It is reference evidence, not incorporated into Helix's MIT implementation.",
        },
    }


class ReconciliationTests(unittest.TestCase):
    def setUp(self):
        self.target = "G" + "A" * 186
        self.original = {
            "Chromosome": "11", "Position": "5227022", "Ref": "G", "Alt": "-",
            "Tags": "68", "DNA": "231", "RNA": "538", "Value": "-0.51",
            "P-Value": "0.00289870880628956",
        }
        self.mave = {
            "accession": "example#1", "hgvs_nt": "n.1=", "score": "-0.51",
            "unique_tags": "68.0", "dna_sequences": "231.0", "rna_sequences": "538.0",
            "p-value": "0.0029",
        }

    def test_all_independent_numeric_fields_are_required(self):
        for field in ["score", "unique_tags", "dna_sequences", "rna_sequences", "p-value"]:
            changed = dict(self.mave, **{field: "999"})
            with self.assertRaises(ValueError):
                reconcile([changed], [self.original], self.target)

    def test_identity_and_bijection_are_required(self):
        for hgvs in ["n.0=", "n.188=", "n.2=", "n.1A>T"]:
            with self.assertRaises(ValueError):
                reconcile([dict(self.mave, hgvs_nt=hgvs)], [self.original], self.target)
        with self.assertRaises(ValueError):
            reconcile([self.mave, self.mave], [self.original], self.target)

    def test_original_deletion_is_not_reclassified_as_wild_type(self):
        report = reconcile([self.mave], [self.original], self.target)
        self.assertEqual(report["equality_to_original_deletion_matches"], 1)
        self.assertEqual(report["deletion_mapping"][0]["source_anchored_deletion_description"], "n.1del")
        self.assertEqual(report["deletion_mapping"][0]["mavedb_literal_hgvs"], "n.1=")
        self.assertNotEqual(report["deletion_mapping"][0]["altered_target_sha256"], digest(self.target.encode()))

    def test_preexisting_qc_can_remove_a_raw_identity_inconsistency(self):
        target = "A" * 187
        originals = [
            dict(self.original, Ref="A", Position="5227022", Tags="1", Value="-1"),
            dict(self.original, Ref="A", Position="5227023", Tags="12", Value="1"),
        ]
        maves = [
            dict(self.mave, accession="example#1", hgvs_nt="n.1=", unique_tags="1", score="-1"),
            dict(self.mave, accession="example#2", hgvs_nt="n.2=", unique_tags="12", score="1"),
        ]
        result = reconcile(maves, originals, target)["source_estimate_inconsistency"]
        self.assertEqual(Decimal(result["lower_bound_mse_over_all_deletion_rows"]), 1)
        qc = result["predefined_minimum_ten_tag_qc"]
        self.assertEqual(qc["retained_deletion_rows"], 1)
        self.assertEqual(qc["repeated_identity_groups"], 0)
        self.assertEqual(Decimal(qc["lower_bound_mse_from_repeated_identity"]), 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Write the derived report, never source files")
    parser.add_argument("--verify", action="store_true", help="Verify the existing report against frozen inputs")
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    if arguments.self_test:
        result = unittest.main(argv=[__file__], exit=False).result
        raise SystemExit(0 if result.wasSuccessful() else 1)
    else:
        report = compute_report()
        destination = ROOT / "research/regulatory-forensics.json"
        if arguments.write:
            destination.write_text(json.dumps(report, indent=2) + "\n")
        if arguments.verify:
            require(json.loads(destination.read_text()) == report, "saved forensic report is stale")
        print(json.dumps({
            "status": report["conclusion"]["status"],
            "matched": report["reconciliation"]["one_to_one_matches"],
            "deletion_mappings": report["reconciliation"]["equality_to_original_deletion_matches"],
            "global_novelty_claimed": False,
        }))
