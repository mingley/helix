use std::{
    fs,
    path::Path,
    process::{Command, Output},
};

fn invoke(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_helix"))
        .current_dir(Path::new(env!("CARGO_MANIFEST_DIR")).join("../.."))
        .args(args)
        .output()
        .unwrap()
}
fn success(output: Output) -> Output {
    assert!(
        output.status.success(),
        "stdout: {}\nstderr: {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}
#[test]
fn full_demo_checks_lean_is_repeatable_and_audits_tampering() {
    let temp = tempfile::tempdir().unwrap();
    let first = temp.path().join("first");
    let second = temp.path().join("second");
    success(invoke(&["demo", "--out", first.to_str().unwrap()]));
    success(invoke(&["demo", "--out", second.to_str().unwrap()]));
    let report: serde_json::Value =
        serde_json::from_slice(&fs::read(first.join("report.json")).unwrap()).unwrap();
    assert_eq!(report["proof"]["status"], "checked");
    assert_eq!(report["training_count"], 29);
    assert_eq!(report["heldout_count"], 25);
    assert_eq!(report["evaluation"]["count"], 25);
    assert_eq!(report["prevalence_baseline"]["count"], 25);
    assert!(report["evaluation"]["brier_score"]
        .as_f64()
        .unwrap()
        .is_finite());
    for name in [
        "report.json",
        "predictions.json",
        "evidence.json",
        "ledger.json",
        "proof.lean",
        "artifacts.json",
    ] {
        assert_eq!(
            fs::read(first.join(name)).unwrap(),
            fs::read(second.join(name)).unwrap(),
            "non-reproducible {name}"
        );
    }
    success(invoke(&["audit", first.to_str().unwrap()]));
    assert!(!invoke(&["demo", "--out", first.to_str().unwrap()])
        .status
        .success());
    fs::write(first.join("predictions.json"), "{}").unwrap();
    let result = invoke(&["audit", first.to_str().unwrap()]);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("checksum mismatch"));
}
#[test]
fn missing_lean_and_exhausted_budget_fail_and_skip_is_explicit() {
    let temp = tempfile::tempdir().unwrap();
    let out = temp.path().join("run");
    let path = out.to_str().unwrap();
    assert!(!invoke(&["demo", "--out", path, "--lean", "/no/such/lean"])
        .status
        .success());
    assert!(!out.exists());
    assert!(!invoke(&["demo", "--out", path, "--max-steps", "1"])
        .status
        .success());
    assert!(!out.exists());
    success(invoke(&["demo", "--out", path, "--skip-proof"]));
    let report: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("report.json")).unwrap()).unwrap();
    assert_eq!(report["proof"]["status"], "explicitly_skipped");
    assert!(report["proof"]["artifact"].is_null());
}

#[test]
fn research_reconstructs_reference_enumerates_mechanisms_and_preserves_commitments() {
    let temp = tempfile::tempdir().unwrap();
    let first = temp.path().join("research-first");
    let second = temp.path().join("research-second");
    success(invoke(&["research", "--out", first.to_str().unwrap()]));
    success(invoke(&["research", "--out", second.to_str().unwrap()]));
    let report: serde_json::Value =
        serde_json::from_slice(&fs::read(first.join("report.json")).unwrap()).unwrap();
    assert_eq!(report["reference_validation"]["transcript_bases"], 628);
    assert_eq!(
        report["reference_validation"]["coding_bases_including_stop"],
        444
    );
    assert_eq!(
        report["reference_validation"]
            ["protein_residues_excluding_stop_including_initiator_methionine"],
        147
    );
    assert_eq!(
        report["saturation"]["all_possible_single_base_substitutions"],
        4824
    );
    assert_eq!(report["saturation"]["coding_substitutions"], 1332);
    assert_eq!(
        report["saturation"]["canonical_splice_core_substitutions"],
        24
    );
    assert_eq!(
        report["saturation"]["selected_clinvar_variants_within_locus"],
        53
    );
    assert_eq!(report["annotation_comparison"]["agreements"], 54);
    assert_eq!(report["coding_annotation_comparison"]["agreements"], 17);
    assert_eq!(report["genome_context"]["gene_loci"], 29949);
    assert_eq!(report["promoter_assay"]["rows"]["all_source_rows"], 623);
    assert_eq!(
        report["promoter_assay"]["rows"]["reference_validated_snv_rows"],
        561
    );
    assert_eq!(report["promoter_assay"]["rows"]["source_equality_rows"], 62);
    assert_eq!(
        report["promoter_assay"]["mapping"]["target_reference_alignment_strand"],
        "+"
    );
    assert_eq!(
        report["promoter_assay"]["mapping"]["hbb_gene_transcription_strand"],
        "-"
    );
    assert_eq!(
        report["promoter_clinvar_overlap"].as_array().unwrap().len(),
        2
    );
    assert_eq!(
        report["statistics"]["coarse_category_information_limit"]
            ["best_possible_in_sample_fit_for_a_deterministic_category_only_classifier"],
        51
    );
    assert_eq!(
        report["scientific_status"]["new_biological_discovery"],
        false
    );
    let p = report["statistics"]["previous_classification_benchmark_post_hoc"]
        ["two_sided_exact_p_value"]
        .as_f64()
        .unwrap();
    assert!((p - 0.375).abs() < 1e-12);
    let artifacts: serde_json::Value =
        serde_json::from_slice(&fs::read(first.join("artifacts.json")).unwrap()).unwrap();
    for name in artifacts["files"].as_object().unwrap().keys() {
        assert_eq!(
            fs::read(first.join(name)).unwrap(),
            fs::read(second.join(name)).unwrap(),
            "non-reproducible research artifact: {name}"
        );
    }
    success(invoke(&["audit", first.to_str().unwrap()]));
    assert!(!invoke(&["research", "--out", first.to_str().unwrap()])
        .status
        .success());

    // A mixed-up report can have a correct file hash and still conflict with
    // the run's ledger commitments. Audit must reject that inconsistency.
    let mut changed = report;
    changed["input_commitments"]["protocol_sha256"] = serde_json::json!("a".repeat(64));
    let bytes = serde_json::to_vec(&changed).unwrap();
    fs::write(first.join("report.json"), &bytes).unwrap();
    let mut updated = artifacts;
    updated["files"]["report.json"] = serde_json::json!(genome_data::sha256(&bytes));
    fs::write(
        first.join("artifacts.json"),
        serde_json::to_vec(&updated).unwrap(),
    )
    .unwrap();
    let rejected = invoke(&["audit", first.to_str().unwrap()]);
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("input commitments differ"));
    let mut inconsistent: serde_json::Value =
        serde_json::from_slice(&fs::read(second.join("report.json")).unwrap()).unwrap();
    inconsistent["statistics"]["previous_classification_benchmark_post_hoc"]
        ["two_sided_exact_p_value"] = serde_json::json!(0.001);
    let bytes = serde_json::to_vec(&inconsistent).unwrap();
    fs::write(second.join("report.json"), &bytes).unwrap();
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(second.join("artifacts.json")).unwrap()).unwrap();
    manifest["files"]["report.json"] = serde_json::json!(genome_data::sha256(&bytes));
    fs::write(
        second.join("artifacts.json"),
        serde_json::to_vec(&manifest).unwrap(),
    )
    .unwrap();
    let rejected = invoke(&["audit", second.to_str().unwrap()]);
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("report statistics differs"));

    let invalid_hypotheses = temp.path().join("unrelated-hypotheses.json");
    fs::write(
        &invalid_hypotheses,
        serde_json::to_vec(&serde_json::json!({
            "schema_version":1,
            "hypotheses":[{"id":"unrelated-1"},{"id":"unrelated-2"},{"id":"unrelated-3"}]
        }))
        .unwrap(),
    )
    .unwrap();
    let rejected = invoke(&[
        "research",
        "--out",
        temp.path().join("invalid").to_str().unwrap(),
        "--hypotheses",
        invalid_hypotheses.to_str().unwrap(),
    ]);
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("five defined hypothesis IDs"));
}

#[test]
fn reference_catalog_fetch_and_summary_are_deterministic_and_hash_bound() {
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("tiny.2bit");
    // Independently packed twoBit records in non-alphabetical contig order.
    let mut bytes = [0x1a412743_u32, 0, 2, 0]
        .into_iter()
        .flat_map(u32::to_be_bytes)
        .collect::<Vec<_>>();
    bytes.extend_from_slice(b"\x04chrZ\x00\x00\x00\x22\x04chrA\x00\x00\x00\x44");
    bytes.extend(
        [7_u32, 1, 2, 2, 1, 1, 4, 0]
            .into_iter()
            .flat_map(u32::to_be_bytes),
    );
    bytes.extend([0x9c, 0x58]); // ACGTCCA with N at [2,4) and mask at [1,5).
    bytes.extend([4_u32, 0, 0, 0].into_iter().flat_map(u32::to_be_bytes));
    bytes.push(0x0e); // TTGA.
    fs::write(&path, &bytes).unwrap();
    let check = |command: &[&str]| {
        let mut args = vec![
            "reference",
            "--file",
            path.to_str().unwrap(),
            "--assembly",
            "GRCh38",
        ];
        args.extend_from_slice(command);
        invoke(&args)
    };
    let report = |command: &[&str]| -> serde_json::Value {
        serde_json::from_slice(&success(check(command)).stdout).unwrap()
    };

    let expected = genome_data::sha256(&bytes);
    let catalog = ["--expect-sha256", expected.as_str(), "catalog"];
    let indexed = report(&catalog);
    assert_eq!(indexed["schema_version"], 1);
    assert_eq!(indexed["operation"], "catalog");
    assert_eq!(indexed["reference"]["declared_assembly"], "GRCh38");
    assert_eq!(indexed["reference"]["byte_order"], "big");
    assert_eq!(indexed["reference"]["sha256"], expected.as_str());
    assert_eq!(indexed["contigs"][0]["name"], "chrA");
    assert_eq!(indexed["contigs"][1]["name"], "chrZ");
    assert_eq!(indexed["contigs"][1]["length"], 7);
    assert_eq!(indexed["contigs"][1]["n_blocks"], 1);
    assert_eq!(indexed["contigs"][1]["mask_blocks"], 1);
    assert_eq!(
        success(check(&catalog)).stdout,
        success(check(&catalog)).stdout
    );

    let extracted = report(&["fetch", "chrZ", "1", "6", "--strand=-", "--masked"]);
    assert_eq!(extracted["operation"], "fetch");
    assert_eq!(extracted["start0"], 1);
    assert_eq!(extracted["end0"], 6);
    assert_eq!(extracted["sequence"], "Ggnng");
    assert!(extracted["reference"]["sha256"].is_null());
    assert_eq!(report(&["fetch", "chrZ", "0", "7"])["sequence"], "ACNNCCA");

    let summary = report(&["summary"]);
    assert_eq!(summary["operation"], "summary");
    assert_eq!(summary["contigs"][0]["name"], "chrA");
    assert_eq!(
        summary["total"],
        serde_json::json!({"a": 3, "c": 3, "g": 1, "t": 2, "n": 2})
    );
    let split = report(&["summary", "--by-soft-mask"]);
    let counts = |a, c, g, t, n| serde_json::json!({"a": a, "c": c, "g": g, "t": t, "n": n});
    assert_eq!(split["operation"], "summary_by_soft_mask");
    assert_eq!(
        split["contigs"],
        serde_json::json!([
            {"name": "chrA", "length": 4, "counts": {
                "unmasked": counts(1, 0, 1, 2, 0), "masked": counts(0, 0, 0, 0, 0)
            }},
            {"name": "chrZ", "length": 7, "counts": {
                "unmasked": counts(2, 1, 0, 0, 0), "masked": counts(0, 2, 0, 0, 2)
            }}
        ])
    );
    assert_eq!(split["total"]["unmasked"], counts(3, 1, 1, 2, 0));
    assert_eq!(split["total"]["masked"], counts(0, 2, 0, 0, 2));

    let rejected = check(&["--expect-sha256", &"0".repeat(64), "catalog"]);
    assert!(!rejected.status.success());
    assert!(rejected.stdout.is_empty());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("SHA-256 mismatch"));
    assert!(!check(&["fetch", "Z", "0", "7"]).status.success());
    assert!(!check(&["fetch", "chrZ", "0", "8"]).status.success());

    let arguments = ["reference", "--file", path.to_str().unwrap()];
    for (extra, reason) in [
        (vec!["atlas-composition"], "declare --assembly GRCh38"),
        (
            vec!["atlas-composition", "--by-soft-mask"],
            "declare --assembly GRCh38",
        ),
        (
            vec!["--assembly", "GRCh37", "atlas-composition"],
            "atlas uses GRCh38 coordinates",
        ),
        (
            vec!["--assembly", "GRCh38", "atlas-composition"],
            "reference lacks primary atlas contig",
        ),
    ] {
        let rejected = invoke(&[arguments.as_slice(), extra.as_slice()].concat());
        assert!(!rejected.status.success());
        assert!(rejected.stdout.is_empty());
        assert!(String::from_utf8_lossy(&rejected.stderr).contains(reason));
    }
}
