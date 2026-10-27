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
