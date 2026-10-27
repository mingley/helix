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
