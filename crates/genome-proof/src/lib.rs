//! A subprocess boundary to the pinned Lean authority. Missing/failed Lean is an error.
use anyhow::{bail, ensure, Context, Result};
use genome_logic::Derivation;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::Path,
    process::{Command, Stdio},
    thread,
    time::{Duration, Instant},
};

pub const LEAN_VERSION: &str = "4.19.0";
const PRELUDE: &str = include_str!("../../../proofs/Helix.lean");

#[derive(Debug, Serialize)]
pub struct CheckedProof {
    pub authority: String,
    pub source_sha256: String,
    pub scope: String,
}
pub fn source(derivation: &Derivation, arithmetic: &str) -> String {
    format!(
        "{PRELUDE}\nnamespace Helix.Generated\n{}\n{arithmetic}\nend Helix.Generated\n",
        derivation.lean()
    )
}
pub fn check(lean: &Path, source: &str) -> Result<CheckedProof> {
    validate_source(source)?;
    let version = Command::new(lean)
        .arg("--version")
        .output()
        .context("Lean unavailable; run scripts/setup.sh")?;
    ensure!(version.status.success(), "Lean version command failed");
    let version = String::from_utf8(version.stdout)?.trim().to_owned();
    ensure!(
        version.starts_with(&format!("Lean (version {LEAN_VERSION},")),
        "expected Lean {LEAN_VERSION}, found {version}"
    );
    let temp = tempfile::tempdir()?;
    let path = temp.path().join("Obligations.lean");
    fs::write(&path, source)?;
    let out = fs::File::create(temp.path().join("stdout"))?;
    let err = fs::File::create(temp.path().join("stderr"))?;
    let mut child = Command::new(lean)
        .arg(&path)
        .stdin(Stdio::null())
        .stdout(out)
        .stderr(err)
        .spawn()?;
    let start = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait()? {
            break status;
        }
        if start.elapsed() > Duration::from_secs(30) {
            child.kill()?;
            child.wait()?;
            bail!("Lean proof check exceeded 30 second budget");
        }
        thread::sleep(Duration::from_millis(20));
    };
    ensure!(
        status.success(),
        "Lean rejected proof:\n{}\n{}",
        fs::read_to_string(temp.path().join("stdout"))?,
        fs::read_to_string(temp.path().join("stderr"))?
    );
    Ok(CheckedProof {
        authority: version,
        source_sha256: format!("{:x}", Sha256::digest(source.as_bytes())),
        scope: "Conditional propositional entailment and probability bounds only; scientific premises are not proved.".into(),
    })
}

/// Audit also consumes saved files. Restrict them to our finite proof language
/// before invoking Lean; arbitrary Lean can execute IO through evaluation commands.
fn validate_source(text: &str) -> Result<()> {
    ensure!(text.len() <= 4_000_000, "proof source exceeds size budget");
    let generated = text
        .strip_prefix(&format!("{PRELUDE}\nnamespace Helix.Generated\n"))
        .and_then(|s| s.strip_suffix("\nend Helix.Generated\n"))
        .context("proof source is not a supported Helix export")?;
    let mut lines = generated.lines().filter(|line| !line.is_empty());
    let conditional = lines.next().context("missing conditional theorem")?;
    ensure!(
        conditional.starts_with("theorem conditional_conclusion "),
        "missing conditional theorem"
    );
    ensure!(
        conditional
            .chars()
            .all(|c| c.is_ascii_alphanumeric()
                || matches!(c, '_' | ' ' | '(' | ')' | ':' | '=' | '→')),
        "unsupported conditional proof syntax"
    );
    for token in conditional
        .split(|c: char| !c.is_ascii_alphanumeric() && c != '_')
        .filter(|s| !s.is_empty())
    {
        let variable = token.as_bytes().first().is_some_and(|b| b"phr".contains(b))
            && token.len() > 1
            && token[1..].bytes().all(|b| b.is_ascii_digit());
        ensure!(
            variable || matches!(token, "theorem" | "conditional_conclusion" | "Prop"),
            "unsupported proof token: {token}"
        );
    }
    for (index, line) in lines.enumerate() {
        // Strict token positions admit only our natural-number decision theorems.
        let tokens: Vec<_> = line.split_whitespace().collect();
        ensure!(tokens.len() == 17, "unsupported arithmetic proof syntax");
        let n: u32 = tokens[3]
            .strip_prefix('(')
            .context("invalid numerator")?
            .parse()?;
        let d: u32 = tokens[7].parse()?;
        ensure!(
            line == format!(
                "theorem probability_{index} : ({n} : Nat) ≤ {d} ∧ (0 : Nat) < {d} := by decide"
            ),
            "unsupported arithmetic proof syntax"
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn missing_authority_cannot_produce_a_checked_artifact() {
        let d = Derivation::derive(vec!["observed".into()], vec![], "observed".into()).unwrap();
        let error = check(Path::new("/no/such/lean"), &source(&d, "")).unwrap_err();
        assert!(error.to_string().contains("Lean unavailable"));
    }
    #[test]
    fn exporter_contains_no_admitted_proof() {
        let d = Derivation::derive(vec!["observed".into()], vec![], "observed".into()).unwrap();
        let text = source(
            &d,
            "theorem probability_0 : (1 : Nat) ≤ 2 ∧ (0 : Nat) < 2 := by decide",
        );
        assert!(text.contains("conditional_conclusion"));
        assert!(!text.contains(":= sorry"));
        assert!(!text.contains("\naxiom "));
    }
    #[test]
    fn lean_rejects_a_false_computational_obligation() {
        let d = Derivation::derive(vec!["observed".into()], vec![], "observed".into()).unwrap();
        let error = check(
            Path::new("lean"),
            &source(
                &d,
                "theorem probability_0 : (3 : Nat) ≤ 2 ∧ (0 : Nat) < 2 := by decide",
            ),
        )
        .unwrap_err();
        assert!(error.to_string().contains("Lean rejected proof"), "{error}");
    }
    #[test]
    fn saved_proofs_cannot_add_io_commands_or_admit_results() {
        let d = Derivation::derive(vec!["observed".into()], vec![], "observed".into()).unwrap();
        for injection in [
            "#eval IO.println \"unsafe\"",
            "axiom injected : False",
            "theorem probability_0 : False := by sorry",
        ] {
            assert!(validate_source(&source(&d, injection)).is_err());
        }
        assert!(validate_source(&source(&d, "")).is_ok());
    }
}
