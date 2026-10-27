use anyhow::{ensure, Context, Result};
use clap::{Parser, Subcommand};
use genome_agents::{OfflineAgent, ResearchAgent, ResearchInput};
use genome_data::{sha256, Dataset};
use genome_evaluation::evaluate;
use genome_inference::{Predictor, PrevalenceBaseline};
use genome_knowledge::{Ledger, Revision};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{
    collections::BTreeMap,
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::Command,
};
mod research;
mod study;
mod study_context;

#[derive(Parser)]
#[command(
    name = "helix",
    version,
    about = "An auditable computational genomics research laboratory"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}
#[derive(Subcommand)]
enum Commands {
    /// Execute the frozen multi-promoter regulatory prediction study.
    Study(study::StudyArgs),
    /// Reconstruct HBB mechanisms, audit evidence, and map primary-genome context.
    Research(research::ResearchArgs),
    /// Validate the frozen public fixture and execute a bounded offline research loop.
    Demo {
        #[arg(long, default_value = "datasets/clinvar-hbb")]
        dataset: PathBuf,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value = "lean")]
        lean: PathBuf,
        /// Explicitly skip formal checking; report remains marked unverified.
        #[arg(long)]
        skip_proof: bool,
        #[arg(long, default_value_t = 1000)]
        max_steps: usize,
    },
    /// Verify artifact integrity, revision ordering, and the saved Lean proof.
    Audit {
        directory: PathBuf,
        #[arg(long, default_value = "lean")]
        lean: PathBuf,
    },
}
#[derive(Serialize, Deserialize)]
struct Artifacts {
    schema_version: u32,
    files: BTreeMap<String, String>,
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Commands::Study(args) => study::run(&args),
        Commands::Research(args) => research::run(&args),
        Commands::Demo {
            dataset,
            out,
            lean,
            skip_proof,
            max_steps,
        } => demo(&dataset, &out, &lean, skip_proof, max_steps),
        Commands::Audit { directory, lean } => audit(&directory, &lean),
    }
}
fn demo(
    directory: &Path,
    out: &Path,
    lean: &Path,
    skip_proof: bool,
    max_steps: usize,
) -> Result<()> {
    ensure!(
        !out.exists(),
        "output already exists; use a new run directory"
    );
    let dataset = Dataset::open(directory)?;
    let targets = dataset.heldout_features();
    let mut research = OfflineAgent.run(ResearchInput {
        training: &dataset.training,
        targets: &targets,
        source_sha256: &dataset.manifest.files["source.json"],
        max_steps,
    })?;
    let baseline = PrevalenceBaseline::fit(&dataset.training)?;
    let baseline_predictions = targets
        .iter()
        .map(|f| baseline.predict(f))
        .collect::<Result<Vec<_>>>()?;
    // Freeze both models before any held-out label is supplied to the evaluator.
    research
        .ledger
        .append("freeze_baseline_before_evaluation", &baseline_predictions)?;
    let probabilities: Vec<_> = research
        .predictions
        .iter()
        .chain(&baseline_predictions)
        .map(|p| p.probability_pathogenic)
        .collect();
    let proof_source = genome_proof::source(
        &research.derivation,
        &genome_logic::probability_obligations(&probabilities),
    );
    let proof = if skip_proof {
        None
    } else {
        Some(genome_proof::check(lean, &proof_source)?)
    };
    research.ledger.append(
        "formal_verification",
        json!({ "checked": proof.is_some(), "artifact": proof }),
    )?;
    let evaluation = evaluate(&research.predictions, dataset.heldout())?;
    let baseline_evaluation = evaluate(&baseline_predictions, dataset.heldout())?;
    research.ledger.append(
        "evaluate_withheld_labels",
        json!({ "model": evaluation, "baseline": baseline_evaluation }),
    )?;
    research.ledger.append("update", json!({
        "accepted_as_biological_truth": false,
        "action": "retain disagreements as unresolved; freeze independent future validation before updating the model",
        "classification_disagreements": evaluation.disagreements.len(),
        "next_questions": research.questions
    }))?;
    Ledger::verify(research.ledger.revisions())?;
    let commit = Command::new("git")
        .args(["rev-parse", "--verify", "HEAD"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_owned());
    let dirty = Command::new("git")
        .args(["status", "--porcelain"])
        .output()
        .ok()
        .map(|o| !o.status.success() || !o.stdout.is_empty());
    let report = json!({
        "schema_version": 1,
        "software": {"version": env!("CARGO_PKG_VERSION"), "commit": commit, "working_tree_dirty": dirty},
        "dataset_manifest_sha256": dataset.manifest_sha256,
        "dataset": dataset.manifest,
        "workflow": "offline deterministic agent; no LLM, private data, or clinical advice",
        "training_count": dataset.training.len(), "heldout_count": targets.len(),
        "model": research.model,
        "evaluation": evaluation, "prevalence_baseline": baseline_evaluation,
        "questions": research.questions, "contradictions": research.contradictions,
        "reasoning": research.derivation,
        "proof": {"status": if proof.is_some() {"checked"} else {"explicitly_skipped"}, "artifact": proof},
        "verification_categories": {
            "rust_contracts": "assembly-indexed SNVs, reference validation, checked source metadata",
            "formal": "implication under stated assumptions and rational probability bounds",
            "statistical": "descriptive accuracy, balanced accuracy, Brier, log-loss, five-bin calibration",
            "independent_replication": "not established; partitions share ClinVar curation",
            "biological_validation": "not established; classification is not a functional assay"
        },
        "limitations": [
            "Small, selected HBB demonstration; no generalization claim or standard genomic model comparison",
            "Disjoint loci prevent direct overlap, not shared submitter/cohort leakage or curation bias",
            "Binary aggregation combines likely and definitive classifications and excludes uncertain/conflicting records",
            "Consequence annotations are existing knowledge, not a discovered biological relationship",
            "No computed GA4GH identifiers, indel normalization, assembly liftover, or causal conclusions",
            "Ledger and hashes detect changes relative to commitments; artifacts are not authenticated or signed"
        ]
    });
    if let Some(parent) = out.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    fs::create_dir(out).context("create a new immutable run directory")?;
    let files = [
        ("report.json", pretty(&report)?),
        (
            "predictions.json",
            pretty(&json!({"model": research.predictions, "baseline": baseline_predictions}))?,
        ),
        ("evidence.json", pretty(&research.evidence)?),
        ("ledger.json", pretty(&research.ledger.revisions())?),
        ("proof.lean", proof_source.into_bytes()),
    ];
    let mut hashes = BTreeMap::new();
    for (name, bytes) in files {
        write_new(&out.join(name), &bytes)?;
        hashes.insert(name.to_owned(), sha256(&bytes));
    }
    write_new(
        &out.join("artifacts.json"),
        &pretty(&Artifacts {
            schema_version: 1,
            files: hashes,
        })?,
    )?;
    println!(
        "Validated {} training / {} held-out public HBB records.",
        dataset.training.len(),
        targets.len()
    );
    println!("Consequence model: accuracy {:.3}, Brier {:.4}; prevalence baseline: accuracy {:.3}, Brier {:.4}.",
        evaluation.accuracy, evaluation.brier_score, baseline_evaluation.accuracy, baseline_evaluation.brier_score);
    println!(
        "Formal check: {}. Independent biological replication: not established.",
        if proof.is_some() {
            "Lean checked conditional reasoning and arithmetic"
        } else {
            "explicitly skipped"
        }
    );
    println!("Saved auditable artifacts to {}", out.display());
    Ok(())
}
fn pretty(value: &impl Serialize) -> Result<Vec<u8>> {
    let mut bytes = serde_json::to_vec_pretty(value)?;
    bytes.push(b'\n');
    Ok(bytes)
}
fn write_new(path: &Path, bytes: &[u8]) -> Result<()> {
    let mut file = fs::OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(path)?;
    file.write_all(bytes)?;
    file.sync_all()?;
    Ok(())
}
fn audit(directory: &Path, lean: &Path) -> Result<()> {
    let root = directory.canonicalize()?;
    let artifacts: Artifacts = serde_json::from_slice(&fs::read(root.join("artifacts.json"))?)?;
    if artifacts.schema_version == 3 {
        return study::audit(&root, &artifacts, lean);
    }
    if artifacts.schema_version == 2 {
        return research::audit(&root, &artifacts, lean);
    }
    let required = [
        "evidence.json",
        "ledger.json",
        "predictions.json",
        "proof.lean",
        "report.json",
    ];
    ensure!(
        artifacts.schema_version == 1
            && artifacts
                .files
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>()
                == required,
        "unsupported/incomplete artifact manifest"
    );
    for (name, expected) in &artifacts.files {
        let path = root.join(name).canonicalize()?;
        ensure!(
            path.starts_with(&root),
            "artifact path escapes run directory"
        );
        ensure!(
            sha256(&fs::read(path)?) == *expected,
            "artifact checksum mismatch: {name}"
        );
    }
    let revisions: Vec<Revision> = serde_json::from_slice(&fs::read(root.join("ledger.json"))?)?;
    Ledger::verify(&revisions)?;
    let report: serde_json::Value = serde_json::from_slice(&fs::read(root.join("report.json"))?)?;
    match report["proof"]["status"].as_str() {
        Some("checked") => {
            let proof = genome_proof::check(lean, &fs::read_to_string(root.join("proof.lean"))?)?;
            ensure!(report["proof"]["artifact"]["source_sha256"].as_str() == Some(proof.source_sha256.as_str()), "proof commitment differs from report");
            println!("Artifact hashes, {} ledger revisions, and Lean proof verified.", revisions.len());
        }
        Some("explicitly_skipped") => println!("Artifact hashes and {} ledger revisions verified. Formal proof was explicitly skipped.", revisions.len()),
        _ => anyhow::bail!("unknown proof verification state"),
    }
    println!("Integrity checks do not authenticate the run or establish biological validity.");
    Ok(())
}
