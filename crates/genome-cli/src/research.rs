//! Integrate independently implemented reference, atlas, evidence and statistics.
use anyhow::{ensure, Context, Result};
use genome_assay::PromoterAssay;
use genome_atlas::Atlas;
use genome_data::{sha256, Dataset};
use genome_inference::{ConsequenceModel, Predictor, PrevalenceBaseline};
use genome_knowledge::{Ledger, Revision};
use genome_logic::{Derivation, Rule};
use genome_mechanism::{MechanisticAnnotation, ReferenceBundle, Region};
use genome_research_stats::{agreement_summary, paired_comparison};
use genome_types::{ClinicalLabel, Consequence, GRCh38, GenomicPosition, Probability, Variant};
use serde::Serialize;
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
    process::Command,
};

#[derive(clap::Args)]
pub struct ResearchArgs {
    #[arg(long, default_value = "datasets/clinvar-hbb")]
    dataset: PathBuf,
    #[arg(long, default_value = "datasets/hbb-reference")]
    reference: PathBuf,
    #[arg(long, default_value = "datasets/genome-atlas")]
    atlas: PathBuf,
    #[arg(long, default_value = "datasets/hbb-promoter-assay")]
    assay: PathBuf,
    #[arg(long, default_value = "research/hbb-hypotheses.json")]
    hypotheses: PathBuf,
    #[arg(long, default_value = "experiments/hbb-mechanism-protocol.md")]
    protocol: PathBuf,
    #[arg(long, default_value = "experiments/hbb-promoter-assay-protocol.md")]
    assay_protocol: PathBuf,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, default_value = "lean")]
    lean: PathBuf,
}

const FILES: [&str; 15] = [
    "annotation-comparison.json",
    "controls.json",
    "evidence-audit.json",
    "genome-atlas.json",
    "hypotheses.json",
    "hypothesis-outcomes.json",
    "ledger.json",
    "promoter-assay.json",
    "promoter-clinvar-overlap.json",
    "proof.lean",
    "report.json",
    "research.html",
    "research.md",
    "saturation.json",
    "statistics.json",
];

#[derive(Serialize)]
struct Comparison {
    id: String,
    supplied_consequence: Consequence,
    recomputed: MechanisticAnnotation,
    agrees_at_coarse_resolution: bool,
    source_binary_classification: ClinicalLabel,
    partition: &'static str,
    interpretation: &'static str,
}

pub fn run(args: &ResearchArgs) -> Result<()> {
    ensure!(
        !args.out.exists(),
        "output already exists; choose a new research run directory"
    );
    let dataset = Dataset::open(&args.dataset)?;
    let reference = ReferenceBundle::open(&args.reference)?;
    let atlas = Atlas::open(&args.atlas)?;
    let assay = PromoterAssay::open(&args.assay, &reference)?;
    let evidence = genome_evidence_audit::audit_dataset(&args.dataset)?;
    let hypothesis_bytes =
        fs::read(&args.hypotheses).context("read structured research hypotheses")?;
    let hypotheses: Value = serde_json::from_slice(&hypothesis_bytes)?;
    let protocol = fs::read(&args.protocol).context("read research protocol")?;
    let assay_protocol = fs::read(&args.assay_protocol).context("read reporter-assay protocol")?;
    ensure!(
        !protocol.is_empty()
            && !assay_protocol.is_empty()
            && hypotheses["schema_version"] == 1
            && hypotheses["hypotheses"]
                .as_array()
                .is_some_and(|items| items.len() >= 3),
        "empty protocol or invalid hypothesis collection"
    );
    let expected_hypotheses: BTreeSet<_> = [
        "HBB-REF-01",
        "HBB-ANNOT-02",
        "HBB-EQUIV-03",
        "HBB-JUNCTION-04",
        "HBB-COVERAGE-05",
    ]
    .into_iter()
    .collect();
    let hypothesis_items = hypotheses["hypotheses"].as_array().unwrap();
    let actual_ids: BTreeSet<_> = hypothesis_items
        .iter()
        .filter_map(|h| h["id"].as_str())
        .collect();
    ensure!(
        hypothesis_items.len() == expected_hypotheses.len() && actual_ids == expected_hypotheses,
        "this finite research protocol requires exactly its five defined hypothesis IDs"
    );
    let sources_commitment = verify_literature_sources()?;
    let reviewed: Value =
        serde_json::from_slice(&fs::read("research/sources/uniprot-P68871.json")?)?;
    ensure!(
        reviewed["sequence"]["value"].as_str() == Some(reference.protein_sequence()),
        "reviewed UniProt sequence disagrees with reference protein"
    );
    ensure!(
        reference.protein_sequence().as_bytes()[63] == b'H'
            && reference.protein_sequence().as_bytes()[92] == b'H',
        "named heme-interaction residues disagree with reference"
    );
    let mut ledger = Ledger::default();
    let commitments = json!({
        "clinvar_manifest_sha256": dataset.manifest_sha256,
        "reference_manifest_sha256": reference.manifest_sha256(),
        "atlas_manifest_sha256": atlas.manifest_sha256(),
        "promoter_assay_manifest_sha256": assay.manifest_sha256(),
        "hypotheses_sha256": sha256(&hypothesis_bytes),
        "protocol_sha256": sha256(&protocol),
        "assay_protocol_sha256": sha256(&assay_protocol),
        "literature_source_manifest_sha256": sources_commitment
    });
    ledger.append("set_research_goal", json!({
        "goal": "Reconstruct HBB DNA-to-RNA-to-protein consequences, enumerate its local SNV mechanism space, reconcile curated assertions, and place the locus in a primary-human annotation atlas",
        "scope": "Computational reference and evidence analysis; not full human-genome understanding or independent biological discovery",
        "input_commitments": commitments
    }))?;
    let model = reference.manifest();
    let locus_length = model.transcript.tx_end - model.transcript.tx_start;
    ensure!(
        atlas.report()["hbb"]["start"].as_u64() == Some(model.transcript.tx_start)
            && atlas.report()["hbb"]["end"].as_u64() == Some(model.transcript.tx_end),
        "atlas HBB interval disagrees with the reference transcript protocol"
    );
    let sequence_checks = json!({
        "genomic_reference_bases": reference.genomic_sequence().len(),
        "transcript_bases": reference.transcript_sequence().len(),
        "coding_bases_including_stop": reference.cds_sequence().len(),
        "protein_residues_excluding_stop_including_initiator_methionine": reference.protein_sequence().len(),
        "transcript_locus_bases": locus_length,
        "exon_intervals_zero_based_half_open": model.transcript.exons,
        "coding_bounds_zero_based_half_open": [model.transcript.cds_start,model.transcript.cds_end],
        "exact_transcript_reconstruction": true,
        "exact_reference_protein_translation": true,
        "reviewed_uniprot_sequence_matches_refseq": true,
        "reference_transcript": model.transcript.accession,
        "reference_protein": model.transcript.protein_accession,
        "strand": model.transcript.strand
    });
    ensure!(reference.transcript_sequence().len() == 628 && reference.cds_sequence().len() == 444
        && reference.protein_sequence().len() == 147 && locus_length == 1608,
        "this protocol and its finite Lean facts require the frozen NM_000518.5 reference quantities");
    ledger.append("validate_reference_reconstruction", &sequence_checks)?;
    let controls = controls(&reference)?;
    ledger.append("validate_known_reference_controls", &controls)?;
    let train_ids: BTreeSet<_> = dataset
        .training
        .iter()
        .map(|e| e.feature.id.as_str())
        .collect();
    let all_labels: BTreeMap<_, _> = dataset
        .training
        .iter()
        .map(|e| (e.observation.id.as_str(), e.observation.label))
        .chain(dataset.heldout().iter().map(|l| (l.id.as_str(), l.label)))
        .collect();
    let source_assertions: BTreeMap<_, _> = dataset
        .training
        .iter()
        .map(|e| (e.observation.id.as_str(), &e.observation))
        .chain(dataset.heldout().iter().map(|l| (l.id.as_str(), l)))
        .collect();
    let comparisons: Vec<_> = dataset.features.iter().map(|feature| {
        let recomputed = reference.annotate(&feature.variant)?;
        Ok(Comparison {
            id: feature.id.clone(),
            supplied_consequence: feature.consequence,
            agrees_at_coarse_resolution: recomputed.consequence == feature.consequence,
            recomputed,
            source_binary_classification: all_labels[feature.id.as_str()],
            partition: if train_ids.contains(feature.id.as_str()) { "training" } else { "previously_exposed_heldout" },
            interpretation: "Reference-derived consequence is compared with existing coarse ClinVar annotation, not measured function or clinical truth",
        })
    }).collect::<Result<_>>()?;
    let agreements = comparisons
        .iter()
        .filter(|c| c.agrees_at_coarse_resolution)
        .count();
    let agreement = agreement_summary(dataset.features.len(), comparisons.len(), agreements)?;
    let coding_rows: Vec<_> = comparisons
        .iter()
        .filter(|c| c.recomputed.region == Region::Coding)
        .collect();
    let coding_agreements = coding_rows
        .iter()
        .filter(|c| c.agrees_at_coarse_resolution)
        .count();
    let coding_agreement =
        agreement_summary(coding_rows.len(), coding_rows.len(), coding_agreements)?;
    ledger.append("compare_reference_derived_annotations", &agreement)?;
    let saturation = reference.saturation_map()?;
    ensure!(
        saturation.len() as u64 == locus_length * 3,
        "incomplete saturation enumeration"
    );
    let unique_keys: BTreeSet<_> = saturation.iter().map(|a| a.variant_key.as_str()).collect();
    ensure!(
        unique_keys.len() == saturation.len(),
        "duplicate saturation variants"
    );
    let mut effect_counts = BTreeMap::<String, usize>::new();
    let mut region_counts = BTreeMap::<String, usize>::new();
    for annotation in &saturation {
        *effect_counts
            .entry(enum_name(&annotation.effect)?)
            .or_default() += 1;
        *region_counts
            .entry(enum_name(&annotation.region)?)
            .or_default() += 1;
    }
    let coding_substitutions = saturation
        .iter()
        .filter(|a| a.cds_position_1based.is_some())
        .count();
    let canonical_splice_substitutions = saturation
        .iter()
        .filter(|a| !a.splice_sites.is_empty())
        .count();
    let observed_locus_variants = comparisons
        .iter()
        .filter(|c| {
            (model.transcript.tx_start..model.transcript.tx_end)
                .contains(&c.recomputed.genomic_position_0based)
        })
        .count();
    let observed_keys: BTreeSet<_> = comparisons
        .iter()
        .map(|c| c.recomputed.variant_key.as_str())
        .collect();
    let saturation_summary = json!({
        "interval_zero_based_half_open": [model.transcript.tx_start, model.transcript.tx_end],
        "interval_scope": "Transcript genomic span including introns; excludes flanking promoter/enhancers",
        "all_possible_single_base_substitutions": saturation.len(),
        "unique_variant_keys": unique_keys.len(),
        "coding_substitutions": coding_substitutions,
        "canonical_splice_core_substitutions": canonical_splice_substitutions,
        "effect_counts": effect_counts, "region_counts": region_counts,
        "selected_clinvar_variants_within_locus": observed_locus_variants,
        "selected_clinvar_fraction_of_enumerated_variants": observed_locus_variants as f64 / saturation.len() as f64,
        "known_control_presence_in_selected_clinvar": controls.as_array().unwrap().iter().map(|c| json!({
            "name":c["name"], "present":observed_keys.contains(c["annotation"]["variant_key"].as_str().unwrap())
        })).collect::<Vec<_>>(),
        "heme_interaction_positions": ([64,93].iter().map(|&position| json!({
            "protein_position_including_initiator_methionine":position,
            "reference_residue": (reference.protein_sequence().as_bytes()[position-1] as char).to_string(),
            "selected_variant_count": comparisons.iter().filter(|c| c.recomputed.protein_position_1based == Some(position as u64)).count(),
            "role":if position == 64 {"distal histidine stabilizes bound ligand"} else {"proximal histidine directly coordinates heme iron"},
            "source":"Reviewed UniProt P68871"
        })).collect::<Vec<_>>()),
        "interpretation": "Exhaustive reference-level possibilities, not predicted pathogenic variants. Absence from this selected ClinVar sample is not absence from biology or ClinVar."
    });
    ledger.append("enumerate_local_snv_mechanisms", &saturation_summary)?;
    let promoter_overlap: Vec<_> = assay.records().iter().filter_map(|measured| {
        comparisons.iter().find(|c|c.recomputed.variant_key == measured.variant_key).map(|clinical|json!({
            "variant_key":measured.variant_key,
            "source_clinvar_accession":clinical.id,
            "curated_binary_clinical_classification":clinical.source_binary_classification,
            "source_clinvar_assertion":source_assertions[clinical.id.as_str()],
            "recomputed_reference_annotation":clinical.recomputed,
            "measured_reporter_record":measured,
            "interpretation":"Exact allele identity links the records, not their endpoints. Reporter coefficients are not clinical probabilities; curation may have used this experiment, so independent confirmation is not established."
        }))
    }).collect();
    let assay_payload = json!({"summary":assay.report(),"records":assay.records(),"source_rows":assay.source_rows()});
    ledger.append("integrate_measured_reporter_endpoint", assay.report())?;
    ledger.append(
        "link_assay_alleles_without_equating_endpoints",
        &promoter_overlap,
    )?;
    let consequence_model = ConsequenceModel::fit(&dataset.training)?;
    let baseline = PrevalenceBaseline::fit(&dataset.training)?;
    let targets = dataset.heldout_features();
    let model_predictions = targets
        .iter()
        .map(|f| consequence_model.predict(f))
        .collect::<Result<Vec<_>>>()?;
    let baseline_predictions = targets
        .iter()
        .map(|f| baseline.predict(f))
        .collect::<Result<Vec<_>>>()?;
    let labels: BTreeMap<_, _> = dataset
        .heldout()
        .iter()
        .map(|l| (l.id.as_str(), l.label))
        .collect();
    let correctness = |predictions: &[genome_inference::Prediction]| -> Vec<bool> {
        predictions
            .iter()
            .map(|p| p.probability_pathogenic.label() == labels[p.id.as_str()])
            .collect()
    };
    ensure!(
        model_predictions
            .iter()
            .zip(&baseline_predictions)
            .all(|(a, b)| a.id == b.id),
        "paired prediction identities differ"
    );
    let paired = paired_comparison(
        &correctness(&model_predictions),
        &correctness(&baseline_predictions),
    )?;
    let mut category_labels = BTreeMap::<Consequence, [usize; 2]>::new();
    for row in &comparisons {
        category_labels.entry(row.supplied_consequence).or_default()
            [usize::from(row.source_binary_classification.positive())] += 1;
    }
    let maximum_fit: usize = category_labels
        .values()
        .map(|counts| counts[0].max(counts[1]))
        .sum();
    let statistics = json!({
        "reference_annotation_agreement": agreement,
        "coding_only_annotation_agreement": coding_agreement,
        "previous_classification_benchmark_post_hoc": paired,
        "coarse_category_information_limit":{
            "category_label_counts_benign_then_pathogenic":category_labels,
            "best_possible_in_sample_fit_for_a_deterministic_category_only_classifier":maximum_fit,
            "selected_records":comparisons.len(),
            "scope":"Conditional finite-sample upper bound using all already-inspected labels, not prospective performance or a biological limit"
        },
        "sample_scope":"Existing heldout has already been inspected. Intervals and p-values are descriptive under unestablished sampling assumptions; no superiority or discovery claim."
    });
    ledger.append("audit_prior_classification_result_post_hoc", &statistics)?;
    ledger.append("audit_evidence_strength_and_selection", &evidence)?;
    ledger.append("map_primary_genome_context", atlas.report())?;
    ledger.append("record_falsifiable_questions", &hypotheses)?;
    let equivalent: Vec<_> = comparisons
        .iter()
        .filter(|c| {
            c.recomputed.protein_position_1based == Some(109)
                && c.recomputed.alternate_amino_acid.as_deref() == Some("K")
        })
        .collect();
    let junction: Vec<_> = comparisons
        .iter()
        .filter(|c| c.recomputed.cds_position_1based == Some(92))
        .collect();
    let hypothesis_outcomes = json!([
        {"id":"HBB-REF-01","status":"reference_consistency_checks_passed","result":sequence_checks,
            "controls":controls,"biological_validation":false},
        {"id":"HBB-ANNOT-02","status":"descriptive_annotation_comparison_completed",
            "all_regions":agreement,"coding_only":coding_agreement,
            "disagreement_ids":comparisons.iter().filter(|c|!c.agrees_at_coarse_resolution).map(|c|&c.id).collect::<Vec<_>>(),
            "downstream_rna_and_function":"No comprehensive matched RNA/protein assay coverage; limited promoter-reporter overlaps are recorded separately"},
        {"id":"HBB-EQUIV-03","status":if equivalent.len()>=2 {"narrow_translation_equivalence_observed"} else {"not_covered_by_selected_records"},
            "alleles":equivalent,"claim_scope":"Same translated amino-acid substitution only; distinct DNA identities, no inferred RNA equivalence or replication"},
        {"id":"HBB-JUNCTION-04","status":"computational_boundary_context_checked_biological_extension_untested",
            "alleles":junction,"additional_splicing_effect_established":false},
        {"id":"HBB-COVERAGE-05","status":"reference_defined_coverage_audit_completed",
            "result":saturation_summary,"absence_is_not_benign":true}
    ]);
    ledger.append(
        "evaluate_computable_hypothesis_components",
        &hypothesis_outcomes,
    )?;
    let derivation = Derivation::derive(
        vec!["The frozen reference, transcript model and protein sequence are correct and describe the selected isoform".into()],
        vec![Rule {
            premise: "The frozen reference, transcript model and protein sequence are correct and describe the selected isoform".into(),
            conclusion: "Reference-level codon changes are a constrained description of sequence consequences, not a sufficient clinical model".into(),
            assumption: "The annotation algorithm follows the standard genetic code and correctly handles strand and exon boundaries; this assumption is tested computationally, not proved by the Lean theorem".into(),
        }],
        "Reference-level codon changes are a constrained description of sequence consequences, not a sufficient clinical model".into()
    )?;
    let ratios = [
        Probability::new(
            u32::try_from(agreements)?,
            u32::try_from(comparisons.len())?,
        )?,
        Probability::new(
            u32::try_from(observed_locus_variants)?,
            u32::try_from(saturation.len())?,
        )?,
        Probability::new(
            u32::try_from(coding_substitutions)?,
            u32::try_from(saturation.len())?,
        )?,
    ];
    let proof_source = genome_proof::mechanism_source(
        &derivation,
        &genome_logic::probability_obligations(&ratios),
    );
    let proof = genome_proof::check_mechanisms(&args.lean, &proof_source)?;
    ledger.append("check_scoped_formal_properties", &proof)?;
    ledger.append("research_update", json!({
        "independent_biological_validation": false,
        "evidence_transition": "No predicted mechanism is promoted to measured or replicated evidence",
        "external_reporter_evidence_integrated": true,
        "next_goal": "Freeze a regulatory-prediction task with held-out experimental units and additional contexts; no fitted predictor or superiority claim is made in this descriptive run"
    }))?;
    Ledger::verify(ledger.revisions())?;
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
        "schema_version":2, "workflow":"reference-grounded-hbb-research-v1",
        "goal":"Understand how a defined human locus converts sequence into transcript/protein consequences and where that mechanism map falls short of clinical explanation",
        "software":{"version":env!("CARGO_PKG_VERSION"),"commit":commit,"working_tree_dirty":dirty},
        "input_commitments":commitments, "reference_validation":sequence_checks,
        "annotation_comparison":agreement, "saturation":saturation_summary,
        "coding_annotation_comparison":coding_agreement,
        "statistics":statistics, "evidence_audit":evidence,
        "promoter_assay":assay.report(), "promoter_clinvar_overlap":promoter_overlap,
        "genome_context":atlas.report(), "hypotheses":hypotheses,
        "hypothesis_outcomes":hypothesis_outcomes,
        "reasoning":derivation, "proof":{"status":"checked","artifact":proof},
        "scientific_status":{
            "reference_computations":"validated against frozen RefSeq sequences",
            "annotation_consistency":"compared with existing ClinVar coarse consequences; not an independent functional endpoint",
            "clinical_analysis":"post hoc on previously exposed selected records",
            "measured_regulatory_endpoint":"Source-estimated log2 reporter-expression coefficients in HEL 92.1.7 plasmid MPRA; not endogenous HBB or clinical risk",
            "independent_validation_of_toolkit_functional_predictions":false,
            "independent_replication":false, "new_biological_discovery":false,
            "whole_human_genome_reverse_engineered":false
        },
        "limitations":[
            "One HBB transcript; local SNV consequences do not model all isoforms, structural variants, regulatory networks or the rest of the genome",
            "Canonical splice motifs and proximity are features; unchanged motifs do not exclude cryptic-splice or regulatory effects",
            "Amino-acid substitution does not establish stability, oxygen affinity, polymerization, expression or disease severity",
            "Reporter coefficients are source regression estimates from a pooled mutagenized library; culture replicates are not independent studies and coefficient p-values have source-acknowledged inflation limits",
            "Noncoding and synonymous annotations do not imply absence of function",
            "No clinical recommendations; no experimental intervention or patient data",
            "Lean checks small definitions, arithmetic and conditional entailment, not the complete Rust implementation or empirical premises",
            "Genome atlas is a source-specific curated annotation and omits some gene classes; uncovered sequence is not shown to be nonfunctional"
        ]
    });
    let narrative = narrative(&report)?;
    let html = include_str!("research.html").replace(
        "__DATA_JSON__",
        &serde_json::to_string(&json!({
            "report":report, "comparisons":comparisons, "controls":controls
        }))?
        .replace('<', "\\u003c")
        .replace('>', "\\u003e")
        .replace('&', "\\u0026"),
    );
    let payloads = [
        ("annotation-comparison.json", super::pretty(&comparisons)?),
        ("controls.json", super::pretty(&controls)?),
        ("evidence-audit.json", super::pretty(&evidence)?),
        ("genome-atlas.json", super::pretty(atlas.report())?),
        ("hypotheses.json", super::pretty(&hypotheses)?),
        (
            "hypothesis-outcomes.json",
            super::pretty(&hypothesis_outcomes)?,
        ),
        ("ledger.json", super::pretty(&ledger.revisions())?),
        ("promoter-assay.json", super::pretty(&assay_payload)?),
        (
            "promoter-clinvar-overlap.json",
            super::pretty(&promoter_overlap)?,
        ),
        ("proof.lean", proof_source.into_bytes()),
        ("report.json", super::pretty(&report)?),
        ("research.html", html.into_bytes()),
        ("research.md", narrative.into_bytes()),
        ("saturation.json", super::pretty(&saturation)?),
        ("statistics.json", super::pretty(&statistics)?),
    ];
    if let Some(parent) = args.out.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    fs::create_dir(&args.out)?;
    let mut hashes = BTreeMap::new();
    for (name, bytes) in payloads {
        super::write_new(&args.out.join(name), &bytes)?;
        hashes.insert(name.into(), sha256(&bytes));
    }
    super::write_new(
        &args.out.join("artifacts.json"),
        &super::pretty(&super::Artifacts {
            schema_version: 2,
            files: hashes,
        })?,
    )?;
    println!("Reference verified: 628 transcript bases, 444 coding bases, 147 protein residues.");
    println!(
        "Enumerated {} SNVs; compared {} existing records, {} coarse-annotation agreements.",
        saturation.len(),
        comparisons.len(),
        agreements
    );
    println!("Prior model comparison is post hoc: exact paired p = {:.3}; independent biological validation not established.", paired.two_sided_exact_p_value);
    println!("Integrated {} source-estimated reporter SNV effects in HEL 92.1.7; {} exact allele overlaps with clinical assertions (different endpoints).",
        assay.records().len(), promoter_overlap.len());
    println!(
        "Genome context: {} curated gene loci across primary chromosomes.",
        atlas.report()["gene_loci"]
    );
    println!("Saved research narrative, interactive explorer, complete mechanism map and audit artifacts to {}", args.out.display());
    Ok(())
}

fn enum_name(value: &impl Serialize) -> Result<String> {
    serde_json::to_value(value)?
        .as_str()
        .map(str::to_owned)
        .context("expected a string enum")
}
fn controls(reference: &ReferenceBundle) -> Result<Value> {
    let mut result = Vec::new();
    for (name, offset, base, alternate, expected) in [
        ("HbS", 5_227_001, "T", "A", "V"),
        ("HbC", 5_227_002, "C", "T", "K"),
    ] {
        let variant = Variant::<GRCh38>::snv(GenomicPosition::new("11", offset)?, base, alternate)?;
        let annotation = reference.annotate(&variant)?;
        ensure!(
            annotation.protein_position_1based == Some(7)
                && annotation.reference_amino_acid.as_deref() == Some("E")
                && annotation.alternate_amino_acid.as_deref() == Some(expected),
            "known HBB control disagrees with sequence reconstruction"
        );
        result.push(json!({
            "name":name, "annotation":annotation,
            "numbering":"RefSeq includes initiator methionine (Glu7); traditional mature beta-globin uses Glu6",
            "source":"Reviewed UniProt P68871 and reference-crosscheck.json; established controls, not new discoveries"
        }));
    }
    Ok(json!(result))
}
fn verify_literature_sources() -> Result<String> {
    let root = Path::new("research/sources")
        .canonicalize()
        .context("literature source directory")?;
    let bytes = fs::read(root.join("manifest.json"))?;
    let manifest: Value = serde_json::from_slice(&bytes)?;
    ensure!(
        manifest["schema_version"] == 1,
        "unsupported literature manifest"
    );
    let sources = manifest["files"]
        .as_object()
        .context("invalid literature source commitments")?;
    ensure!(!sources.is_empty(), "no literature sources");
    for (name, commitment) in sources {
        ensure!(
            Path::new(name).components().count() == 1,
            "literature paths must be basenames"
        );
        let path = root.join(name).canonicalize()?;
        ensure!(path.starts_with(&root), "literature source escapes root");
        ensure!(
            Some(sha256(&fs::read(path)?).as_str()) == commitment.as_str(),
            "literature source checksum mismatch: {name}"
        );
    }
    let retrievals = manifest["retrievals"]
        .as_array()
        .context("missing literature retrieval provenance")?;
    ensure!(!retrievals.is_empty(), "no literature retrieval metadata");
    for retrieval in retrievals {
        ensure!(
            retrieval["url"]
                .as_str()
                .is_some_and(|s| s.starts_with("https://")),
            "invalid literature retrieval URL"
        );
        if let Some(retained) = retrieval["retained_in"]
            .as_str()
            .or_else(|| retrieval["file"].as_str())
        {
            ensure!(
                sources.contains_key(retained),
                "retrieval points to uncommitted representation"
            );
        } else {
            ensure!(retrieval["retained_representation"].as_str().is_some_and(|s|s.contains("not redistributed"))
                && retrieval["original_response_sha256"].as_str().is_some_and(|s|s.len()==64),
                "unretained literature retrieval needs an explicit retention limit and original commitment");
        }
    }
    Ok(sha256(&bytes))
}

fn narrative(report: &Value) -> Result<String> {
    let paired = &report["statistics"]["previous_classification_benchmark_post_hoc"];
    let atlas = &report["genome_context"];
    Ok(format!(
        "# HBB mechanism research run\n\n\
         Goal: reconstruct a reference-grounded DNA-to-RNA-to-protein mechanism map, audit its evidence, and place the locus in a source-specific human-genome context.\n\n\
         ## What was computed\n\n\
         - NM_000518.5: 628-base transcript; 444-base CDS including stop; NP_000509.1: 147 residues including initiator methionine.\n\
         - All {} single-base substitutions over the 1608-base transcript locus, including introns. {} are coding and {} alter canonical splice-core positions.\n\
         - {} / {} supplied coarse ClinVar consequences agree with sequence-derived annotations. Agreement checks annotation consistency, not measured function.\n\
         - Primary-chromosome atlas: {} curated gene loci, including {} protein-coding loci, from a frozen RefSeq annotation. This is not the full gene repertoire or functional genome.\n\n\
         ## Measured regulatory endpoint\n\n\
         The separately produced Kircher2019 / MaveDB HBB fragment provides {} reference-validated SNV effect estimates across 187 target bases, plus {} literal equality rows retained outside SNV summaries. Its target maps to the plus genomic strand; native HBB is minus-strand. The fragment contains 50 bases of transcribed 5-prime UTR and 137 upstream bases.\n\n\
         Scores are source-estimated log2 reporter-expression coefficients in HEL 92.1.7 plasmid MPRA, not endogenous HBB abundance or clinical risk. Random library constructs can contain several changes; the additive-regression assumptions and source-acknowledged p-value inflation remain explicit. Three cultures use the same library, not three independent studies. {} exact alleles overlap the clinical fixture; promoter-clinvar-overlap.json links their identities without equating endpoints or claiming independent clinical confirmation.\n\n\
         ## Evidence changes the interpretation\n\n\
         The old clinical benchmark is already exposed and this analysis is post hoc. Model-only correct: {}; baseline-only correct: {}; exact paired p-value: {:.3}. There is no established superiority claim. Selected heldout review-status counts: {}. Disjoint loci can share cohorts and nearby sequence. See evidence-audit.json for missing disease contexts and selection exclusions.\n\n\
         ## Mechanistic limits\n\n\
         HbS and HbC controls verify distinct sequence/codon changes. Protein positions include initiator methionine: Glu7 corresponds to the traditional mature-chain Glu6. Sequence changes alone do not measure oxygen affinity, polymerization, stability, expression, splicing or disease severity. Unchanged canonical splice motifs do not rule out cryptic-splice effects. Synonymous, intronic and unannotated do not mean nonfunctional.\n\n\
         ## What the proofs establish\n\n\
         Lean checks finite complement definitions, coordinate/arithmetic identities, ratio bounds and entailment under stated assumptions. It does not verify the Rust implementation or empirical premises. No prediction is promoted to measured or replicated evidence.\n\n\
         ## Next falsifiable work\n\n\
         Review hypotheses.json, hypothesis-outcomes.json and the source-backed research notes. Freeze a prediction task with meaningful held-out units and new cellular/experimental contexts before fitting or claiming discovery. Current results are a mechanism map and measured-evidence inventory, not reverse engineering of the complete human genome.\n\n\
         Open research.html for the local interactive explorer. Full per-variant outputs are annotation-comparison.json and saturation.json. All consumed source manifests and the protocol are committed by input hashes in report.json; artifact hashes and the revision ledger are verifiable with the CLI audit command.\n",
        report["saturation"]["all_possible_single_base_substitutions"], report["saturation"]["coding_substitutions"],
        report["saturation"]["canonical_splice_core_substitutions"], report["annotation_comparison"]["agreements"],
        report["annotation_comparison"]["compared_targets"], atlas["gene_loci"], atlas["protein_coding_gene_loci"],
        report["promoter_assay"]["rows"]["reference_validated_snv_rows"], report["promoter_assay"]["rows"]["source_equality_rows"],
        report["promoter_clinvar_overlap"].as_array().context("missing assay overlap list")?.len(),
        paired["left_only_correct"], paired["right_only_correct"], paired["two_sided_exact_p_value"].as_f64().context("missing paired statistic")?,
        report["evidence_audit"]["heldout"]["review_status_counts"]
    ))
}

pub fn audit(root: &Path, artifacts: &super::Artifacts, lean: &Path) -> Result<()> {
    ensure!(
        artifacts
            .files
            .keys()
            .map(String::as_str)
            .collect::<Vec<_>>()
            == FILES,
        "unsupported/incomplete research artifact manifest"
    );
    for (name, commitment) in &artifacts.files {
        let path = root.join(name).canonicalize()?;
        ensure!(
            path.starts_with(root),
            "research artifact escapes directory"
        );
        ensure!(
            sha256(&fs::read(path)?) == *commitment,
            "research artifact checksum mismatch: {name}"
        );
    }
    let revisions: Vec<Revision> = serde_json::from_slice(&fs::read(root.join("ledger.json"))?)?;
    Ledger::verify(&revisions)?;
    let report: Value = serde_json::from_slice(&fs::read(root.join("report.json"))?)?;
    ensure!(
        report["schema_version"] == 2
            && report["workflow"] == "reference-grounded-hbb-research-v1"
            && report["proof"]["status"] == "checked",
        "unsupported research report or missing proof"
    );
    let proof =
        genome_proof::check_mechanisms(lean, &fs::read_to_string(root.join("proof.lean"))?)?;
    ensure!(
        report["proof"]["artifact"]["source_sha256"].as_str() == Some(proof.source_sha256.as_str()),
        "research proof commitment differs"
    );
    let goal = revisions
        .iter()
        .find(|r| r.kind == "set_research_goal")
        .context("missing research goal commitment")?;
    ensure!(
        goal.payload["input_commitments"] == report["input_commitments"],
        "report input commitments differ from ledger"
    );
    let formal = revisions
        .iter()
        .find(|r| r.kind == "check_scoped_formal_properties")
        .context("missing formal-check revision")?;
    ensure!(
        formal.payload == report["proof"]["artifact"],
        "report proof artifact differs from ledger"
    );
    for (field, event) in [
        ("reference_validation", "validate_reference_reconstruction"),
        (
            "annotation_comparison",
            "compare_reference_derived_annotations",
        ),
        ("saturation", "enumerate_local_snv_mechanisms"),
        ("promoter_assay", "integrate_measured_reporter_endpoint"),
        (
            "promoter_clinvar_overlap",
            "link_assay_alleles_without_equating_endpoints",
        ),
        ("statistics", "audit_prior_classification_result_post_hoc"),
        ("evidence_audit", "audit_evidence_strength_and_selection"),
        ("genome_context", "map_primary_genome_context"),
        ("hypotheses", "record_falsifiable_questions"),
        (
            "hypothesis_outcomes",
            "evaluate_computable_hypothesis_components",
        ),
    ] {
        let matches: Vec<_> = revisions.iter().filter(|r| r.kind == event).collect();
        ensure!(
            matches.len() == 1 && matches[0].payload == report[field],
            "report {field} differs from its unique ledger event"
        );
    }
    for (filename, field) in [
        ("evidence-audit.json", "evidence_audit"),
        ("genome-atlas.json", "genome_context"),
        ("hypotheses.json", "hypotheses"),
        ("hypothesis-outcomes.json", "hypothesis_outcomes"),
        ("statistics.json", "statistics"),
        ("promoter-clinvar-overlap.json", "promoter_clinvar_overlap"),
    ] {
        let payload: Value = serde_json::from_slice(&fs::read(root.join(filename))?)?;
        ensure!(
            payload == report[field],
            "artifact {filename} differs from report"
        );
    }
    let assay: Value = serde_json::from_slice(&fs::read(root.join("promoter-assay.json"))?)?;
    ensure!(
        assay["summary"] == report["promoter_assay"],
        "assay artifact summary differs from report"
    );
    ensure!(
        report["coding_annotation_comparison"]
            == report["statistics"]["coding_only_annotation_agreement"],
        "coding-only comparison differs across report sections"
    );
    println!(
        "Research artifact hashes, {} ledger revisions, report bindings, and scoped Lean proof verified.",
        revisions.len()
    );
    println!(
        "Integrity and conditional reasoning do not establish independent biological validation."
    );
    Ok(())
}
