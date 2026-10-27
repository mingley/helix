//! A hash-pinned, primary-chromosome RefSeq context atlas.
//!
//! Annotation coverage is never interpreted as the complete functional genome.
use anyhow::{ensure, Context, Result};
use genome_types::{Assembly, GRCh38};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::Path,
};

const RNA_BIOTYPES: [&str; 16] = [
    "RNase_MRP_RNA",
    "RNase_P_RNA",
    "Y_RNA",
    "antisense_RNA",
    "lncRNA",
    "miRNA",
    "misc_RNA",
    "ncRNA",
    "rRNA",
    "scRNA",
    "snRNA",
    "snoRNA",
    "telomerase_RNA",
    "vault_RNA",
    "tRNA",
    "SRP_RNA",
];
const SOURCE_BASE: &str = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/";
const MITO_URL: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=gb&retmode=xml";
const SOURCE_NAMES: [&str; 7] = [
    "mitochondrial-reference.xml",
    "ncbiRefSeqCurated.txt.gz",
    "ncbiRefSeqCurated.sql",
    "ncbiRefSeqLink.txt.gz",
    "ncbiRefSeqLink.sql",
    "chromInfo.txt.gz",
    "chromInfo.sql",
];
fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn primary() -> Vec<String> {
    (1..=22)
        .map(|n| format!("chr{n}"))
        .chain(["chrX".into(), "chrY".into(), "chrM".into()])
        .collect()
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub url: String,
    pub last_modified: Option<String>,
    pub bytes: u64,
    pub sha256: String,
    pub file: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    pub schema_version: u32,
    pub assembly: String,
    pub annotation: String,
    pub coordinate_system: String,
    pub retrieved_at: String,
    pub sources: Vec<Source>,
    pub primary_chromosomes: Vec<String>,
    pub files: BTreeMap<String, String>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Chromosome {
    pub name: String,
    pub length: u64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Gene {
    pub gene_id: String,
    pub symbol: String,
    pub biotype: String,
    pub chromosome: String,
    pub strand: String,
    pub start: u64,
    pub end: u64,
    pub transcript_count: u64,
    pub transcript_intervals: Vec<[u64; 2]>,
    pub exon_intervals: Vec<[u64; 2]>,
    pub cds_intervals: Vec<[u64; 2]>,
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Catalog {
    chromosomes: Vec<Chromosome>,
    genes: Vec<Gene>,
    retained_transcript_count: u64,
    excluded_transcript_count: u64,
    excluded_contig_count: usize,
    excluded_contigs: BTreeMap<String, u64>,
}

pub struct Atlas {
    pub manifest: Manifest,
    manifest_sha256: String,
    catalog: Catalog,
    report: Value,
}
impl Atlas {
    pub fn open(directory: &Path) -> Result<Self> {
        let root = directory.canonicalize().context("genome atlas directory")?;
        let manifest_bytes = fs::read(root.join("manifest.json"))?;
        let manifest: Manifest = serde_json::from_slice(&manifest_bytes)?;
        ensure!(
            manifest.schema_version == 1
                && manifest.assembly == "GRCh38/hg38"
                && manifest.annotation == "ncbiRefSeqCurated",
            "unsupported atlas schema, assembly or annotation"
        );
        ensure!(
            manifest.coordinate_system == "zero-based half-open"
                && manifest.primary_chromosomes == primary(),
            "unsupported coordinate system or chromosome set"
        );
        let expected: BTreeSet<_> = SOURCE_NAMES
            .into_iter()
            .chain(["genes.json", "summary.json"])
            .map(str::to_owned)
            .collect();
        ensure!(
            manifest.files.keys().cloned().collect::<BTreeSet<_>>() == expected,
            "atlas manifest must pin exactly the supported inputs"
        );
        let mut bytes = BTreeMap::new();
        for (name, expected_hash) in &manifest.files {
            let path = root.join(name).canonicalize()?;
            ensure!(
                path.starts_with(&root),
                "atlas source escapes dataset directory"
            );
            let payload = fs::read(path)?;
            ensure!(
                sha256(&payload) == *expected_hash,
                "atlas checksum mismatch: {name}"
            );
            bytes.insert(name.as_str(), payload);
        }
        ensure!(
            manifest.sources.len() == SOURCE_NAMES.len(),
            "incomplete atlas source metadata"
        );
        let mut seen_sources = BTreeSet::new();
        for source in &manifest.sources {
            ensure!(
                SOURCE_NAMES.contains(&source.file.as_str()) && seen_sources.insert(&source.file),
                "unsupported or duplicate atlas source"
            );
            ensure!(
                source.url
                    == if source.file == "mitochondrial-reference.xml" {
                        MITO_URL.to_owned()
                    } else {
                        format!("{SOURCE_BASE}{}", source.file)
                    },
                "atlas source URL does not match authoritative table"
            );
            ensure!(
                source.bytes == bytes[source.file.as_str()].len() as u64
                    && manifest.files[&source.file] == source.sha256,
                "atlas source commitment metadata mismatch"
            );
        }
        let catalog: Catalog = serde_json::from_slice(&bytes["genes.json"])?;
        validate_catalog(&catalog)?;
        let report = summarize(&catalog)?;
        let saved: Value = serde_json::from_slice(&bytes["summary.json"])?;
        ensure!(
            report == saved,
            "atlas summary differs from normalized catalog"
        );
        Ok(Self {
            manifest,
            manifest_sha256: sha256(&manifest_bytes),
            catalog,
            report,
        })
    }
    pub fn report(&self) -> &Value {
        &self.report
    }
    pub fn genes(&self) -> &[Gene] {
        &self.catalog.genes
    }
    pub fn chromosomes(&self) -> &[Chromosome] {
        &self.catalog.chromosomes
    }
    pub fn manifest_sha256(&self) -> &str {
        &self.manifest_sha256
    }
}

/// Return a canonical union; overlapping isoforms contribute each base once.
pub fn merge_intervals(intervals: &[[u64; 2]]) -> Result<Vec<[u64; 2]>> {
    let mut sorted = intervals.to_vec();
    sorted.sort_unstable();
    let mut result: Vec<[u64; 2]> = Vec::new();
    for [start, end] in sorted {
        ensure!(start < end, "invalid half-open interval");
        if let Some(previous) = result.last_mut().filter(|previous| start <= previous[1]) {
            previous[1] = previous[1].max(end);
        } else {
            result.push([start, end]);
        }
    }
    Ok(result)
}
fn coverage(intervals: &[[u64; 2]]) -> Result<u64> {
    merge_intervals(intervals)?
        .into_iter()
        .try_fold(0u64, |sum, [start, end]| {
            sum.checked_add(end - start).context("coverage overflow")
        })
}
fn contains(parent: &[[u64; 2]], interval: [u64; 2]) -> bool {
    parent
        .iter()
        .any(|p| p[0] <= interval[0] && interval[1] <= p[1])
}
fn validate_catalog(catalog: &Catalog) -> Result<()> {
    ensure!(
        catalog
            .chromosomes
            .iter()
            .map(|c| c.name.clone())
            .collect::<Vec<_>>()
            == primary(),
        "atlas primary chromosome set/order mismatch"
    );
    let sizes: BTreeMap<_, _> = catalog
        .chromosomes
        .iter()
        .map(|c| (c.name.as_str(), c.length))
        .collect();
    ensure!(
        catalog
            .chromosomes
            .iter()
            .zip(GRCh38::CHROMOSOME_LENGTHS)
            .all(|(chromosome, expected)| chromosome.length == expected),
        "GRCh38 primary chromosome lengths mismatch"
    );
    ensure!(
        sizes["chr11"] == 135_086_622 && sizes["chrM"] == 16569,
        "GRCh38 reference lengths mismatch"
    );
    ensure!(
        !catalog.genes.is_empty(),
        "empty atlas is not a genome context map"
    );
    let mut identities = BTreeSet::new();
    let mut transcript_count = 0u64;
    for gene in &catalog.genes {
        let size = sizes
            .get(gene.chromosome.as_str())
            .context("non-primary gene chromosome")?;
        ensure!(
            gene.start < gene.end && gene.end <= *size,
            "gene interval outside chromosome"
        );
        ensure!(
            matches!(gene.strand.as_str(), "+" | "-")
                && !gene.symbol.trim().is_empty()
                && !gene.biotype.trim().is_empty(),
            "missing gene metadata"
        );
        ensure!(
            !gene.gene_id.is_empty() && gene.gene_id.bytes().all(|b| b.is_ascii_digit()),
            "invalid NCBI Gene identity"
        );
        ensure!(
            identities.insert((&gene.chromosome, &gene.gene_id, &gene.strand)),
            "duplicate gene/chromosome/strand group"
        );
        ensure!(
            gene.transcript_count > 0
                && !gene.transcript_intervals.is_empty()
                && !gene.exon_intervals.is_empty(),
            "gene has no curated transcript or exon"
        );
        for intervals in [
            &gene.transcript_intervals,
            &gene.exon_intervals,
            &gene.cds_intervals,
        ] {
            ensure!(
                merge_intervals(intervals)? == *intervals,
                "gene intervals are not canonical merged unions"
            );
            ensure!(
                intervals
                    .iter()
                    .all(|[start, end]| gene.start <= *start && *end <= gene.end),
                "interval outside gene envelope"
            );
        }
        ensure!(
            gene.transcript_intervals.first().unwrap()[0] == gene.start
                && gene.transcript_intervals.last().unwrap()[1] == gene.end,
            "gene envelope differs from transcript extents"
        );
        ensure!(
            gene.exon_intervals
                .iter()
                .all(|i| contains(&gene.transcript_intervals, *i)),
            "exon outside transcript coverage"
        );
        ensure!(
            gene.cds_intervals
                .iter()
                .all(|i| contains(&gene.exon_intervals, *i)),
            "CDS outside exon coverage"
        );
        transcript_count = transcript_count
            .checked_add(gene.transcript_count)
            .context("transcript count overflow")?;
    }
    ensure!(
        transcript_count == catalog.retained_transcript_count,
        "retained transcript count mismatch"
    );
    ensure!(
        catalog.excluded_contig_count == catalog.excluded_contigs.len()
            && catalog
                .excluded_contigs
                .keys()
                .all(|c| !sizes.contains_key(c.as_str())),
        "invalid excluded contig metadata"
    );
    ensure!(
        catalog.excluded_contigs.values().sum::<u64>() == catalog.excluded_transcript_count,
        "excluded transcript count mismatch"
    );
    Ok(())
}
fn summarize(catalog: &Catalog) -> Result<Value> {
    let mut chromosomes = Vec::new();
    let mut total_length = 0u64;
    let (mut total_span, mut total_exon, mut total_cds) = (0u64, 0u64, 0u64);
    for chromosome in &catalog.chromosomes {
        let genes: Vec<_> = catalog
            .genes
            .iter()
            .filter(|g| g.chromosome == chromosome.name)
            .collect();
        let spans: Vec<_> = genes
            .iter()
            .flat_map(|g| g.transcript_intervals.iter().copied())
            .collect();
        let exons: Vec<_> = genes
            .iter()
            .flat_map(|g| g.exon_intervals.iter().copied())
            .collect();
        let cds: Vec<_> = genes
            .iter()
            .flat_map(|g| g.cds_intervals.iter().copied())
            .collect();
        let (span_coverage, exon_coverage, cds_coverage) =
            (coverage(&spans)?, coverage(&exons)?, coverage(&cds)?);
        ensure!(
            cds_coverage <= exon_coverage
                && exon_coverage <= span_coverage
                && span_coverage <= chromosome.length,
            "invalid coverage nesting"
        );
        total_length += chromosome.length;
        total_span += span_coverage;
        total_exon += exon_coverage;
        total_cds += cds_coverage;
        chromosomes.push(json!({"name": chromosome.name, "length": chromosome.length, "gene_loci": genes.len(),
            "transcript_count": genes.iter().map(|g| g.transcript_count).sum::<u64>(),
            "transcript_span_union_bases": span_coverage, "exon_union_bases": exon_coverage, "cds_union_bases": cds_coverage}));
    }
    let hbb: Vec<_> = catalog
        .genes
        .iter()
        .filter(|g| g.gene_id == "3043")
        .collect();
    ensure!(hbb.len() == 1, "HBB locus is not unique");
    let hbb = hbb[0];
    ensure!(
        hbb.symbol == "HBB" && hbb.chromosome == "chr11" && hbb.strand == "-",
        "HBB context mismatch"
    );
    let gap = |gene: &Gene| {
        hbb.start
            .saturating_sub(gene.end)
            .max(gene.start.saturating_sub(hbb.end))
    };
    let mut neighbors: Vec<_> = catalog
        .genes
        .iter()
        .filter(|g| g.chromosome == "chr11" && g.gene_id != "3043")
        .collect();
    neighbors.sort_by_key(|gene| (gap(gene), gene.start, &gene.gene_id));
    let neighbors: Vec<_> = neighbors.into_iter().take(12).map(|gene| json!({
        "gene_id":gene.gene_id, "symbol":gene.symbol, "biotype":gene.biotype,
        "start":gene.start, "end":gene.end, "strand":gene.strand, "gap_bases":gap(gene),
        "genomic_side":if gap(gene)==0 {"overlapping"} else if gene.end <= hbb.start {"lower_coordinate"} else {"higher_coordinate"}
    })).collect();
    let mut biotypes = BTreeMap::<&str, usize>::new();
    for gene in &catalog.genes {
        *biotypes.entry(gene.biotype.as_str()).or_default() += 1;
    }
    let coding = catalog
        .genes
        .iter()
        .filter(|g| g.biotype == "protein_coding")
        .count();
    Ok(json!({
        "schema_version":1, "assembly":"GRCh38/hg38", "annotation":"UCSC ncbiRefSeqCurated with ncbiRefSeqLink, frozen snapshot",
        "scope":"Primary chromosomes chr1-22, chrX, chrY, chrM only; curated RefSeq annotation is not a complete gene repertoire or functional map.",
        "gene_loci":catalog.genes.len(), "distinct_gene_ids":catalog.genes.iter().map(|g| &g.gene_id).collect::<BTreeSet<_>>().len(),
        "protein_coding_gene_loci":coding, "other_biotype_gene_loci":catalog.genes.len()-coding,
        "noncoding_rna_gene_loci":catalog.genes.iter().filter(|g| RNA_BIOTYPES.contains(&g.biotype.as_str())).count(),
        "pseudogene_loci":catalog.genes.iter().filter(|g| g.biotype.contains("pseudogene")).count(),
        "unclassified_biotype_gene_loci":catalog.genes.iter().filter(|g| g.biotype != "protein_coding" && !RNA_BIOTYPES.contains(&g.biotype.as_str()) && !g.biotype.contains("pseudogene")).count(),
        "gene_biotype_counts":biotypes,
        "retained_transcript_count":catalog.retained_transcript_count, "excluded_transcript_count":catalog.excluded_transcript_count,
        "excluded_contig_count":catalog.excluded_contig_count,
        "totals":{"length":total_length,"transcript_span_union_bases":total_span,"exon_union_bases":total_exon,"cds_union_bases":total_cds},
        "chromosomes":chromosomes, "hbb":hbb, "hbb_neighbors":neighbors,
        "interpretation":[
            "Intervals are zero-based, half-open; coverage is the length of merged unions, never a sum over overlapping isoforms.",
            "Gene loci group NCBI Gene ID, chromosome and strand; genes in X/Y pseudoautosomal regions can occupy more than one locus.",
            "Transcript spans include introns; exon coverage includes coding and untranslated exons; CDS coverage is coding parts of annotated exons.",
            "Other biotypes include noncoding RNAs and pseudogenes; the exact biotype counts are reported separately.",
            "Sequence outside these annotations is unannotated by this dataset, not demonstrated nonfunctional.",
            "Physical proximity and gene annotation alone do not establish regulatory interactions or disease mechanisms.",
            "The source table includes 13 mitochondrial coding models, represented by protein accessions; their GeneIDs come from NC_012920.1 CDS qualifiers. It omits mitochondrial RNA genes."
        ]
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../datasets/genome-atlas")
    }
    #[test]
    fn overlap_is_counted_once_and_invalid_intervals_fail() {
        assert_eq!(
            merge_intervals(&[[4, 7], [1, 5], [7, 9], [20, 21], [1, 3]]).unwrap(),
            vec![[1, 9], [20, 21]]
        );
        assert_eq!(coverage(&[[0, 10], [2, 5], [9, 12]]).unwrap(), 12);
        assert_eq!(coverage(&[]).unwrap(), 0);
        assert!(coverage(&[[10, 10]]).is_err());
    }
    #[test]
    fn frozen_primary_atlas_recomputes_and_identifies_hbb() {
        let atlas = Atlas::open(&fixture()).unwrap();
        assert_eq!(atlas.chromosomes().len(), 25);
        assert_eq!(atlas.report()["hbb"]["gene_id"], "3043");
        assert!(atlas.genes().len() > 20_000);
        assert_eq!(atlas.manifest_sha256().len(), 64);
    }
    #[test]
    fn checksum_mismatch_fails_before_catalog_interpretation() {
        let temp = tempfile::tempdir().unwrap();
        let root = fixture();
        for entry in fs::read_dir(root).unwrap() {
            let entry = entry.unwrap();
            if entry.file_type().unwrap().is_file() {
                fs::copy(entry.path(), temp.path().join(entry.file_name())).unwrap();
            }
        }
        fs::write(temp.path().join("genes.json"), "{}").unwrap();
        assert!(Atlas::open(temp.path())
            .err()
            .unwrap()
            .to_string()
            .contains("atlas checksum mismatch"));
    }
}
