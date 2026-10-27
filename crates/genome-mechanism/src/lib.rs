//! Reference-grounded SNV mechanisms for a validated HBB transcript.
//!
//! Sequence changes are computed from the genomic reference and transcript model,
//! without clinical labels or precomputed consequence annotations. Molecular
//! consequences do not establish clinical pathogenicity or measured function.
use anyhow::{bail, ensure, Context, Result};
use genome_types::{Consequence, GRCh38, GenomicPosition, Variant};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs,
    path::{Component, Path},
};

const MAX_SEQUENCE_LENGTH: u64 = 100_000;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranscriptModel {
    pub accession: String,
    pub protein_accession: String,
    pub strand: String,
    pub tx_start: u64,
    pub tx_end: u64,
    pub cds_start: u64,
    pub cds_end: u64,
    pub exons: Vec<[u64; 2]>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub file: String,
    pub url: String,
    pub retrieved_at: String,
    pub sha256: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceManifest {
    pub schema_version: u32,
    pub assembly: String,
    pub gene: String,
    pub chromosome: String,
    pub sequence_start: u64,
    pub sequence_end: u64,
    pub transcript: TranscriptModel,
    pub sources: Vec<Source>,
    pub files: BTreeMap<String, String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Region {
    Coding,
    FivePrimeUtr,
    ThreePrimeUtr,
    Intron,
    Upstream,
    Downstream,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MolecularEffect {
    Synonymous,
    Missense,
    StartLost,
    StopGained,
    StopLost,
    SpliceDonor,
    SpliceAcceptor,
    Untranslated,
    Intronic,
    Upstream,
    Downstream,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SpliceSiteKind {
    Donor,
    Acceptor,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BoundarySide {
    Exonic,
    Intronic,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpliceProximity {
    pub intron_index_1based: usize,
    pub kind: SpliceSiteKind,
    pub side: BoundarySide,
    /// One is the first base on the specified side of the exon/intron boundary.
    pub distance_from_boundary_1based: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpliceAnnotation {
    /// Intron index in transcript 5' to 3' order, starting at one.
    pub intron_index_1based: usize,
    pub kind: SpliceSiteKind,
    /// Position within the two-base motif, in transcript orientation.
    pub motif_position_1based: usize,
    pub reference_motif: String,
    pub alternate_motif: String,
    pub canonical_reference: bool,
    pub disrupts_canonical_motif: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MechanisticAnnotation {
    pub variant_key: String,
    pub genomic_position_0based: u64,
    pub reference: String,
    pub alternate: String,
    pub transcript_accession: String,
    pub protein_accession: String,
    pub strand: String,
    pub consequence: Consequence,
    pub effect: MolecularEffect,
    pub region: Region,
    pub transcript_position_1based: Option<u64>,
    pub cds_position_1based: Option<u64>,
    /// Includes initiator methionine. Traditional mature beta-globin numbering
    /// differs by one; this field never silently uses mature-protein numbering.
    pub protein_position_1based: Option<u64>,
    pub reference_codon: Option<String>,
    pub alternate_codon: Option<String>,
    pub reference_amino_acid: Option<String>,
    pub alternate_amino_acid: Option<String>,
    pub splice_sites: Vec<SpliceAnnotation>,
    /// Positional flags within two bases on either side of a junction. Exonic
    /// flags retain the coding consequence and do not assert a splicing defect.
    pub splice_proximity: Vec<SpliceProximity>,
    pub notes: Vec<String>,
}

/// Opened only after integrity, exon reconstruction and protein checks succeed.
pub struct ReferenceBundle {
    manifest: ReferenceManifest,
    manifest_sha256: String,
    genomic: Vec<u8>,
    transcript: Vec<u8>,
    cds: Vec<u8>,
    protein: Vec<u8>,
    transcript_positions: BTreeMap<u64, usize>,
    cds_positions: BTreeMap<u64, usize>,
}

impl ReferenceBundle {
    pub fn open(directory: &Path) -> Result<Self> {
        let root = directory
            .canonicalize()
            .context("reference bundle directory")?;
        let raw_manifest = fs::read(root.join("manifest.json"))?;
        let manifest: ReferenceManifest = serde_json::from_slice(&raw_manifest)?;
        ensure!(
            manifest.schema_version == 1
                && manifest.assembly == "GRCh38"
                && manifest.gene == "HBB"
                && manifest.chromosome == "11",
            "unsupported reference schema, assembly, gene, or chromosome"
        );
        ensure!(
            manifest.sequence_start < manifest.sequence_end
                && manifest.sequence_end <= 135_086_622
                && manifest.sequence_end - manifest.sequence_start <= MAX_SEQUENCE_LENGTH,
            "invalid reference interval or sequence budget exceeded"
        );
        ensure!(
            manifest.transcript.accession == "NM_000518.5"
                && manifest.transcript.protein_accession == "NP_000509.1"
                && manifest.transcript.strand == "-",
            "this engine requires the versioned HBB minus-strand reference transcript"
        );
        for name in [
            "genomic.fasta",
            "transcript.fasta",
            "protein.fasta",
            "annotation.json",
        ] {
            ensure!(
                manifest.files.contains_key(name),
                "reference manifest does not pin {name}"
            );
        }
        let mut contents = BTreeMap::new();
        for (name, commitment) in &manifest.files {
            let mut components = Path::new(name).components();
            ensure!(
                matches!(components.next(), Some(Component::Normal(_)))
                    && components.next().is_none(),
                "reference file must be a basename"
            );
            let path = root.join(name).canonicalize()?;
            ensure!(path.starts_with(&root), "reference file escapes bundle");
            let bytes = fs::read(path)?;
            ensure!(
                digest(&bytes) == *commitment,
                "reference checksum mismatch for {name}"
            );
            contents.insert(name.as_str(), bytes);
        }
        ensure!(
            !manifest.sources.is_empty(),
            "reference provenance is required"
        );
        for source in &manifest.sources {
            ensure!(
                source.url.starts_with("https://")
                    && !source.retrieved_at.is_empty()
                    && manifest.files.get(&source.file) == Some(&source.sha256),
                "source metadata does not match a pinned reference file"
            );
        }
        validate_annotation_source(&manifest, &contents["annotation.json"])?;
        let (genomic_header, genomic) = fasta(&contents["genomic.fasta"], false)?;
        let (transcript_header, transcript) = fasta(&contents["transcript.fasta"], false)?;
        let (protein_header, protein) = fasta(&contents["protein.fasta"], true)?;
        let expected_genomic_id = format!(
            "NC_000011.10:{}-{}",
            manifest.sequence_start + 1,
            manifest.sequence_end
        );
        ensure!(
            genomic_header.split_whitespace().next() == Some(expected_genomic_id.as_str()),
            "genomic FASTA accession/interval mismatch"
        );
        ensure!(
            transcript_header.split_whitespace().next()
                == Some(manifest.transcript.accession.as_str()),
            "transcript FASTA accession mismatch"
        );
        ensure!(
            protein_header.split_whitespace().next()
                == Some(manifest.transcript.protein_accession.as_str()),
            "protein FASTA accession mismatch"
        );
        ensure!(
            genomic.len() as u64 == manifest.sequence_end - manifest.sequence_start,
            "genomic FASTA length mismatch"
        );
        Self::from_sequences(
            manifest,
            digest(&raw_manifest),
            genomic,
            transcript,
            protein,
        )
    }

    fn from_sequences(
        manifest: ReferenceManifest,
        manifest_sha256: String,
        genomic: Vec<u8>,
        transcript: Vec<u8>,
        protein: Vec<u8>,
    ) -> Result<Self> {
        let model = &manifest.transcript;
        ensure!(
            model.tx_start >= manifest.sequence_start
                && model.tx_end <= manifest.sequence_end
                && model.tx_start < model.tx_end
                && model.cds_start >= model.tx_start
                && model.cds_end <= model.tx_end
                && model.cds_start < model.cds_end,
            "transcript/CDS bounds outside reference"
        );
        ensure!(!model.exons.is_empty(), "no exons");
        let mut previous_end = None;
        let mut tx_positions = Vec::new();
        let mut coding_positions = Vec::new();
        for &[start, end] in &model.exons {
            ensure!(
                model.tx_start <= start && start < end && end <= model.tx_end,
                "invalid exon interval"
            );
            if let Some(previous_end) = previous_end {
                ensure!(
                    start >= previous_end + 4,
                    "exons overlap or intron lacks separate core motifs"
                );
            }
            previous_end = Some(end);
            tx_positions.extend(start..end);
            coding_positions.extend(start.max(model.cds_start)..end.min(model.cds_end));
        }
        ensure!(
            model.exons.first().unwrap()[0] == model.tx_start
                && model.exons.last().unwrap()[1] == model.tx_end,
            "terminal exons disagree with transcript bounds"
        );
        tx_positions.reverse();
        coding_positions.reverse();
        let assemble = |positions: &[u64]| -> Result<Vec<u8>> {
            positions
                .iter()
                .map(|position| {
                    let index = usize::try_from(position - manifest.sequence_start)?;
                    let base = *genomic
                        .get(index)
                        .context("reference interval shorter than annotation")?;
                    complement(base)
                })
                .collect()
        };
        ensure!(
            assemble(&tx_positions)? == transcript,
            "spliced genomic exons disagree with transcript FASTA"
        );
        let cds = assemble(&coding_positions)?;
        ensure!(
            !cds.is_empty() && cds.len().is_multiple_of(3) && cds.starts_with(b"ATG"),
            "CDS is incomplete, out of frame, or lacks canonical start"
        );
        let translated = translate(&cds)?;
        ensure!(
            translated.last() == Some(&b'*') && !translated[..translated.len() - 1].contains(&b'*'),
            "CDS must contain one terminal stop and no internal stops"
        );
        ensure!(
            translated[..translated.len() - 1] == protein,
            "translated CDS disagrees with reference protein"
        );
        let cds_first = coding_positions.first().context("empty CDS")?;
        let cds_last = coding_positions.last().context("empty CDS")?;
        let first_tx = tx_positions
            .iter()
            .position(|p| p == cds_first)
            .context("CDS start not exonic")?;
        let last_tx = tx_positions
            .iter()
            .position(|p| p == cds_last)
            .context("CDS end not exonic")?;
        ensure!(
            transcript[first_tx..=last_tx] == cds,
            "CDS is not a contiguous interval of spliced transcript"
        );
        Ok(Self {
            manifest,
            manifest_sha256,
            genomic,
            transcript,
            cds,
            protein,
            transcript_positions: tx_positions
                .into_iter()
                .enumerate()
                .map(|(i, p)| (p, i))
                .collect(),
            cds_positions: coding_positions
                .into_iter()
                .enumerate()
                .map(|(i, p)| (p, i))
                .collect(),
        })
    }

    pub fn manifest(&self) -> &ReferenceManifest {
        &self.manifest
    }
    pub fn manifest_sha256(&self) -> &str {
        &self.manifest_sha256
    }
    pub fn genomic_sequence(&self) -> &str {
        std::str::from_utf8(&self.genomic).unwrap()
    }
    pub fn transcript_sequence(&self) -> &str {
        std::str::from_utf8(&self.transcript).unwrap()
    }
    pub fn cds_sequence(&self) -> &str {
        std::str::from_utf8(&self.cds).unwrap()
    }
    pub fn protein_sequence(&self) -> &str {
        std::str::from_utf8(&self.protein).unwrap()
    }

    pub fn annotate(&self, variant: &Variant<GRCh38>) -> Result<MechanisticAnnotation> {
        ensure!(
            variant.position().chromosome() == self.manifest.chromosome,
            "variant chromosome differs from reference"
        );
        let position = variant.position().offset();
        ensure!(
            (self.manifest.sequence_start..self.manifest.sequence_end).contains(&position),
            "variant outside bundled reference interval"
        );
        ensure!(
            self.base(position)? == variant.reference() as u8,
            "variant reference allele disagrees with bundled reference"
        );
        let model = &self.manifest.transcript;
        let transcript_index = self.transcript_positions.get(&position).copied();
        let coding_index = self.cds_positions.get(&position).copied();
        let region = if coding_index.is_some() {
            Region::Coding
        } else if transcript_index.is_some() {
            if position >= model.cds_end {
                Region::FivePrimeUtr
            } else {
                Region::ThreePrimeUtr
            }
        } else if (model.tx_start..model.tx_end).contains(&position) {
            Region::Intron
        } else if position >= model.tx_end {
            Region::Upstream
        } else {
            Region::Downstream
        };
        let splice_sites = self.splice_annotations(position, variant.alternate() as u8)?;
        let splice_proximity = self.splice_proximity(position);
        let (mut consequence, mut effect) = match region {
            Region::Coding => (Consequence::Other, MolecularEffect::Missense),
            Region::FivePrimeUtr | Region::ThreePrimeUtr => {
                (Consequence::Utr, MolecularEffect::Untranslated)
            }
            Region::Intron => (Consequence::Intron, MolecularEffect::Intronic),
            Region::Upstream => (Consequence::Other, MolecularEffect::Upstream),
            Region::Downstream => (Consequence::Other, MolecularEffect::Downstream),
        };
        let mut reference_codon = None;
        let mut alternate_codon = None;
        let mut reference_amino_acid = None;
        let mut alternate_amino_acid = None;
        let mut protein_position = None;
        let mut notes = vec![
            "Sequence-derived consequence only; pathogenicity, penetrance and clinical actionability are not inferred.".into(),
            "Protein coordinates include the initiator methionine; traditional mature beta-globin numbering differs by one.".into(),
        ];
        if let Some(index) = coding_index {
            let codon_start = index / 3 * 3;
            let original = &self.cds[codon_start..codon_start + 3];
            let mut changed = original.to_vec();
            changed[index % 3] = complement(variant.alternate() as u8)?;
            let original_amino_acid = amino_acid(original)?;
            let changed_amino_acid = amino_acid(&changed)?;
            (consequence, effect) = if original_amino_acid == changed_amino_acid {
                (Consequence::Synonymous, MolecularEffect::Synonymous)
            } else if original_amino_acid == b'*' {
                (Consequence::Stop, MolecularEffect::StopLost)
            } else if changed_amino_acid == b'*' {
                (Consequence::Stop, MolecularEffect::StopGained)
            } else if codon_start == 0 {
                (Consequence::Other, MolecularEffect::StartLost)
            } else {
                (Consequence::Missense, MolecularEffect::Missense)
            };
            reference_codon = Some(String::from_utf8(original.to_vec())?);
            alternate_codon = Some(String::from_utf8(changed)?);
            reference_amino_acid = Some((original_amino_acid as char).to_string());
            alternate_amino_acid = Some((changed_amino_acid as char).to_string());
            protein_position = Some(index as u64 / 3 + 1);
            if effect == MolecularEffect::Synonymous {
                notes.push("An unchanged translated codon does not rule out effects on splicing, RNA stability or expression.".into());
            }
            if effect == MolecularEffect::StartLost {
                notes.push("Loss of the annotated ATG is recorded; alternate initiation is not modeled. The legacy coarse consequence enum has no start-loss category, so the detailed start_lost effect maps to 'other'.".into());
            }
        }
        if let Some(site) = splice_sites.first() {
            consequence = Consequence::Splice;
            effect = match site.kind {
                SpliceSiteKind::Donor => MolecularEffect::SpliceDonor,
                SpliceSiteKind::Acceptor => MolecularEffect::SpliceAcceptor,
            };
            notes.push("Only the intronic two-base donor/acceptor core is assessed. Motif disruption is not a measured splicing defect; branch points, exonic enhancers and cryptic sites are unmodeled.".into());
        }
        if region != Region::Coding && splice_sites.is_empty() {
            notes.push("No expression or splicing effect is inferred for noncoding sequence outside the assessed splice cores.".into());
        }
        if splice_proximity
            .iter()
            .any(|site| site.side == BoundarySide::Exonic)
        {
            notes.push("This variant is within two exonic bases of a splice junction. This positional flag can coexist with a coding consequence; it does not establish a splicing defect.".into());
        }
        Ok(MechanisticAnnotation {
            variant_key: variant.key(),
            genomic_position_0based: position,
            reference: variant.reference().to_string(),
            alternate: variant.alternate().to_string(),
            transcript_accession: model.accession.clone(),
            protein_accession: model.protein_accession.clone(),
            strand: model.strand.clone(),
            consequence,
            effect,
            region,
            transcript_position_1based: transcript_index.map(|p| p as u64 + 1),
            cds_position_1based: coding_index.map(|p| p as u64 + 1),
            protein_position_1based: protein_position,
            reference_codon,
            alternate_codon,
            reference_amino_acid,
            alternate_amino_acid,
            splice_sites,
            splice_proximity,
            notes,
        })
    }

    /// Enumerate the three non-reference SNVs at every position in the annotated
    /// transcript locus [tx_start, tx_end), including introns and untranslated
    /// exons. Bundled upstream/downstream flanks are excluded from this map.
    pub fn saturation_map(&self) -> Result<Vec<MechanisticAnnotation>> {
        let model = &self.manifest.transcript;
        let mut annotations = Vec::with_capacity((model.tx_end - model.tx_start) as usize * 3);
        for position in model.tx_start..model.tx_end {
            let reference = self.base(position)? as char;
            for alternate in ['A', 'C', 'G', 'T'] {
                if reference == alternate {
                    continue;
                }
                let variant = Variant::<GRCh38>::snv(
                    GenomicPosition::new(&self.manifest.chromosome, position)?,
                    &reference.to_string(),
                    &alternate.to_string(),
                )?;
                annotations.push(self.annotate(&variant)?);
            }
        }
        Ok(annotations)
    }

    fn base(&self, position: u64) -> Result<u8> {
        let offset = position
            .checked_sub(self.manifest.sequence_start)
            .context("position precedes reference")?;
        self.genomic
            .get(usize::try_from(offset)?)
            .copied()
            .context("position outside reference")
    }

    fn splice_annotations(&self, position: u64, alternate: u8) -> Result<Vec<SpliceAnnotation>> {
        let exons = &self.manifest.transcript.exons;
        let mut annotations = Vec::new();
        for (index, adjacent) in exons.windows(2).enumerate() {
            let [left, right] = adjacent else {
                unreachable!()
            };
            // Minus-strand transcription proceeds from the right exon to the left.
            let motifs = [
                (SpliceSiteKind::Donor, [right[0] - 2, right[0]], b"GT"),
                (SpliceSiteKind::Acceptor, [left[1], left[1] + 2], b"AG"),
            ];
            for (kind, [start, end], canonical) in motifs {
                if !(start..end).contains(&position) {
                    continue;
                }
                let original = vec![
                    complement(self.base(end - 1)?)?,
                    complement(self.base(start)?)?,
                ];
                let mut changed = original.clone();
                let motif_index = usize::try_from(end - 1 - position)?;
                changed[motif_index] = complement(alternate)?;
                let canonical_reference = original == canonical;
                annotations.push(SpliceAnnotation {
                    intron_index_1based: exons.len() - 1 - index,
                    kind,
                    motif_position_1based: motif_index + 1,
                    reference_motif: String::from_utf8(original)?,
                    alternate_motif: String::from_utf8(changed.clone())?,
                    canonical_reference,
                    disrupts_canonical_motif: canonical_reference && changed != canonical,
                });
            }
        }
        Ok(annotations)
    }

    fn splice_proximity(&self, position: u64) -> Vec<SpliceProximity> {
        let exons = &self.manifest.transcript.exons;
        let mut flags = Vec::new();
        for (index, adjacent) in exons.windows(2).enumerate() {
            let left = adjacent[0];
            let right = adjacent[1];
            let intervals = [
                (
                    SpliceSiteKind::Donor,
                    BoundarySide::Intronic,
                    right[0] - 2,
                    right[0],
                    true,
                ),
                (
                    SpliceSiteKind::Donor,
                    BoundarySide::Exonic,
                    right[0],
                    right[0] + 2,
                    false,
                ),
                (
                    SpliceSiteKind::Acceptor,
                    BoundarySide::Intronic,
                    left[1],
                    left[1] + 2,
                    false,
                ),
                (
                    SpliceSiteKind::Acceptor,
                    BoundarySide::Exonic,
                    left[1] - 2,
                    left[1],
                    true,
                ),
            ];
            for (kind, side, start, end, reverse) in intervals {
                if (start..end).contains(&position) {
                    flags.push(SpliceProximity {
                        intron_index_1based: exons.len() - 1 - index,
                        kind,
                        side,
                        distance_from_boundary_1based: if reverse {
                            end - position
                        } else {
                            position - start + 1
                        },
                    });
                }
            }
        }
        flags
    }
}

fn validate_annotation_source(manifest: &ReferenceManifest, bytes: &[u8]) -> Result<()> {
    let annotation: serde_json::Value = serde_json::from_slice(bytes)?;
    ensure!(
        annotation["genome"] == "hg38" && annotation["chrom"] == "chr11",
        "raw annotation assembly/chromosome mismatch"
    );
    let model = &manifest.transcript;
    let records = annotation["ncbiRefSeq"]
        .as_array()
        .context("raw annotation has no RefSeq records")?;
    let matches: Vec<_> = records
        .iter()
        .filter(|record| record["name"].as_str() == Some(model.accession.as_str()))
        .collect();
    ensure!(
        matches.len() == 1,
        "raw annotation must identify one matching transcript"
    );
    let record = matches[0];
    ensure!(
        record["name2"] == "HBB"
            && record["chrom"] == "chr11"
            && record["strand"] == "-"
            && record["txStart"].as_u64() == Some(model.tx_start)
            && record["txEnd"].as_u64() == Some(model.tx_end)
            && record["cdsStart"].as_u64() == Some(model.cds_start)
            && record["cdsEnd"].as_u64() == Some(model.cds_end)
            && record["exonCount"].as_u64() == Some(model.exons.len() as u64)
            && record["cdsStartStat"] == "cmpl"
            && record["cdsEndStat"] == "cmpl",
        "transcript manifest disagrees with raw RefSeq annotation"
    );
    let coordinates = |field: &str| -> Result<Vec<u64>> {
        record[field]
            .as_str()
            .context("invalid raw exon coordinate list")?
            .trim_end_matches(',')
            .split(',')
            .map(|s| s.parse().context("invalid exon coordinate"))
            .collect()
    };
    let starts = coordinates("exonStarts")?;
    let ends = coordinates("exonEnds")?;
    ensure!(
        starts.len() == ends.len() && starts.len() == model.exons.len(),
        "raw exon count mismatch"
    );
    let exons: Vec<_> = starts
        .into_iter()
        .zip(ends)
        .map(|(start, end)| [start, end])
        .collect();
    ensure!(
        exons == model.exons,
        "transcript exons disagree with raw RefSeq annotation"
    );
    Ok(())
}

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn fasta(bytes: &[u8], protein: bool) -> Result<(String, Vec<u8>)> {
    let text = std::str::from_utf8(bytes)?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .and_then(|line| line.strip_prefix('>'))
        .context("FASTA needs one header")?;
    let mut sequence = Vec::new();
    for line in lines {
        ensure!(
            !line.starts_with('>'),
            "multiple FASTA records are unsupported"
        );
        sequence.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
    }
    let alphabet: &[u8] = if protein {
        b"ACDEFGHIKLMNPQRSTVWY"
    } else {
        b"ACGT"
    };
    ensure!(
        !sequence.is_empty() && sequence.iter().all(|b| alphabet.contains(b)),
        "empty or ambiguous FASTA sequence"
    );
    Ok((header.to_owned(), sequence))
}

fn complement(base: u8) -> Result<u8> {
    match base {
        b'A' => Ok(b'T'),
        b'C' => Ok(b'G'),
        b'G' => Ok(b'C'),
        b'T' => Ok(b'A'),
        _ => bail!("ambiguous DNA base"),
    }
}

fn amino_acid(codon: &[u8]) -> Result<u8> {
    ensure!(codon.len() == 3, "codon must contain three bases");
    let mut index = 0;
    for base in codon {
        let value = match base {
            b'T' => 0,
            b'C' => 1,
            b'A' => 2,
            b'G' => 3,
            _ => bail!("ambiguous codon"),
        };
        index = index * 4 + value;
    }
    Ok(b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[index])
}

fn translate(sequence: &[u8]) -> Result<Vec<u8>> {
    ensure!(sequence.len().is_multiple_of(3), "partial codon");
    sequence.chunks_exact(3).map(amino_acid).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../datasets/hbb-reference")
    }
    fn bundle() -> ReferenceBundle {
        ReferenceBundle::open(&fixture()).unwrap()
    }
    fn variant(position: u64, reference: &str, alternate: &str) -> Variant<GRCh38> {
        Variant::snv(
            GenomicPosition::new("11", position).unwrap(),
            reference,
            alternate,
        )
        .unwrap()
    }

    #[test]
    fn standard_code_has_all_codons_and_correct_landmarks() {
        assert_eq!(
            translate(b"ATGGTGCACCTGACTCCTGAGGAGAAGTAA").unwrap(),
            b"MVHLTPEEK*"
        );
        assert_eq!(amino_acid(b"TGA").unwrap(), b'*');
        assert_eq!(amino_acid(b"TGG").unwrap(), b'W');
        assert!(amino_acid(b"NNN").is_err());
        assert!(translate(b"ATGA").is_err());
        let mut count = 0;
        for a in b"ACGT" {
            for b in b"ACGT" {
                for c in b"ACGT" {
                    amino_acid(&[*a, *b, *c]).unwrap();
                    count += 1;
                }
            }
        }
        assert_eq!(count, 64);
    }

    #[test]
    fn reference_reconstructs_transcript_and_protein() {
        let reference = bundle();
        assert_eq!(reference.transcript_sequence().len(), 628);
        assert_eq!(reference.cds_sequence().len(), 444);
        assert_eq!(reference.protein_sequence().len(), 147);
        assert!(reference.cds_sequence().starts_with("ATGGTGCAT"));
        assert!(reference.protein_sequence().starts_with("MVHLTPEEK"));
    }

    #[test]
    fn sickle_and_hbc_are_reverse_strand_glutamate_changes() {
        let reference = bundle();
        let sickle = reference.annotate(&variant(5_227_001, "T", "A")).unwrap();
        assert_eq!(sickle.effect, MolecularEffect::Missense);
        assert_eq!(sickle.cds_position_1based, Some(20));
        assert_eq!(sickle.protein_position_1based, Some(7));
        assert_eq!(sickle.reference_codon.as_deref(), Some("GAG"));
        assert_eq!(sickle.alternate_codon.as_deref(), Some("GTG"));
        assert_eq!(sickle.reference_amino_acid.as_deref(), Some("E"));
        assert_eq!(sickle.alternate_amino_acid.as_deref(), Some("V"));
        let hbc = reference.annotate(&variant(5_227_002, "C", "T")).unwrap();
        assert_eq!(hbc.cds_position_1based, Some(19));
        assert_eq!(hbc.alternate_codon.as_deref(), Some("AAG"));
        assert_eq!(hbc.alternate_amino_acid.as_deref(), Some("K"));
    }
    #[test]
    fn codons_split_across_exons_and_equivalent_protein_changes_are_preserved() {
        let reference = bundle();
        let terminal_exonic = reference.annotate(&variant(5_226_929, "C", "A")).unwrap();
        assert_eq!(terminal_exonic.cds_position_1based, Some(92));
        assert_eq!(terminal_exonic.protein_position_1based, Some(31));
        assert_eq!(terminal_exonic.reference_codon.as_deref(), Some("AGG"));
        assert_eq!(terminal_exonic.alternate_codon.as_deref(), Some("ATG"));
        assert_eq!(terminal_exonic.effect, MolecularEffect::Missense);
        assert_eq!(terminal_exonic.splice_proximity.len(), 1);
        assert_eq!(
            terminal_exonic.splice_proximity[0].side,
            BoundarySide::Exonic
        );
        assert_eq!(
            terminal_exonic.splice_proximity[0].distance_from_boundary_1based,
            1
        );
        let one = reference.annotate(&variant(5_225_714, "G", "T")).unwrap();
        let two = reference.annotate(&variant(5_225_714, "G", "C")).unwrap();
        assert_eq!(one.protein_position_1based, Some(109));
        assert_eq!(one.reference_codon.as_deref(), Some("AAC"));
        assert_eq!(one.alternate_codon.as_deref(), Some("AAA"));
        assert_eq!(two.alternate_codon.as_deref(), Some("AAG"));
        assert_eq!(one.alternate_amino_acid.as_deref(), Some("K"));
        assert_eq!(one.alternate_amino_acid, two.alternate_amino_acid);
    }
    #[test]
    fn clinical_disagreement_variant_is_still_a_missense_sequence_change() {
        let reference = bundle();
        let annotation = reference.annotate(&variant(5_226_658, "T", "A")).unwrap();
        assert_eq!(annotation.cds_position_1based, Some(233));
        assert_eq!(annotation.protein_position_1based, Some(78));
        assert_eq!(annotation.reference_amino_acid.as_deref(), Some("H"));
        assert_eq!(annotation.alternate_amino_acid.as_deref(), Some("L"));
        assert_eq!(annotation.effect, MolecularEffect::Missense);
    }

    #[test]
    fn canonical_minus_strand_donor_and_acceptor_are_flagged() {
        let reference = bundle();
        let exons = &reference.manifest.transcript.exons;
        for adjacent in exons.windows(2) {
            for position in [
                adjacent[1][0] - 1,
                adjacent[1][0] - 2,
                adjacent[0][1],
                adjacent[0][1] + 1,
            ] {
                let base = reference.base(position).unwrap() as char;
                let alternate = if base == 'A' { "C" } else { "A" };
                let annotation = reference
                    .annotate(&variant(position, &base.to_string(), alternate))
                    .unwrap();
                assert_eq!(annotation.consequence, Consequence::Splice);
                assert_eq!(annotation.splice_sites.len(), 1);
                assert!(annotation.splice_sites[0].canonical_reference);
                assert!(annotation.splice_sites[0].disrupts_canonical_motif);
                assert!(matches!(
                    annotation.splice_sites[0].reference_motif.as_str(),
                    "GT" | "AG"
                ));
            }
        }
    }

    #[test]
    fn saturation_is_complete_unique_and_does_not_read_clinical_labels() {
        let reference = bundle();
        let annotations = reference.saturation_map().unwrap();
        assert_eq!(
            annotations.len(),
            (reference.manifest.transcript.tx_end - reference.manifest.transcript.tx_start)
                as usize
                * 3
        );
        let keys: std::collections::BTreeSet<_> =
            annotations.iter().map(|a| &a.variant_key).collect();
        assert_eq!(keys.len(), annotations.len());
        assert_eq!(
            annotations
                .iter()
                .filter(|a| a.region == Region::Coding)
                .count(),
            444 * 3
        );
        assert_eq!(
            annotations
                .iter()
                .filter(|a| !a.splice_sites.is_empty())
                .count(),
            2 * 2 * 2 * 3
        );
        assert!(annotations
            .iter()
            .any(|a| a.effect == MolecularEffect::Synonymous));
        assert!(annotations
            .iter()
            .any(|a| a.effect == MolecularEffect::StartLost));
        assert!(annotations
            .iter()
            .any(|a| a.effect == MolecularEffect::StopGained));
        assert!(annotations
            .iter()
            .any(|a| a.effect == MolecularEffect::StopLost));
        assert!(annotations.iter().any(|a| a.region == Region::FivePrimeUtr));
        assert!(annotations
            .iter()
            .any(|a| a.region == Region::ThreePrimeUtr));
        assert!(!annotations.iter().any(|a| a.region == Region::Upstream));
        assert!(!annotations.iter().any(|a| a.region == Region::Downstream));
        for (position, expected) in [
            (reference.manifest.transcript.tx_end, Region::Upstream),
            (
                reference.manifest.transcript.tx_start - 1,
                Region::Downstream,
            ),
        ] {
            let base = reference.base(position).unwrap() as char;
            let alternate = if base == 'A' { "C" } else { "A" };
            assert_eq!(
                reference
                    .annotate(&variant(position, &base.to_string(), alternate))
                    .unwrap()
                    .region,
                expected
            );
        }
    }

    #[test]
    fn wrong_reference_and_unsupported_coordinates_fail() {
        let reference = bundle();
        assert!(reference.annotate(&variant(5_227_001, "C", "A")).is_err());
        assert!(reference.annotate(&variant(1, "A", "T")).is_err());
        let other =
            Variant::<GRCh38>::snv(GenomicPosition::new("12", 5_227_001).unwrap(), "T", "A")
                .unwrap();
        assert!(reference.annotate(&other).is_err());
    }

    #[test]
    fn corrupted_bundle_cannot_be_loaded() {
        let temp = tempfile::tempdir().unwrap();
        for entry in fs::read_dir(fixture()).unwrap() {
            let entry = entry.unwrap();
            if entry.file_type().unwrap().is_file() {
                fs::copy(entry.path(), temp.path().join(entry.file_name())).unwrap();
            }
        }
        fs::write(temp.path().join("protein.fasta"), ">NP_000509.1\nAAAA\n").unwrap();
        assert!(ReferenceBundle::open(temp.path())
            .err()
            .unwrap()
            .to_string()
            .contains("checksum mismatch"));
    }
}
