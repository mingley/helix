//! Inspect an explicitly selected reference with deterministic JSON output.
use anyhow::{ensure, Result};
use genome_atlas::{composition::AtlasComposition, Atlas};
use genome_data::twobit::{BaseCounts, ByteOrder, Contig, Strand, TwoBitReference};
use serde::Serialize;
use std::{
    io::{self, BufWriter, Write},
    path::PathBuf,
};

// The twoBit format has no assembly label. This digest identifies the
// reference fetched by scripts/setup_reference.sh when no label is supplied.
const FROZEN_HG38_SHA256: &str = "1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0";

#[derive(clap::Args)]
pub struct ReferenceArgs {
    /// Local UCSC twoBit v0 file.
    #[arg(long)]
    file: PathBuf,
    /// Assembly label declared by the caller.
    #[arg(long, value_enum)]
    assembly: Option<AssemblyLabel>,
    /// Verify the file against this SHA-256 before producing output.
    #[arg(long, value_parser = parse_sha256)]
    expect_sha256: Option<String>,
    #[command(subcommand)]
    command: ReferenceCommand,
}

#[derive(Clone, Copy, Serialize, clap::ValueEnum)]
enum AssemblyLabel {
    #[serde(rename = "GRCh37")]
    #[value(name = "GRCh37", alias = "grch37")]
    Grch37,
    #[serde(rename = "GRCh38")]
    #[value(name = "GRCh38", alias = "grch38")]
    Grch38,
}

#[derive(Clone, Copy, Serialize, clap::ValueEnum)]
enum Orientation {
    #[serde(rename = "+")]
    #[value(name = "+", alias = "forward")]
    Forward,
    #[serde(rename = "-")]
    #[value(name = "-", alias = "reverse")]
    Reverse,
}

impl From<Orientation> for Strand {
    fn from(orientation: Orientation) -> Self {
        match orientation {
            Orientation::Forward => Self::Forward,
            Orientation::Reverse => Self::Reverse,
        }
    }
}

#[derive(clap::Subcommand)]
enum ReferenceCommand {
    /// List literal contig names, lengths, and the number of N and mask blocks.
    Catalog,
    /// Decode a zero-based half-open interval.
    Fetch {
        contig: String,
        start0: u64,
        end0: u64,
        /// Sequence orientation.
        #[arg(long, value_enum, default_value = "+", allow_hyphen_values = true)]
        strand: Orientation,
        /// Preserve the reference's soft masks as lowercase bases.
        #[arg(long)]
        masked: bool,
    },
    /// Count A/C/G/T/N across every indexed contig.
    Summary,
    /// Count bases by curated coding, exonic and transcript coverage on the GRCh38 primaries.
    AtlasComposition {
        /// Directory containing the verified RefSeq atlas.
        #[arg(long, default_value = "datasets/genome-atlas")]
        atlas: PathBuf,
    },
}

#[derive(Serialize)]
struct ReferenceFile {
    format: &'static str,
    format_version: u32,
    file_bytes: u64,
    byte_order: ByteOrder,
    declared_assembly: Option<AssemblyLabel>,
    sha256: Option<String>,
}

#[derive(Serialize)]
struct Report<'a> {
    schema_version: u32,
    reference: ReferenceFile,
    #[serde(flatten)]
    result: ReferenceResult<'a>,
}

#[derive(Serialize)]
#[serde(tag = "operation", rename_all = "snake_case")]
enum ReferenceResult<'a> {
    Catalog {
        contigs: Vec<Contig<'a>>,
    },
    Fetch {
        contig: &'a str,
        start0: u64,
        end0: u64,
        strand: Orientation,
        masked: bool,
        sequence: String,
    },
    Summary {
        contigs: Vec<ContigCounts>,
        total: BaseCounts,
    },
    AtlasComposition(Box<AtlasComposition>),
}

#[derive(Serialize)]
struct ContigCounts {
    name: String,
    length: u64,
    counts: BaseCounts,
}

fn parse_sha256(value: &str) -> Result<String, String> {
    if value.len() == 64 && value.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        Ok(value.to_ascii_lowercase())
    } else {
        Err("expected 64 hexadecimal digits".into())
    }
}

pub fn run(args: &ReferenceArgs) -> Result<()> {
    let mut reference = TwoBitReference::open(&args.file)?;
    let uses_atlas = matches!(&args.command, ReferenceCommand::AtlasComposition { .. });
    if uses_atlas {
        ensure!(
            !matches!(args.assembly, Some(AssemblyLabel::Grch37)),
            "the primary RefSeq atlas uses GRCh38 coordinates"
        );
    }
    let sha256 = if uses_atlas || args.expect_sha256.is_some() {
        let actual = reference.sha256()?;
        if let Some(expected) = args.expect_sha256.as_deref() {
            ensure!(
                actual == expected,
                "reference SHA-256 mismatch: expected {expected}, found {actual}"
            );
        }
        if uses_atlas && args.assembly.is_none() {
            ensure!(
                actual == FROZEN_HG38_SHA256,
                "declare --assembly GRCh38 for a compatible reference with a different digest"
            );
        }
        Some(actual)
    } else {
        None
    };
    let source = ReferenceFile {
        format: "UCSC twoBit",
        format_version: 0,
        file_bytes: reference.file_bytes(),
        byte_order: reference.byte_order(),
        declared_assembly: args.assembly,
        sha256,
    };
    let result = match &args.command {
        ReferenceCommand::Catalog => ReferenceResult::Catalog {
            contigs: reference.contigs().collect(),
        },
        ReferenceCommand::Fetch {
            contig,
            start0,
            end0,
            strand,
            masked,
        } => ReferenceResult::Fetch {
            contig,
            start0: *start0,
            end0: *end0,
            strand: *strand,
            masked: *masked,
            sequence: reference.fetch(contig, *start0, *end0, (*strand).into(), *masked)?,
        },
        ReferenceCommand::Summary => {
            let contigs = reference
                .contigs()
                .map(|contig| (contig.name.to_owned(), contig.length))
                .collect::<Vec<_>>();
            let mut results = Vec::with_capacity(contigs.len());
            let mut total = BaseCounts::default();
            for (name, length) in contigs {
                let counts = reference.base_counts(&name)?;
                total += counts;
                results.push(ContigCounts {
                    name,
                    length,
                    counts,
                });
            }
            ReferenceResult::Summary {
                contigs: results,
                total,
            }
        }
        ReferenceCommand::AtlasComposition { atlas } => {
            let atlas = Atlas::open(atlas)?;
            ReferenceResult::AtlasComposition(Box::new(atlas.sequence_composition(&mut reference)?))
        }
    };
    let report = Report {
        schema_version: 1,
        reference: source,
        result,
    };
    let mut output = BufWriter::new(io::stdout().lock());
    serde_json::to_writer(&mut output, &report)?;
    output.write_all(b"\n")?;
    output.flush()?;
    Ok(())
}
