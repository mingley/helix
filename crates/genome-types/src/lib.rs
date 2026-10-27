//! Assembly-indexed, zero-based, half-open genomic contracts.
//!
//! ```compile_fail
//! use genome_types::{Variant, GRCh37, GRCh38};
//! fn accepts_grch38(_: Variant<GRCh38>) {}
//! fn cannot_mix_assemblies(v: Variant<GRCh37>) { accepts_grch38(v); }
//! ```
use anyhow::{bail, ensure, Result};
use serde::{Deserialize, Serialize};
use std::{fmt, marker::PhantomData};

mod sealed {
    pub trait Sealed {}
}

pub trait Assembly: sealed::Sealed + Clone + fmt::Debug + PartialEq + Eq {
    const NAME: &'static str;
    /// Primary chromosomes 1..22, X, Y, MT; assembly-specific reference bounds.
    const CHROMOSOME_LENGTHS: [u64; 25];
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GRCh38;
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GRCh37;
impl sealed::Sealed for GRCh38 {}
impl sealed::Sealed for GRCh37 {}
impl Assembly for GRCh38 {
    const NAME: &'static str = "GRCh38";
    const CHROMOSOME_LENGTHS: [u64; 25] = [
        248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
        138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
        83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569,
    ];
}
impl Assembly for GRCh37 {
    const NAME: &'static str = "GRCh37";
    const CHROMOSOME_LENGTHS: [u64; 25] = [
        249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
        141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753,
        81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16569,
    ];
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(bound = "")]
pub struct GenomicPosition<A: Assembly> {
    chromosome: String,
    offset: u64,
    #[serde(skip)]
    assembly: PhantomData<A>,
}
impl<A: Assembly> GenomicPosition<A> {
    pub fn new(chromosome: &str, offset: u64) -> Result<Self> {
        let chromosome = chromosome.strip_prefix("chr").unwrap_or(chromosome);
        let (chromosome, index) = match chromosome {
            "X" => ("X".to_owned(), 22),
            "Y" => ("Y".to_owned(), 23),
            "M" | "MT" => ("MT".to_owned(), 24),
            number => match number.parse::<usize>() {
                Ok(n @ 1..=22) => (n.to_string(), n - 1),
                _ => bail!("unsupported chromosome: {chromosome}"),
            },
        };
        ensure!(
            offset < A::CHROMOSOME_LENGTHS[index],
            "coordinate outside {} chromosome {chromosome}",
            A::NAME
        );
        Ok(Self {
            chromosome,
            offset,
            assembly: PhantomData,
        })
    }
    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }
    pub fn offset(&self) -> u64 {
        self.offset
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(bound = "")]
pub struct Variant<A: Assembly> {
    position: GenomicPosition<A>,
    reference: char,
    alternate: char,
}
impl<A: Assembly> Variant<A> {
    pub fn snv(position: GenomicPosition<A>, reference: &str, alternate: &str) -> Result<Self> {
        fn base(s: &str) -> Result<char> {
            let s = s.to_ascii_uppercase();
            match s.as_str() {
                "A" | "C" | "G" | "T" => Ok(s.chars().next().unwrap()),
                _ => bail!("only unambiguous biallelic SNVs are supported: {s}"),
            }
        }
        let reference = base(reference)?;
        let alternate = base(alternate)?;
        ensure!(
            reference != alternate,
            "reference and alternate must differ"
        );
        Ok(Self {
            position,
            reference,
            alternate,
        })
    }
    pub fn position(&self) -> &GenomicPosition<A> {
        &self.position
    }
    pub fn reference(&self) -> char {
        self.reference
    }
    pub fn alternate(&self) -> char {
        self.alternate
    }
    /// A local display key, explicitly NOT a GA4GH computed identifier.
    pub fn key(&self) -> String {
        format!(
            "{}:{}:{}:{}>{}",
            A::NAME,
            self.position.chromosome,
            self.position.offset,
            self.reference,
            self.alternate
        )
    }
    pub fn interval(&self) -> (u64, u64) {
        (self.position.offset, self.position.offset + 1)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Consequence {
    Missense,
    Synonymous,
    Splice,
    Stop,
    Intron,
    Utr,
    Other,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ClinicalLabel {
    Benign,
    Pathogenic,
}
impl ClinicalLabel {
    pub fn positive(self) -> bool {
        self == Self::Pathogenic
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BiologicalContext {
    pub gene: String,
    pub population: String,
    pub endpoint: String,
}
impl BiologicalContext {
    pub fn hbb_germline() -> Self {
        Self {
            gene: "HBB".into(),
            population: "human germline; ClinVar aggregate assertions".into(),
            endpoint: "binary clinical classification, not a measured molecular effect".into(),
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct FeatureRecord {
    pub id: String,
    pub variant: Variant<GRCh38>,
    pub consequence: Consequence,
    pub context: BiologicalContext,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Probability {
    numerator: u32,
    denominator: u32,
}
impl Probability {
    pub fn new(numerator: u32, denominator: u32) -> Result<Self> {
        ensure!(
            denominator > 0 && numerator <= denominator,
            "invalid probability {numerator}/{denominator}"
        );
        Ok(Self {
            numerator,
            denominator,
        })
    }
    pub fn value(self) -> f64 {
        self.numerator as f64 / self.denominator as f64
    }
    pub fn numerator(self) -> u32 {
        self.numerator
    }
    pub fn denominator(self) -> u32 {
        self.denominator
    }
    pub fn label(self) -> ClinicalLabel {
        if u64::from(self.numerator) * 2 >= u64::from(self.denominator) {
            ClinicalLabel::Pathogenic
        } else {
            ClinicalLabel::Benign
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn canonical_snv_and_inter_residue_coordinates() {
        let v = Variant::<GRCh38>::snv(GenomicPosition::new("chr11", 5227001).unwrap(), "t", "a")
            .unwrap();
        assert_eq!(v.key(), "GRCh38:11:5227001:T>A");
        assert_eq!(v.interval(), (5227001, 5227002));
    }
    #[test]
    fn rejects_invalid_representations() {
        for (r, a) in [("A", "A"), ("AT", "A"), ("N", "T"), ("G", "<DEL>")] {
            assert!(Variant::<GRCh38>::snv(GenomicPosition::new("11", 1).unwrap(), r, a).is_err());
        }
        assert!(GenomicPosition::<GRCh38>::new("23", 1).is_err());
        assert!(GenomicPosition::<GRCh38>::new("11", 135_086_622).is_err());
        assert!(GenomicPosition::<GRCh38>::new("11", 135_006_516).is_ok());
        assert!(GenomicPosition::<GRCh37>::new("11", 135_006_516).is_err());
        assert_eq!(
            GenomicPosition::<GRCh38>::new("chr01", 0)
                .unwrap()
                .chromosome(),
            "1"
        );
        assert!(Probability::new(1, 0).is_err());
        assert!(Probability::new(2, 1).is_err());
    }
}
