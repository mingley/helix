# Independent variant-transfer source

Siraj et al., *Functional dissection of complex trait variants at
single-nucleotide resolution* (Nature, 2026), DOI
[10.1038/s41586-026-10121-6](https://doi.org/10.1038/s41586-026-10121-6),
provides a separate experiment with directly documented 200-base reference and
alternative reporter inserts and measurements in K562 and HepG2. Source selection
uses published design, units and availability; no S4 numeric outcome cell has
been decoded at this source-preparation stage.

The exact publisher-linked Supplementary Table 2 is **40,795,798 bytes**, SHA256
`4b6e20d1789021ee8f11f57af366a65680d25496de9351a9e781e9424716d824`.
Its single sheet, `S2.Variant Oligo Seqs`, has a title in row 1 and exact row-2
headers `Variant`, `Allele 1 Oligo`, `Allele 2 Oligo`. The raw workbook is cached
outside Git at `/workspace/.tools/variant-transfer-v1/S2.xlsx`.

The DNA-only inventory contains 365,762 unique IDs with no duplicate or missing
ID/construct cells. Exactly 336,475 nominal primary-chromosome SNV pairs have
two 200-base inserts, one difference at zero-based index 99, and bases matching
the ID's Ref and Alt alleles. Four of these contain ambiguous bases elsewhere;
336,471 pairs meet the A/C/G/T-only contract. Other lengths, indels and ambiguous
DNA remain in the inventory rather than being silently repaired.

The file named `siraj2026-S2-header-before-DNA.json` records a header inspection
and **one first DNA-row preview**. It predates the full DNA inventory, rather
than claiming that no DNA row was previewed. No reporter outcome is in S2.

## Assay and endpoint contract

The Methods describe 200-base genomic inserts upstream of a minimal TATA
promoter and GFP, with a 20-base barcode in GFP's 3′ UTR. Main complex-trait
libraries were measured in K562, SK-N-SH, HepG2 and A549; eQTL libraries in K562,
SK-N-SH, HepG2 and HCT116; saturation mutagenesis in K562 and HepG2. Each library
had at least five transfection replicates. These are separate cultured-cell
experiments, not donor-level replication or endogenous gene-expression effects.

Main-panel S4 allelic effects are the DESeq2 **Alt-versus-Ref by RNA-versus-DNA
interaction in delta log2 fold-change units**, following source summit-shift
normalization. They are not the same estimator as the earlier Agarwal published
mean log2(RNA/DNA) activity. The source excludes either-allele DNA counts below
30 or zero RNA counts and selects the library with highest documented plasmid
support when a variant occurs in several libraries. A transfer test must retain
this assay/estimator shift, source processing and missingness explicitly. It
must not select only source-significant or source-active variants.

The exact official S4 workbook is the publisher's `MOESM6_ESM.xlsx`,
224,405,074 bytes, SHA256
`723c4a0e1cfdb2c9c6ae55b3fa9206855dd95e0241319b9fc8d0831d3b42f9c2`.
It was acquired only after the committed
[semantic protocol](../../../experiments/variant-transfer-v1.json). Exact S4
headers and source metadata require a separately committed binding before numeric outcome parsing. No
satmut, shrunken S7, meta-analysis or diplotype endpoint is an automatic fallback.

The single sheet is `S4.MPRA results`; row 2 columns A, AN and AO are exactly
`Variant`, `log2Skew in HEPG2` and `log2Skew in K562`. Columns AX and AY are
`log2Skew.SE_HEPG2` and `log2Skew SE in K562`. ID-column-only extraction finds
304,278 rows and unique IDs, with no duplicates or blank IDs. The archive's ten
members pass bounded path/size/CRC/SHA checks; total uncompressed size is
1,078,863,882 bytes. All numeric effect columns remain uninspected at this stage.

Neither the S4 legend nor the available final-workbook metadata defines an
accepted missing token. Intermediate R exporters do not establish final Excel
encoding. The [missing-token review](missing-token-evidence.json) therefore
specifies **no accepted text or blank token**. A required NA, blank, formula,
error, boolean, numeric string or other nonordinary numeric cell must stop the
source gate, rather than being silently excluded or imputed.

The separately authorized official S30 design sheet resolves technical-control
identity. Its bundle is 19,928,224 bytes, SHA256
`3e62a3876b70e9c75951df26ed09ba552c3ba6f9f3b318934c1b01a5a42403b1`.
Only `S30.MPRA Controls` A:E fields `ID`, `Control Type`, `Source`, `Libraries`
and `Sequence` were interpreted; all other bundled scientific tables remain
uninterpreted. S30 contains 2,835 rows (99 positive, 104 allelic, 2,632 negative
controls), including 2,830 full 200-base sequences and five shorter sequences.
Either-allele full 200-base exact/RC matching flags 41 reference-valid external
variant IDs as reused technical controls. Their full components must be excluded
after all required unions. Source-ID suffixes/prefixes are not guessed or stripped
to manufacture an identity match. The complete metadata/duplicate/control ledger
is retained, including repeated S30 control IDs and the five shorter inserts.

Saturation mutagenesis is specifically unsuitable as an unexamined substitute:
Methods prose describes baseline-minus-substitution posterior activity, while
the pinned public `satmut.R` defines substitution-minus-baseline. Its empirical
Bayes processing also differs from the main-panel estimator. This disagreement
is retained; outcome agreement must never decide its sign.

## Assembly validation and exposure guards

The Methods and public analysis code operate on hg19 variant coordinates and
explicitly lift external hg38 annotations into hg19. New source coordinates
are never compared numerically with the old hg38 panel. The strict mapper uses
the official UCSC `hg19ToHg38.over.chain.gz`, SHA256
`5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d`.
For a centered SNV at one-based position p, the candidate source insert is
the zero-based half-open interval `[p-100,p+100)`.

Every accepted 200-base insert must map wholly within a single aligned chain
block, have one unique coordinate/strand mapping before DNA validation, and
match the entire frozen hg38 reference sequence exactly. Reverse chain mappings
use the corresponding reverse complement. Eight synthetic controls check
forward/reverse coordinates, whole-interval rejection, overlapping-block search,
inclusive guard boundaries and reverse-complement canonicalization. This is a
bounded assembly-consistency check, not proof of genome-wide sequence uniqueness.

Among 336,471 A/C/G/T-only SNV pairs, 332,719 match hg38 exactly after unique
whole-insert mapping. Exclusions retain 1,040 inserts with no single-block
mapping, 2,700 full-insert reference mismatches and 12 nonprimary mapped contigs;
none had an ambiguous whole-insert mapping. No assembly mismatch is repaired.

`design-mapping-v1` is the initial rowwise inventory. Its prefiltered
`eligible-pairs` file **must not define the final protocol cohort**.
`design-mapping-v2/all-reference-validated-pairs.jsonl.gz` retains all 332,719
validated pairs **before** prior-locus/exact-match exclusion. The final protocol
must form all required transitive 1-Mb start-bin, actual interval-overlap,
parent, exact/RC and 90%-identity
components on that full population, then exclude whole components touched by
an explored locus or homologous sequence. Adding homology after rowwise
prefiltering would lose required exclusion edges.

The upper bound using the frozen coordinate rule contains 2,828 components, of which
2,599 contain an old-locus 10-kb hit. At most 4,373 variants in 229 components
survive before the remaining relations and source-class/coverage gates.
This is not a completed homology check or final admissible cohort. Mandatory
additional unions/exposure flags can merge or remove those survivors; they
cannot restore an excluded component. Prior HBB/Kircher and Agarwal outcomes
remain explored evidence.

The earlier `coordinate-upper-bound.json` uses touched-bin unions and is a
superseded preparation calculation, not a frozen-v1 bound. The corrected
`coordinate-upper-bound-frozen-rule.json` follows `floor(start/1Mb)` plus
positive interval overlaps. Five controls distinguish crossed-bin boundaries,
positive overlaps, adjacency, within-bin grouping and transitive overlap.
Both preparations happen to leave the same 4,373/229 upper bound; the different
component counts are retained rather than relabeled.

## Distribution and licenses

The publisher's article states **CC BY 4.0**, with the usual exception for
separately credited material. S2 contains only its design sheet and standard
OpenXML properties; no contrary license or credit was found in its package.
Publisher terms and source attribution are retained in the methods metadata.
The separate Zenodo record **15297965 is CC BY-NC 4.0**. Its license is not
silently replaced with the article's license, and no Zenodo response ZIP was
downloaded for this preparation.

The public analysis repository `julirsch/finemapped_mpra` is pinned at commit
`1e38d1e3d4414dfa3a62b22110eba0136eda197d`. Its **GPL-3.0** license accompanies
the five retained README/license/plain-R source files under `producer-code/`.
They are inspected provenance, not imported into the laboratory's implementation.
No notebooks, rendered figures, assay payloads or numeric fixtures were acquired
from that repository. Its README describes a preprint-era code release, so the
commit is not authenticated proof of the final producer pipeline.

Alternative metadata candidates were considered without acquiring new response
payloads: Tewhey 2016 offers 150-base paired constructs and HepG2 but has a 2018
correction whose scope remains unresolved; Ulirsch 2016 offers K562 but uses a
quantile-normalized endpoint and lacked an established reuse license here.
Neither is substituted according to observed performance. Broad variant-effect
prediction and saturation analysis are established prior art; this source alone
does not establish a novel biological mechanism.

The mapping summaries describe that processor’s DNA-only I/O. The separate S4
acquisition was already underway during the later mapping replay, after the
semantic freeze; no S4 numeric effect was interpreted by either mapping run.
The [source semantic review](source-semantics-review.json) records the exact
units, sign, published QC/library-selection basis and remaining gates.
