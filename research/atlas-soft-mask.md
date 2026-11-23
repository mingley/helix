# Does the nuclear GC ordering persist within the hg38 soft mask?

The earlier primary-chromosome scan found that, on each of the 24 nuclear
chromosomes, GC among called bases decreases from coding sequence to other
exons, to transcript spans outside exons, to sequence outside curated spans.
Each class can contain a different fraction of soft-masked bases. This study
asks whether the same ordering holds within each mask stratum and when all
classes receive the same masked fraction.
[Original counts and annotation definitions](genome-atlas.md#reference-sequence-composition-by-annotation)

## Results

The full four-class GC ordering holds on 24 of 24 nuclear chromosomes when
each class uses its own observed mixture. It holds on 20 of 24 among unmasked
bases, 23 of 24 among masked bases, and 23 of 24 after giving all classes
their chromosome's observed masked fraction. Every comparison had called
bases in both strata. The original computations counted 93,455 curated
transcript models on 25 primary chromosomes.
[Exact comparisons and source identities](atlas-composition/soft-mask/gc-ordering.json)

| Bases compared | Chromosomes in strict order | Exceptions |
| --- | ---: | --- |
| All called bases, observed separately in each class | 24/24 | None |
| Unmasked called bases | 20/24 | On chr2, chr10, chr15 and chr22, sequence outside transcript spans has higher GC than the portions of transcript spans outside exons. |
| Masked called bases | 23/24 | On chr19, transcript spans outside exons have higher GC than exons outside CDS. |
| Classes standardized to each chromosome's masked fraction | 23/24 | The same two classes reverse on chr19. |

On chr19, exons outside CDS exceed transcript spans outside exons by
1.242 percentage points in the observed mixture. Standardizing both to the
chromosome's 59.715% masked share reverses the difference to −0.112 percentage
points. The fractions of called bases that are masked differ between the two
classes, 33.87% and 56.89%, respectively.
[Per-chromosome nucleotide counts](atlas-composition/soft-mask/summary.json)

Across all 25 primaries, the annotation classes also differ substantially in
how much called sequence lies inside the mask.

| Annotation class | Masked fraction of called bases | GC among unmasked calls | GC among masked calls |
| --- | ---: | ---: | ---: |
| Coding sequence | 1.83% | 52.15% | 59.96% |
| Exon outside coding sequence | 22.75% | 45.31% | 45.20% |
| Transcript span outside exons | 47.53% | 40.23% | 42.23% |
| Outside curated transcript spans | 58.42% | 39.32% | 40.76% |

### A common range across all 24 nuclear chromosomes

With one masked fraction shared by every chromosome and class, all 72
adjacent GC comparisons hold strictly inside an open interval. Its endpoints
are approximately **40.0726% and 54.5171% masked**. At the lower endpoint nonexonic positions inside
chr15 transcript spans have the same GC fraction as positions outside those
spans. At the upper endpoint exons outside CDS and transcript
spans outside exons tie on chr19. Both bounds are reported as reduced
rational numbers in the data.

Five of the 72 comparisons change sign as the common fraction moves from
zero to one. The other 67 favor the original ordering for that entire range.

| Chromosome | Higher class in the original comparison | Lower class | Masked fraction at tie | Original order holds |
| --- | --- | --- | ---: | --- |
| chr2 | Transcript span outside exons | Outside spans | 32.2734% | Above tie |
| chr10 | Transcript span outside exons | Outside spans | 34.0301% | Above tie |
| chr15 | Transcript span outside exons | Outside spans | 40.0726% | Above tie |
| chr19 | Exon outside CDS | Transcript span outside exons | 54.5171% | Below tie |
| chr22 | Transcript span outside exons | Outside spans | 7.1398% | Above tie |

The table percentages are rounded. The independent checker evaluates all
four order definitions, each exact tie, both endpoints of the unit interval
and one interior point between consecutive boundaries. Each GC contrast is
linear between boundaries, so these checks exhaust its possible sign changes.
[Exact rational bounds](atlas-composition/soft-mask/gc-ordering.json) ·
[Independent calculation](atlas-composition/verify_soft_mask_ordering.py)

## Inputs and comparisons

These comparisons were specified before the nuclear mask-stratified scan.
We use the same frozen UCSC hg38 twoBit, RefSeq atlas and four disjoint
annotation classes as the original scan. The twoBit mask blocks define the
masked stratum. An unresolved N can belong to either stratum. Both Ns and
called bases contribute to coverage. GC denominators contain only A, C, G
and T. The Rust report includes chrM alongside the 24 nuclear chromosomes.
The ordering analysis uses chr1–22, chrX and chrY.

For each nuclear chromosome, compare the same three adjacent class pairs
under four definitions of GC fraction.

1. Combine masked and unmasked counts within each class to recover the earlier result.
2. Use only unmasked called bases in each class.
3. Use only masked called bases in each class.
4. Give each class the chromosome's overall masked fraction among all called
   bases, computed across its four disjoint classes. Its standardized GC
   fraction is `(1 − w) × GC_unmasked + w × GC_masked`.

For each definition, count complete strict four-class orderings and identify
every reversal, tie or undefined contrast. A stratum with zero called bases
has an undefined GC fraction. At a weight of exactly zero or one, only the
stratum with positive weight is needed.

A second comparison varies one shared masked fraction `w` from zero to one
across all four classes and all 24 nuclear chromosomes. For each adjacent
pair the GC difference at weight `w` is the linear interpolation of the two
stratum-specific differences. Intersect the 72 strict inequalities and
report the exact interval of weights and the contrasts that constrain its
endpoints. A tie excludes that weight. Arithmetic uses whole base counts
and exact rational numbers.

## Source checks

An independent case-preserving Python reader joins the original raw genePred
intervals to the reference. Compare all 1,000 chromosome, class, stratum and
base counts, then sum the two strata and reconcile the original 500 counts.
The original UCSC sequence-size report supplies a separate check on called
lowercase, called uppercase and unresolved bases across all 455 contigs.
The retained UCSC chrM API response and GenBank feature record provide a
base-by-base mitochondrial control.
[UCSC assembly build record](https://github.com/ucscGenomeBrowser/kent/blob/fa1ec8f7a7071301796d3050f7bb566679b55bbe/src/hg/makeDb/doc/hg38/hg38.txt#L335-L380)

The independent replay reproduced all 1,000 chromosome-by-class nucleotide
counts split by mask status. Summing the two strata recovered every count
in the original atlas. On the 25 primary chromosomes the image contains
1,531,174,404 masked called bases and 2,413 masked Ns. The single exonic
N previously traced to chr21 is among those masked Ns. Across all 455
contigs, the Rust summary reproduces the original UCSC build totals of
1,588,630,985 masked calls, 1,460,684,798 unmasked calls and 159,970,322 Ns.
One hosted Rust run took 1.44 seconds with peak memory below 96 MiB for
the primary atlas, including hashing, loading the atlas and counting.
[Independent raw-source replay](atlas-composition/soft-mask/source-check.json) ·
[Full-image comparison](atlas-composition/soft-mask/ucsc-control.json) ·
[Timed run](https://github.com/mingley/helix/actions/runs/45262329700/job/138000081478)

The mixed-case UCSC chrM sequence has five masked runs covering 418 bases.
Exactly one masked base lies in a coding sequence, the first A of ND1 at
`[3306,3307)`. GenBank accounts for 370 other masked bases as RNA (214 tRNA
and 156 rRNA) and 45 in its D-loop feature. The remaining two are between TRNL1 and ND1,
at `[3304,3306)`. Three positions are in both TRNI and TRNQ, so counting
unique tRNA coordinates matters here. The mitochondrial N is uppercase at
`[3106,3107)`. Its GenBank feature records that it preserves historical
genome annotation numbering.
[Mitochondrial feature intersections and source hashes](atlas-composition/soft-mask/mitochondrial-control.json)

UCSC built the mask from RepeatMasker and Tandem Repeats Finder. A change in
the order across these strata describes composition in this particular
assembly and annotation. It cannot identify a repeat family or establish a
causal effect of repeats on GC.
[UCSC hg38 download notes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/README.txt)

## Replay

Compute the exact comparisons directly from the retained integer counts.

```bash
counts=research/atlas-composition/soft-mask/summary.json
cargo run --release --locked -p genome-atlas --example soft_mask_ordering -- "$counts" > ordering.json
python3 research/atlas-composition/verify_soft_mask_ordering.py "$counts" ordering.json
```

To regenerate the counts from hg38, run the Rust and raw-source checks from
the repository root.

```bash
source scripts/env.sh
bash scripts/setup_reference.sh
reference="$HELIX_TOOLS_DIR/reference/hg38.2bit"
cargo run --release --locked -- reference --file "$reference" --assembly GRCh38 \
  --expect-sha256 1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0 \
  atlas-composition --by-soft-mask > atlas-soft-mask.json
python3 research/atlas-composition/verify_soft_mask.py "$reference" atlas-soft-mask.json
python3 research/atlas-composition/mitochondrial_soft_mask.py --report atlas-soft-mask.json
```
