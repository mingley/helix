# Does the nuclear GC ordering persist within the hg38 soft mask?

The earlier primary-chromosome scan found that, on each of the 24 nuclear
chromosomes, GC among called bases decreases from coding sequence to other
exons, to transcript spans outside exons, to sequence outside curated spans.
Each class can contain a different fraction of soft-masked bases. This study
asks whether the same ordering holds within each mask stratum and when all
classes receive the same masked fraction.
[Original counts and annotation definitions](genome-atlas.md#reference-sequence-composition-by-annotation)

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

UCSC built the mask from RepeatMasker and Tandem Repeats Finder. A change in
the order across these strata describes composition in this particular
assembly and annotation. It cannot identify a repeat family or establish a
causal effect of repeats on GC.
[UCSC hg38 download notes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/README.txt)
