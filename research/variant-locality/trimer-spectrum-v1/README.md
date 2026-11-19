# Exact counts for sequences sharing one trimer spectrum

The panel model uses 32 reverse-complement canonical 3-mer frequencies and GC
fraction. An exact count under stricter matching conditions gives a lower bound
on how many different sequences this representation groups together. Here we
fix both terminal dinucleotides and the complete vector of all 64 oriented
3-mer counts. The counts fix sequence length. Counting the first two bases
and the last base of every 3-mer then gives the total GC count, so every
sequence counted here has the same 33 model features.

Two prior Agarwal reporter designs were selected for the reference-coordinate
audit, one in each orientation. Their 200-base inserts are verified against
the retained reference responses and core SHA-256 values. The negative-strand
reference is reverse complemented before counting.

| Design ID | First and last two bases | Exact class size, rounded |
| --- | --- | ---: |
| WTC11_seq1000_F | GT and CT | 1.17225 × 10^76 |
| HepG2_DNasePeakNoPromoter12142 | AT and CT | 1.55336 × 10^77 |

The exact integers are stored as decimal strings. All counts include the
original insert. In both examples, reverse complementing a member changes its
fixed prefix, so every counted member belongs to a different reverse-complement
orbit. The model's full equivalence class contains both orientations of each
orbit. These are counts over all 200-base strings satisfying the constraints.
The two verified designs supply the reference spectra for the counting
calculation. [Exact counts](examples.json) ·
[Reference coordinate check](../../context-source-alternatives/design-coordinate-crosscheck.json)

## Counts across the frozen panel

The exact count ran on all 55,276 retained 200-base sequences. For 54,983
records (99.47%), the lower bound is at least 10^70 distinct
reverse-complement orbits per record. Every record has at least 2.6091 × 10^24
such orbits, with the displayed number rounded down. Terminal dinucleotides
alone rule out reverse-complement pairs within the counted class for 51,284
records.

| Panel statistic | Exact fixed-end class size | Lower bound on reverse-complement orbits |
| --- | ---: | ---: |
| Minimum | 2.60911 × 10^24 | 2.60911 × 10^24 |
| 25th percentile | 1.79706 × 10^77 | 1.71644 × 10^77 |
| Median | 8.60568 × 10^78 | 8.19540 × 10^78 |
| 75th percentile | 2.63805 × 10^80 | 2.50801 × 10^80 |
| Maximum | 2.94894 × 10^86 | 2.94894 × 10^86 |

These counts describe possible strings under the feature constraints for each
input. The panel contains only the 55,276 original designs. The exact summary
retains the complete integers and records each nearest-rank percentile.
[Panel summary](panel-summary.json) ·
[Recorded scan](https://github.com/mingley/helix/actions/runs/44610553987)

## Counting argument

Represent a 3-mer `abc` by a directed edge from `ab` to `bc`. Each observed
copy gives an edge of the same type. A linear string with the specified ends
traces an Euler trail from its first dinucleotide to its last, using each edge
once. Denote the final vertex by `t`, the edge multiplicities by `m_e`, and the
outdegrees, including loops, by `d_v`. Initially distinguish copies of each
edge with labels.

In a labelled trail, the last exits at vertices other than `t` form a tree
directed toward `t`. Along any chosen last exit to a non-root vertex, the
last exit at that next vertex must occur later, which rules out cycles.
Conversely, pick such a tree, reserve its edge as the last exit at each
non-root vertex, and order all other outgoing edges. Following the lists from
the fixed first vertex gives a unique labelled trail. The degree balances force
it to stop at `t`. Any remaining edges would imply a remaining non-root last
exit path into `t`, which contradicts stopping there. The same argument
covers a closed trail whose start equals `t`.

Let `T_t` be the number of labelled trees directed toward `t`. The directed
matrix tree theorem computes it from the out-Laplacian minor with the row and
column of `t` removed. Loops cancel in that minor. They still contribute to
degrees and outgoing orderings. Dividing the labelled count by
`product_e(m_e!)` removes the independent permutations of copies of each
edge type. The exact string count is

```text
N = T_t * product_v((d_v - 1 + [v = t])!) / product_e(m_e!)
```

There is no rotation factor because the strings are linear and the beginning
is fixed. The program uses integer Bareiss elimination for the determinant.
For an arbitrary panel record, the count gives at least `ceil(N/2)`
reverse-complement orbits, and exactly `N` among the counted strings when the
fixed ends rule out a reverse-complement pair.

The result extends a standard way to think about k-mer-preserving
rearrangements. Altschul and Erickson used Euler walks to preserve dinucleotide
or trinucleotide usage, and uShuffle implements general k-let preservation.
The numbers above apply the counting formula to Helix's feature contract.
[Altschul and Erickson, 1985](https://pubmed.ncbi.nlm.nih.gov/3870875/) ·
[Jiang et al., 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2375906/)

## Constructive controls

A two-base anchor lets us change longer sequence context while preserving
every oriented trimer. For any two-base string `D`, swapping `A` and `B` in
`P + D + A + D + B + D + S` preserves the spectrum. Appending each segment
and its trailing anchor starts and ends at the same graph vertex `D`. Any
candidate should be checked for exact identity and reverse complementation
as well as feature equality.

As a synthetic illustration, put 49 fixed two-base blocks after an initial
`TAC`, append `AC` after every block and finish with `T`. Use three copies of
each possible block and one extra `AA`. Each block contributes the same four
3-mers regardless of its position. The result has length 200, fixed GC count
98, and exactly `49! / (4! × (3!)^15)` distinct block permutations. Roughly
5.39 × 10^49 strings therefore share one feature vector. All begin and end
in `T`, so reverse complementation takes them outside the family. A shorter
anchor swap is checked in the tests and changes just two sites, three bases
apart. These constructions make model equality testable on sequence pairs.
Reporter measurements would be needed to judge the biological consequence.

## Replay

Run from the repository root with Node.js 18 or later.

```bash
node --test research/variant-locality/trimer-spectrum-v1/test_count.cjs
node research/variant-locality/trimer-spectrum-v1/audit.cjs --verify
node research/variant-locality/trimer-spectrum-v1/audit.cjs --panel | diff -u research/variant-locality/trimer-spectrum-v1/panel-summary.json -
```

The tests enumerate all 21,824 A/C/G/T strings of length three through seven
and compare the count against direct grouping by ends and trimer multiset.
The replay of the two examples checks coordinates, orientation and DNA hashes.
A separate subset-based determinant computation with prime factorial
valuations reproduces both example counts.

The optional panel run verifies the frozen 55,276-record file by SHA-256 and
computes an exact count per input record. It reports nearest-rank quantiles and
counts above declared orbit thresholds. The parser loads the panel JSON,
including previously explored specificity fields, and the calculation uses
only `sequence`. The distribution weights each retained panel record equally.
