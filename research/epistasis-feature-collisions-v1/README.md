# Feature-indistinguishable constructs in the S18 design

This is an **outcome-free analysis of published design DNA**, not evidence of a
measured interaction or of a new biological mechanism. Among the 7,542 previously
audited S18 four-construct windows, **282 windows in 78 variant-pair identifiers**
contain distinct constructs with exactly the same 32 canonical 3-mer counts and
GC count. Those windows contain **506 within-window equal-feature construct-pair
occurrences**, corresponding to **480 distinct unordered exact DNA-sequence
pairs** in their supplied orientation. None of the 506 pairs are whole-sequence
reverse complements. These are design constructs, not independent experiments;
windows, variant identifiers, and sequence pairs are different counting units.

The existing activity model divides those integer counts by fixed denominators
198 and 200 for 200-base sequences. Therefore equal counts give identical model
input features. **Any deterministic predictor that uses only these features and
the same cell-context input must make identical predictions** for each such
pair, even if it is nonlinear. A model supplied with other sequence features,
position, construct identity or other information need not be indistinguishable
on these constructs. This result is stronger than the previously established
zero-interaction identity for the particular *linear* 3-mer-plus-GC model, but
it is a specific limitation of the representation, not evidence that the actual
constructs have the same or different activity. K-mer-spectrum collisions and
their information loss are established concepts; this audit quantifies them on
this selected source and does not claim priority.

## Exact square classification

Let `F` be the 33-integer feature vector and write
`v = F(RR)`, `a = F(AR) − F(RR)`, `b = F(RA) − F(RR)`. All accepted source
edits are at least three bases apart. No 3-mer contains both edit sites, so
the exact four vectors are `v`, `v+a`, `v+b`, and `v+a+b`: GC counts also add
site by site. All vector operations are integer and require no fitted values.
The code directly recomputes all four complete vectors and checks this identity
for every audited window.

This structure exhausts the possible within-window feature collisions:

| Exact delta relation | Equal-feature construct pairs | Windows | Pair identifiers |
| --- | --- | ---: | ---: |
| `a = 0`, `b ≠ 0` | RR–AR and RA–AA | 104 | 32 |
| `a ≠ 0`, `b = 0` | RR–RA and AR–AA | 120 | 32 |
| `a = b ≠ 0` | AR–RA | 46 | 10 |
| `a = −b ≠ 0` | RR–AA | 12 | 4 |
| `a = b = 0` | all six distinct pairs | 0 | 0 |
| none of these | none | 7,260 | 2,078 |

The pair-identifier counts describe identifiers represented in each row; they
are not a general rule that different classifications cannot share a pair ID.
All 282 windows with a collision had passed the separate stringent hg38
reference rule; all 7,542 reference inserts matched hg19. The hg38 rule is
reported as context, not used to select a more favorable outcome subset.

For a concrete nonzero-delta example, the S18 `var1/left` window for
`chr17:27410544:G:A` and `chr17:27410548:G:A` has equal AR and RA feature
vectors, even though each differs from RR. For another, the `var1/left` window
for `chr19:1456544:A:G` and `chr19:1456653:G:T` has equal RR and AA feature
vectors although neither single-edit delta is zero. These coordinates and
alleles are the source's hg19 identifiers. The exact 200-base sequences and
source row numbers are bound through the S18 workbook and per-window ledger.

For a possible future test, suppose a matched equal-feature pair has observed
activities `y1` and `y2` in one cell, on one verified scale, and both are given
equal evaluation weight. Every deterministic feature-only predictor gives the
same value `t` to both, so its observed pairwise mean squared error satisfies

```text
((t - y1)^2 + (t - y2)^2) / 2
  = (t - (y1 + y2)/2)^2 + (y1 - y2)^2/4
  ≥ (y1 - y2)^2/4.
```

This is a finite-sample bound for those observed labels, including their
measurement error; it is not automatically an estimate of irreducible
biological variation or held-out performance. Pair occurrences here can share
constructs or sequences, so their individual bounds cannot simply be summed as
if all observations were distinct and independently sampled. No values of `y1`
or `y2` were inspected for this analysis.

## Inputs and reproducibility

The program checks SHA-256 bindings for the original S18-only design workbook,
the prior square-design audit, the two reference-validation ledgers, and the
prior S18 row reader. It checks each source square and its identity before
computing features. It retains all previously eligible windows; it reads no S20
workbook, reporter measurements, significance calls, or model coefficients.
The separately published source and related outcomes have already been explored
in this repository, so this is not a claim of an entirely unexposed experiment.

[`results/ledger.jsonl`](results/ledger.jsonl) records each source identity,
reference status, both complete integer deltas, classification and all matching
construct pairs, with a whole-sequence reverse-complement flag.
[`results/summary.json`](results/summary.json) gives counts and feature order;
[`results/manifest.json`](results/manifest.json) binds the implementation and
both outputs. The [independent check](../epistasis-feature-collisions-independent/README.md)
uses a separate source reader and feature calculation.

From the repository root, with the research environment activated, choose a new
output path:

```sh
source scripts/research-env.sh
python -m unittest discover -s research/epistasis-feature-collisions-v1 -p 'test_*.py'
python research/epistasis-feature-collisions-v1/audit.py --out /tmp/s18-feature-collision-replay
cmp research/epistasis-feature-collisions-v1/results/ledger.jsonl /tmp/s18-feature-collision-replay/ledger.jsonl
cmp research/epistasis-feature-collisions-v1/results/summary.json /tmp/s18-feature-collision-replay/summary.json
cmp research/epistasis-feature-collisions-v1/results/manifest.json /tmp/s18-feature-collision-replay/manifest.json
```

The script refuses an existing output directory. Synthetic controls check a
zero-feature edit on either or both sites, equal and opposite nonzero deltas,
absence of a collision, reverse-complement invariance, and a nearby interaction
that violates the separated-edit identity. They are not outcome validation.

Siraj's pair designs were selected using trait and genomic annotations, not
sampled uniformly from human genetic variation. The counts are properties of
these designs and cannot be interpreted as population frequencies. The
availability, matched normalization, replication and uncertainty of the
corresponding reporter measurements remain separate source gates. A future test
of these forced equalities needs a fixed cohort and appropriate measurement
model; a statistical discrepancy would be evidence about the stated reporter
assay and features, not on its own endogenous regulation or causality.
