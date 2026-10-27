# Metadata-only context-study source review

The strongest current route is the independently produced Agarwal2025 panel being
reviewed by the other source agent, provided its exact DNA joins and data-reuse
scope are verified. These older studies are valuable, bounded backups with clear
matched-context designs, but their currently observed reuse terms are weaker.
No quantitative activity matrix, count table, or score file was requested here.

## Inoue2017: clean two-context design, permission/correction questions

- DOI `10.1101/gr.212092.116`, PMID27831498, PMC5204343, GEO **GSE83894**.
- GEO identifies 2,236 candidate liver enhancers plus204 controls assayed in HepG2
  with integration-competent versus integration-defective lentiviral components.
- The downloaded **column-description metadata only** explicitly labels
  `MT1/2/3` as nonintegrating and `WT1/2/3` as integrating RNA/DNA ratio replicates.
  The processed1.5MB matrix itself was not retrieved.
- The direct design file is identified in GSM2221119 metadata; raw sequence and
  exact `Region` joins still require validation after a frozen acquisition plan.
- Source genome build is GRCh37. Do not relabel it GRCh38.
- Article terms are CC-BY-NC4.0. Reviewed GEO metadata did not establish a separate
  CC0/CC-BY per-file data grant. NCBI imposes no restrictions itself on molecular
  data but explicitly acknowledges possible third-party intellectual-property
  claims; hosting alone is not a substitute for source terms.
- A2018corrigendum, DOI`10.1101/gr.237321.118`, exists. The retrieved official XML
  has bibliographic/frontmatter only. Its substantive correction is unresolved
  and must be checked before relying on the old processing interpretation.

## Klein2020: richer contexts, shared sequence ancestry

- DOI `10.1038/s41592-020-0965-y`, PMID33046894, PMC7727316, GEO **GSE142696**.
- GEO identifies2,440 designs in HepG2 across nine MPRA configurations, separate
  forward/reverse constructs, and nested192/354/678bp versions.
- Most named designs have three DNA/RNA replicate labels; length/orientation
  comparisons have four. Source Methods specifies technical-run pooling and two
  first-replicate exclusions during averaging. These are not extra independent
  biological sequence families.
- The study explicitly reuses the same2,440×171bp array from Inoue2017. Combining
  their records requires grouping by shared sequence ancestry, not study label.
- All processed means/replicates are small: the nine-design means file is122KB,
  its replicate file257KB, orientation means40KB, and length means48KB compressed.
  Exact FASTA/coordinate files are separately identified. None was requested.
- Methods describes log2(RNAcpm/DNAcpm) plus median-based replicate normalization.
  Resolve exact processed units/order from source code and file metadata before
  freezing transformations; these are not unnormalized molecular abundances.
- The PMC manuscript gives academic-research text/data-mining terms, not a general
  CC-BY data grant. MPRAflow's Apache2.0 software license does not license the
  processed experimental data.

`candidate-review.json` contains exact candidate endpoints, context/replicate
units, direct download URLs and remaining gates. Metadata snapshots and their
hashes are indexed in `manifest.json`. Published abstracts/series descriptions can
contain authors' aggregate findings; no per-sequence quantitative outcome table
was acquired or analyzed for selection.

Before either backup is used, freeze data terms, sequence identity/grouping,
control exclusions, endpoint transformations, replicate handling and split rules.
Neither study's design alone proves a complete endogenous regulatory mechanism.
