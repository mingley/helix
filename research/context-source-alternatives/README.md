# Metadata-first matched-cell-context source discovery

**Subsequent frozen acquisition:** protocol
`a089a6668ad73e943797df405b3ba510b1377190b416da1bd8df12bc368b77c2`
authorized downloading S7 and S9. Both complete-file MD5s passed. Their schemas
failed the primary endpoint gate: S7 contains design only, while S9 explicitly
contains **Element Specificity Scores**, not observed mean log2(RNA/DNA).
No model was fitted and no target was substituted. See
`datasets/context-panel/raw/semantic-gate.json` and its request receipts.
The metadata-only account below describes the earlier discovery phase.

No processed quantitative outcome table or model weight was acquired during this
discovery phase. This directory contains publication/method metadata, public
repository metadata and code references, and **ZIP directory metadata only**.
`candidates.json` records the current evidence and unresolved schema requirements.

## Strongest candidate: Agarwal et al.2025 joint lentiMPRA library

- *Massively parallel characterization of transcriptional regulatory elements*,
  [DOI10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9),
  PMID39814889, PMC11903340.
- A shared **60,000-cCRE joint library** was tested in HepG2, K562 and WTC11.
  This is distinct from the larger collection of different cell-specific
  libraries and is the relevant same-input context experiment.
- Official joint accessions: HepG2 **ENCSR405QCT**, K562 **ENCSR203UFY**,
  WTC11 **ENCSR336MKI**.
- WTC11 is a human iPSC line, RRID:CVCL_Y803, passages43–49.
- Three independent lentiviral infections per cell context used the same
  constructed joint library. They are within-study replicates, not separate
  studies or independent library construction.
- The methods define per-replicate log2(RNA/DNA), centering by that replicate's
  full-library median and averaging across three replicates. This is relative
  library-normalized activity, not absolute transcription per DNA copy.
- ESS further transforms values by cell-level normalization and element-level
  contrasts. It is a different endpoint and must not be substituted for activity.

The article and author Zenodo data record explicitly grant **CC BY4.0**.
This is a distinct experimental study from Kircher2019, with overlapping
investigators and technology. It is not independent laboratory replication.

The source common panel was selected using earlier activity bins and high
across-cell variance. A study on this fixed panel must retain that ascertainment;
it cannot estimate genome-wide activity or cell-specificity prevalence.
MOI also differs: HepG2=50, K562=10, WTC11=10. Integration, cell state and technical
context limit a purely causal interpretation of a cell-label contrast.

## Bounded acquisition options, pending a protocol freeze

The paper links [Zenodo10558183](https://zenodo.org/records/10558183).
Its metadata expose:

| Artifact | Declared bytes | Acquisition status |
| --- | ---: | --- |
| `Table S7 - joint lib design.xlsx` | 8,301,566 | Not acquired |
| `final_joint_dump.zip` | 3,967,659,598 | Whole archive not acquired |

The workbook title is not proof that it contains only design metadata. The
final article calls its joint/ESS table Supplementary Table10, whereas this
older deposit uses preprint numbering. It may contain quantitative outcomes;
do not fetch it before the acquisition/analysis freeze.

The large ZIP supports normal HTTPS byte-range requests. Only its 22-byte
end-of-central-directory record and its exact 39,943-byte central directory were
retrieved. No member payload was read.

The directory identifies **30 small TSV members**:

```text
final_joint_dump/scores/final/HepG2/joint_results_testaug/testfold_{1..10}.tsv
final_joint_dump/scores/final/K562/joint_results_testaug/testfold_{1..10}.tsv
final_joint_dump/scores/final/WTC11/joint_results_testaug/testfold_{1..10}.tsv
```

Together these occupy **2,296,453 compressed / 6,193,852 uncompressed bytes**.
Every member's path, compression method, source CRC32 and local-header offset is
retained in `joint-archive-index.json`. This allows bounded extraction after a
freeze without downloading gigabytes of model checkpoints.

However, their **observed-target schema remains unverified**. The pinned author
CNN prediction writer outputs ID and model-prediction columns; it joins observed
labels only in a subsequent workflow stage. Therefore these fold files must not
be treated as measured activity merely because the directory is named scores.
The workbook may be the actual observation/design source. Acquisition must fail
closed if no exact oligo-sequence join and measured-mean activity are available;
predictions and ESS are not acceptable fallback targets.

The archive's outer MD5 is recorded from Zenodo metadata but cannot be claimed
verified without downloading the whole archive. After a permitted partial
extraction, verify each member against its ZIP CRC and retain a new SHA-256 plus
the exact range/source receipt. The workbook has its own complete-file MD5.

## Access and alternative-study limits

The three official ENCODE metadata JSON requests returned 403 Forbidden here;
no access bypass or outcome download was attempted. The author-published Zenodo
record is an independent documented distribution route.

Tewhey2016 and Ulirsch2016 publication metadata were found, but a bounded
matched-input multi-context matrix, exact schema and clear data license were
not established in this pass. Full-text services returned500/429; they are not
silently promoted to ready alternatives. Tewhey2016 has a2018 correction that
also needs review. Broad Sherwood/GTEx name matches were insufficient for a
specific downloadable endpoint.

## Provenance and source code

`retrieval.json` records hashes, sizes, URLs, extraction descriptions and times.
Article XML was used only to extract methods, availability and licensing
metadata; no linked processed/supplementary quantitative table was retrieved.

Author code/README references are retained from:

- visze/sequence_cnn_models, source ref
  `7638fce137db0445123efe8d7e2c35e248fafc5f`.
- autosome-ru/human_legnet, source ref
  `a20dbc7d0b08daae3af1aebcfb8a71d9f198ab5d`.

Both repositories declare MIT; their license notices are retained. These files
are reference evidence, not instructions or a request to run pretrained models.
