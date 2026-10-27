# Source normalization code trace

The strongest result is a concrete, externally specified operation to test:
the final-paper modeling table S12 may be a **per-cell sample z-score of the
source S11 ESS values, computed over exactly the S12 element subset**. This is
a hypothesis from downstream public code. We have not executed that code,
retrieved new assay tables, or validated the identity against source values
in this investigation. The original producer's activity-to-ESS construction
remains unidentified.

This directory records code and metadata acquired on 2026-10-27. Raw numerical
workbook acquisition and subsequent checks require the separately frozen
recovery protocol. No workbook, TSV, Parquet, checkpoint or archive payload
was downloaded here. Existing repository outcome knowledge was not used to
choose among new table values. Public source code and documentation can
contain numeric constants and previously reported assertions; those are
identified as external assertions, not our measurements.

## The operation specified by downstream code

Repository: [TaykhoomDalal/LentiMPRA-Enhancer-Processing](https://github.com/TaykhoomDalal/LentiMPRA-Enhancer-Processing).
Pinned head: `d2751a6ddff67ce4406288e9f6c3d932b0c12cdd`.
The checked source file is
[process.py](dalal/LentiMPRA/process.py), SHA-256
`a904d02c65072be4dd328d700958011c5c804e43799a7f4fc1cefd551559b6cb`
(32,951 bytes).

Its `validate_joint_ess` function, lines 633–661, specifies two distinct checks:

1. S11 `all_cell_types_ESS` contains three scores whose row sum is below
   `1e-9` in absolute value. It also checks multiples of `1/18000`.
2. For each cell, select S11 ESS values with identifiers present in S12,
   compute their mean and sample standard deviation (`ddof=1`), and compare
   `(S11_ESS - mean) / sample_sd` to the corresponding S12 score, with maximum
   absolute error at most `1e-12`.

The function is an assertion in a downstream verification program. Its
existence does not establish that it passes on the original files or that
the study authors used that sequence of operations. If subsequently
verified on official bytes, it would establish a numerical file relation.
It would still leave the historical creation of S11 ESS unresolved.

GitHub reports that the repository's initial commit,
`ce0779aa266f2f581b63007b5107010ddd74443d`, dates to
2026-08-21T19:37:31Z. The same function appears there at lines 637–665, with
an identical Python AST to the checked head. That initial code snapshot has
SHA-256 `5cf57ac415a3682d9cb4d89cd7ce9e81a67a5ae9f48bebf8b9a13ad3dc857eba`
(33,521 bytes). The reported head commit date is 2026-08-25T12:26:24Z.
These dates precede this workspace's 2026-10-27 investigation according to
public GitHub metadata. Commit dates are not proof of scientific independence,
an authenticated historical timestamp, or linkage to the original laboratory
pipeline. The head commit is reported unsigned.

The limitation is explicit in the downstream repository itself:
[common.py](dalal/scripts/common.py) says S11 ESS values are copied by exact
source identifier because their upstream normalization is unpublished.
Its [README](dalal/README.md) likewise says that normalization is unavailable.
We do not substitute its identifier aliases or missingness policies for our
frozen source selection and exact-identifier requirements.

## File identities predicted before acquisition

[common.py](dalal/scripts/common.py) supplies these expected identities. They
are downstream assertions awaiting verification against the official archive,
not publisher signatures or independently observed source-byte identities.

The official archive URL is
<https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-08430-9/MediaObjects/41586_2024_8430_MOESM5_ESM.zip>.
Its asserted size is 342,234,089 bytes and asserted SHA-256 is
`aefb1974e9ebe5c196a89a00f4e868c0abaf69e3d3952e79c74bfb42f217f6a4`.
All members below have prefix `2023-03-03628-s5/`.

| Member basename | Asserted bytes | Asserted SHA-256 |
| --- | ---: | --- |
| `2023-03-03628C-Table S10 - joint lib design.xlsx` | 8,301,566 | `f9ec4c6c2d230b0ef507a39e5d5c73b3b36f90ccc5758c33b442c7561e20e7e5` |
| `2023-03-03628C-Table S11 - jointdata.xlsx` | 53,557,805 | `9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e` |
| `2023-03-03628C-Table S12 - folds and performance.xlsx` | 3,564,742 | `0f9ff22d4b46fab36d419399e6ab743dc08a0813f812dc7bd96f90eb7e0f2b8a` |

The asserted S10/S12 identities equal the already committed preprint-associated
S7/S9 identities. That is a testable byte-identity prediction; filenames and
table numbers alone cannot establish it.

The downstream code identifies S11 activity sheets as
`HepG2_replicate_level`, `K562_replicate_level`, `WTC11_replicate_level`,
the corresponding `<cell>_summary_data` sheets, and `all_cell_types_summary`.
The combined activity headers are `name` and `<cell> [log2(rna/dna)]`.
`all_cell_types_ESS` uses `name` and `ElementSpecificityScore_<cell>`.
S12 uses `Fold assignment, joint library` with `Sequence ID` and
`<cell>, Element Specificity Score`. These are code-derived expected schemas,
not independently inspected workbook headers in this trace.

## What the producer-linked repositories establish

The article names MPRAflow 1.0 `count.nf` for per-replicate activity generation.
The public [shendurelab/MPRAflow](https://github.com/shendurelab/MPRAflow)
repository is therefore a producer-linked upstream pipeline, but the paper
does not pin a Git commit. Current public tags begin at v2.1; no v1.0 tag was
found in the acquired tag metadata.

We inspected a 2020 snapshot at
`6b136aacf0cfa5d46023668237a379d9d9edd05b`, whose own configuration says
version **2.2**, plus current head
`87d8b4c623274011841d3580f4b32a2bcd5bf54b`. Neither is certified as the
study's stated 1.0 implementation. In the historical
[merge_label.py](mpraflow-2020/src/merge_label.py), lines 133–154 show
pseudocount-adjusted, total-count-normalized RNA/DNA ratios and `log2` output
per element/replicate. `count.nf` defaults to ten observed barcodes.

A specific caution prevents treating generic pipeline code as proof of the
article's downstream normalization. In
[make_master_tables.R](mpraflow-2020/src/make_master_tables.R), lines 40–48
calculate median-normalized `nrm_reps`, but save `all_reps` instead. Lines
53–57 aggregate `all_reps` by condition alone, using `log2(mean(ratio))`.
The normalized object is not used in that saved average. This file is
byte-identical at the inspected current head. Its path history reports its
last change in 2019. These are observations about checked source code;
they do not show which postprocessing the study actually executed or
invalidate measurements made through an unobserved pipeline.

The paper-linked `visze/sequence_cnn_models` head
`7638fce137db0445123efe8d7e2c35e248fafc5f` supplies training and label-join
code, not an observed activity-to-ESS constructor. The checked
`regression_input.smk` joins supplied labels to sequence IDs and divides
supplied folds. Removing `_F` in that workflow is not a rule we import into
the frozen exact-ID analysis. The coauthor-associated `autosome-ru/MPRA-MNIST`
head `644590f095e33fa8ab484550f5395e8fc94e6105` similarly selects supplied
`*_Specificity_Score` or `*_log2` columns in its Agarwal loader; it does not
establish their original construction. We did not fetch its datasets or
notebooks, which may include numerical outputs.

## Additional final-table lead, with incomplete provenance

The downstream
[GAME evaluator](https://github.com/de-Boer-Lab/GAME-Agarwal-MPRA-joint-library-evaluator)
at `8f95a6ce026e0b1e09eba2b13c1421f35e16f583` names final S11 jointdata and
S10 design in its [data README](game/agarwal_joint_library/evaluator_data/README.md).
Its available [filtering script](game/agarwal_joint_library/evaluator_data/make_56k_filtered_file.py)
uses `dropna()` and exact-name intersections. It reads an unavailable local
intermediate `all_cell_type_measured.xlsx`; it neither constructs ESS nor
proves how that intermediate came from S11. The README also names a
preparation script absent from the checked tree and describes adapter trimming
that is not performed in the checked filtering script.

Its untouched candidate workbooks are recorded in `trace.json` with exact
Git blob IDs and sizes. A downstream copy is not independent experimental
validation or authenticated producer history. The generic 2026 MPRA QC
repository was also inspected as a negative lead; no linkage to the study's
ESS construction was established.

## Integrity and reuse

`retrieval.json` records URLs, retrieval times, response sizes, SHA-256 values
and failed searches. `manifest.json` inventories every other evidence file.
`trace.json` separates observed code from unverified source predictions.
Verify this snapshot offline from the repository root:

```bash
python3 research/context-normalization-trace/verify_evidence.py
```

Dalal processing code is MIT licensed; its retained `NOTICE.md` explicitly
keeps the underlying Agarwal article/data terms separate (CC BY 4.0 as stated
by that source). MPRAflow is Apache-2.0; the retained sequence-CNN code is MIT.
GitHub reports no repository license for GAME, MPRA-MNIST or the checked QC
fork; do not infer a new reuse grant or relabel their snapshots as Helix code.
Their public code is retained here as attributed research evidence. Original
copyright and license notices remain in the snapshots. None of these licenses
authenticates the historical normalization recipe.
