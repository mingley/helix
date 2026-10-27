# LentiMPRA enhancer parquet processing

This repository converts the supplementary workbooks from Agarwal et al.,
*Massively parallel characterization of transcriptional regulatory elements*,
into four analysis-ready Parquet files.

It:

1. downloads and verifies the 342 MB supplementary archive;
2. extracts the six workbooks needed for the large-scale and joint libraries;
3. removes the constant synthesis adapters while preserving every designed
   construct;
4. attaches full-precision replicate measurements, published targets, folds,
   coordinates, and nullable activity calls; and
5. validates row identity, sequence geometry, replicate filtering, active
   counts, and the published target transformations.

No source workbooks or generated Parquet files are committed.

## Repository layout

```text
LentiMPRA/
|-- HF_DATASET_CARD.md
|-- PROCESSING_NOTES.md
|-- process.py
|-- README.md
|-- lentimpra-hepg2.parquet
|-- lentimpra-k562.parquet
|-- lentimpra-wtc11.parquet
|-- lentimpra-joint.parquet
|-- published-performance-large-scale.parquet
`-- published-performance-joint.parquet
scripts/
|-- build_hf_release.py
|-- common.py
|-- create_environment.sh
|-- download_data.py
|-- refresh_manifest.py
|-- release_staged_outputs.py
|-- run_all.sh
`-- verify_outputs.py
environment.yml
environment-lock-linux-64.txt
manifest.json
```

`manifest.json` is the single machine-readable authority for source identities,
the closed output allowlist, ordered Arrow schemas, scientific contracts,
fixed source-derived column fingerprints and row test vectors, environment
identity, processing revision, and Hugging Face payload. The task directory
owns assay processing and documentation; repository-wide download,
verification, release, and manifest maintenance stay under `scripts/`.

## Outputs

There are two dataset families:

- **Cell-specific large-scale libraries:** 164,307 HepG2, 243,780 K562, and
  75,542 WTC11 constructs with the same 23-column schema.
- **Joint library:** 60,000 constructs with separate HepG2, K562, and WTC11
  measurements, source element-specificity scores, three published targets,
  and one shared 10-fold assignment.

Every design row is retained. Null measurements record assay, QC, or modeling
exclusions rather than being converted to inactive labels.

## Hugging Face

The four public Parquets are published together as
[`Taykhoom/lentimpra-enhancer-activity`](https://huggingface.co/datasets/Taykhoom/lentimpra-enhancer-activity).

```python
from datasets import load_dataset

hepg2 = load_dataset(
    "Taykhoom/lentimpra-enhancer-activity",
    "cell-specific",
    split="hepg2",
)
joint = load_dataset(
    "Taykhoom/lentimpra-enhancer-activity",
    "joint",
    split="joint",
)
performance = load_dataset(
    "Taykhoom/lentimpra-enhancer-activity",
    "published-performance",
    split="large_scale",
)
```

The cell-line names are biological library splits, not train/test partitions.

The `published-performance` config preserves the paper's held-out-fold
Pearson correlations from Supplementary Tables S8 and S12. It has one row per
cell line, test fold, and published model: 180 large-scale rows and 150 joint
rows. These are aggregate metrics from the paper, not row-level predictions.

Mean Pearson correlation across the ten published folds:

| Library | Cell line | Biochemical | MPRAnn | MPRALegNet | Sei | Enformer |
|---|---|---:|---:|---:|---:|---:|
| Large-scale | HepG2 | 0.7320 | 0.7665 | 0.8333 | 0.7914 | 0.8306 |
| Large-scale | K562 | 0.7235 | 0.8082 | 0.8665 | 0.8150 | 0.8497 |
| Large-scale | WTC11 | 0.7124 | 0.7572 | 0.7946 | 0.8212 | 0.8335 |
| Joint | HepG2 | 0.6610 | 0.7255 | 0.7795 | 0.7891 | 0.8191 |
| Joint | K562 | 0.5738 | 0.6780 | 0.7465 | 0.7457 | 0.7811 |
| Joint | WTC11 | 0.6057 | 0.7139 | 0.7668 | 0.7820 | 0.8011 |

The large-scale sheet also reports a universal-feature biochemical model,
with mean correlations of 0.7365, 0.7404, and 0.7458 for HepG2, K562, and
WTC11 respectively.

## Reproducibility

The source archive is pinned by URL, byte size, and SHA-256:

```text
aefb1974e9ebe5c196a89a00f4e868c0abaf69e3d3952e79c74bfb42f217f6a4
```

The processing uses Supplementary Tables S3, S4, S8, S10, S11, and S12. The
exact table relationships and recovered published transformations are
documented in
[LentiMPRA/PROCESSING_NOTES.md](LentiMPRA/PROCESSING_NOTES.md).

## Environment and normal build

```bash
bash scripts/create_environment.sh
bash scripts/run_all.sh
```

`run_all.sh` invokes every Python command through the named
`lentimpra-enhancer-processing` environment and verifies its complete
`linux-64` artifact URL/SHA-256 set against
`environment-lock-linux-64.txt`, including the installed `conda-meta`
records.
`environment.yml` remains the portable direct-dependency specification. There
is no `requirements-lock.txt` because no dependency is installed with pip.

Normal processing builds all six Parquets in fresh repository-local staging,
independently verifies them against the committed root manifest, and only then
installs the manifest-matching bytes. It never updates `manifest.json` and
cannot bless changed output hashes or schemas.

To reuse the already downloaded archive:

```bash
bash scripts/run_all.sh \
  --reuse-archive /path/to/41586_2024_8430_MOESM5_ESM.zip
```

An alternate data root is propagated through download, processing, and
verification:

```bash
bash scripts/run_all.sh \
  --data-root /path/to/lentimpra-data \
  --reuse-archive /path/to/41586_2024_8430_MOESM5_ESM.zip
```

## Explicit manifest refresh

An intentional scientific or processing change uses a separate, nonpublishing
candidate workflow:

```bash
mamba run --name lentimpra-enhancer-processing \
  python scripts/refresh_manifest.py \
  --candidate-manifest /review/path/lentimpra-manifest-candidate.json
```

The command rebuilds in fresh staging and independently verifies the candidate,
including a complete source-workbook reconstruction of every output column,
but does not replace `manifest.json` or published Parquets. A normal refresh
preserves the tracked manifest's fixed source-derived per-column fingerprints
and declarative row vectors, so coordinated producer-and-oracle drift is still
rejected. Exact source identifier mapping is used for published S11 ESS values
because their upstream normalization is unavailable. After review, the
maintainer explicitly replaces the root manifest in the same change as the
processing update. Until then, normal commands reject the stale manifest
contract.

Only a reviewed scientific contract change may rebaseline those fixed values,
using the deliberately narrow break-glass option:

```bash
mamba run --name lentimpra-enhancer-processing \
  python scripts/refresh_manifest.py \
  --candidate-manifest /review/path/lentimpra-manifest-candidate.json \
  --maintainer-break-glass-rebaseline-source-output-semantics-after-review
```

## Hugging Face release

Build a verified local Hugging Face checkout:

```bash
mamba run --name lentimpra-enhancer-processing \
  python scripts/build_hf_release.py \
  --output /path/to/lentimpra-enhancer-activity
```

The builder validates an exact seven-file allowlist—one dataset card and six
Parquets—in a fresh sibling staging directory, then transactionally replaces
the target so stale files cannot mix with the new release. It never publishes
`manifest.json`. Stable `main` links in the card point back to the processing
repository, task documentation, and root manifest without embedding a commit.
The target must not overlap the processing repository. An existing target is
accepted only when it is empty and disposable or is the clean root of a Git
repository containing exactly the validated seven-file payload plus `.git`
and `.gitattributes`. Hugging Face's standard managed attributes template is
accepted when its effective Parquet attributes remain
`filter=lfs diff=lfs merge=lfs -text`; missing rules and later conflicting
overrides are rejected. The attributes file remains Git metadata outside the
Hub payload and is preserved byte-for-byte across replacement; arbitrary
non-empty directories are refused. Interrupted publication validates staged
payload bytes before promotion or backup deletion and restores a valid backup
when staging is invalid.

## Source details

Each source oligo is 230 nt: a constant 15 nt adapter, the 200 nt assayed
insert, and a second constant 15 nt adapter. The Parquets store only the
assayed insert.

The three large-scale files preserve orientation-specific measurements. Their
published S8 targets are orientation-averaged and standardized, so both
orientations of one element share a fold and target. The joint file preserves
the released per-cell activity and element-specificity measurements and the
published S12 three-target modeling subset.

## Citation

Agarwal, V. et al. *Massively parallel characterization of transcriptional
regulatory elements*. Nature (2025).
https://doi.org/10.1038/s41586-024-08430-9

## License

The processing code is MIT licensed. The article and supplementary data are
CC BY 4.0; see [NOTICE.md](NOTICE.md).
