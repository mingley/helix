# Recovered context-source evidence

Numerical reproduction: **true**
(7008 recorded checks). Literal published-activity binding:
**true**. The source S11-to-S12
sample-zscore relationship is compatible; upstream ESS history identified:
**false**.

These flags are separate claims. No global scientific-success flag is implied.
This is source recovery and computational reanalysis of one selected experiment,
not an independent study, new biological mechanism, unseen-cell prediction,
clinical validation or absolute transcription calibration. Original published
values and all frozen earlier results remain unchanged.

Verify byte integrity and cross-artifact commitments from the repository root
without scientific Python dependencies:

```bash
python3 scripts/package_context_activity.py --verify research/results/context-endpoint-recovery-v1
```

The raw source is retained separately in repository
`datasets/context-endpoint-recovery/raw/jointdata.xlsx` (SHA-256 `9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e`,
53,557,805 bytes). The verifier checks it but does not interpret it. The
342 MB official outer archive is not duplicated in this snapshot; acquisition,
member CRC and full-archive receipts remain attached to the source evidence.
The roughly 800 MB hg38 cache remains separate for full reference revalidation.

`publication.json` inventories every other published file by SHA-256 and byte
size. Current result and source manifests are checked; historical manifests
retain their original superseded or failed status. Existing compressed JSON files
retain all source values, fits, predictions or bootstrap draws; no scientific
outputs are reduced to a summary for this publication. `packaging-implementation.py`
and `packaging-helpers.py` preserve the packaging code. `history/`, when present,
retains supplied earlier attempts and stale input bindings as historical evidence.
Their filenames and statuses do not promote them to successful endpoints.
`addenda/`, when present, contains separately versioned post-hoc diagnostics and
does not replace the original 36-route audit or change the activity endpoint.

Source attribution: Agarwal et al., *Massively parallel characterization of
transcriptional regulatory elements*, Nature (2025),
[DOI 10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9).
Retain the article/supplementary-data CC BY 4.0 attribution and source-version
receipts; retained third-party code keeps its own license. Repository code
licensing does not relicense source data. Hashes establish consistency with
commitments, not publisher-signed authenticity or empirical truth.


`primary/` retains the complete original arithmetic grid, source tokens,
ID/missingness ledgers, code, parameters, derived arrays and every comparison.
`independent-recovery.tar.gz` retains every separately parsed/reconstructed
output, including the complete large summary, under archive root
`independent-python/`. It is deterministic (sorted members, zero timestamps)
and can be unpacked with `tar -xzf independent-recovery.tar.gz`.
`comparison/` binds both manifests and preserves all numerical checks.
The matching S11-to-S12 operation certifies a numerical relation between source
tables; it does not identify the unpublished activity-to-S11-ESS pipeline.

To replay source arithmetic, use fresh output directories and the repository's
pinned scientific environment. The exact source export and binding are retained:

```bash
source scripts/research-env.sh
python3 scripts/audit_recovered_context.py \
  --binding research/results/context-endpoint-recovery-v1/primary/schema-binding.json \
  --source research/results/context-endpoint-recovery-v1/primary/source-values.json.gz \
  --out runs/recovered-context-replay
```
