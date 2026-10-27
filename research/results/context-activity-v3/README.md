# Published reporter activity: frozen v3

Numerical reproduction: **true**
(258 recorded checks). Literal published-activity binding:
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
python3 scripts/package_context_activity.py --verify research/results/context-activity-v3
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


`core/` retains every schema-5 Rust fit, scaler, prediction, common-activity
identity, fold/group assignment, paired bootstrap and source receipt.
`independent-python/` retains the complete independent implementation and outputs.
`comparison/` retains every source/numerical check; `figures/` retains standalone
figures with numeric tables and the exact parent-result commitment.
`recovery-publication.json` binds the separate complete recovery snapshot;
verify that snapshot separately to inspect its full original arithmetic outputs.
The retained earlier input-binding activity payload is byte-identical to the current activity payload. Earlier inventory/36-route validation bindings were superseded; that hardening did not replace measured values.


Rust result code commit: `035ca5731fe594622cd2c3c27cd3bb51ec65c924`; working-tree dirty flag:
`False`. These are the values recorded at the
original execution, not claims about the current checkout.

Audit the Rust evidence and replay only the saved numerical comparison:

```bash
source scripts/research-env.sh
cargo run --release --frozen -- audit research/results/context-activity-v3/core
python3 scripts/replicate_context_activity.py \
  --compare-rust research/results/context-activity-v3/core/study.json.gz \
  --independent-run research/results/context-activity-v3/independent-python \
  --out runs/context-activity-published-comparison
```

The replay output directory must be new. Numerical agreement does not establish
normalization history, biological replication or a unique regulatory mechanism.
All five models and all declared comparisons remain available regardless of
the sign of the primary C-versus-A improvement.
