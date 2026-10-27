#!/usr/bin/env python3
"""Package recovered-source and activity evidence; verify with the standard library."""
import argparse
import itertools
import json
from pathlib import Path
import shlex
import shutil

from package_context_study import Files, deterministic_tar, relative_name, require


RECOVERY_SHA = "00175e41cf9fa0cbf177ce8fc7ef143eb508aec8b11d058f7acb4458284dc868"
ACTIVITY_SHA = "64e14037fb1f42ddc1ecaca1d23bb8096fabac1f1e892d105eaa20f525fbbb20"
ENDPOINT = "published_joint_log2_rna_dna_activity"
RAW_SHA = "9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e"
RAW_BYTES = 53_557_805
MAX_FILE_BYTES = 100_000_000
ROOT = Path(__file__).resolve().parents[1]
RAW_PATH = ROOT / "datasets/context-endpoint-recovery/raw/jointdata.xlsx"
GATES = ("archive_integrity_verified", "empirical_activity_schema_bound", "coverage_complete",
         "design_identity_verified", "replicate_aggregation_status_recorded")
CORE_FILES = {"activity.json.gz", "archive-validation.json", "common-activity.json.gz",
              "extraction-receipt.json", "features.json.gz", "findings.md", "grouping.json",
              "inputs.json", "ledger.json", "protocol.json", "recovery-summary.json",
              "resampling.json.gz", "source-binding.json", "source-export.json.gz",
              "source-ledger.jsonl.gz", "source-manifest.json", "study.json.gz", "summary.json"}


def checked_manifest(files, name, schema=None):
    document = files.json(name)
    if schema is not None:
        require(document.get("schema_version") == schema, "unexpected manifest schema: " + str(files.root))
    expected = document["files"]
    require(isinstance(expected, dict) and expected, "empty evidence manifest")
    require(set(files.files) == set(expected) | {name}, "manifest inventory differs: " + str(files.root))
    for filename, commitment in expected.items():
        relative_name(filename)
        actual = files.files.get(filename)
        if isinstance(commitment, str):
            good = actual is not None and actual["sha256"] == commitment
        else:
            good = actual == commitment
        require(good, "artifact hash/size differs: " + str(files.root) + "/" + filename)
    return document


def matched_comparison(files, protocol_sha):
    manifest = checked_manifest(files, "manifest.json", 1)
    comparison = files.json("comparison.json")
    require(manifest.get("protocol_sha256") == comparison.get("protocol_sha256") == protocol_sha,
            "comparison protocol differs")
    require(comparison.get("matched") is True and comparison.get("failure_count") == 0
            and comparison.get("failures") == [] and comparison.get("checks")
            and all(check.get("matched") is True for check in comparison["checks"]),
            "complete numerical reproduction has not passed")
    return comparison


def check_raw_source(path=RAW_PATH):
    require(path.is_file() and not path.is_symlink(), "restore the frozen raw workbook: " + str(path))
    # Streaming avoids loading the 53.5 MB workbook; it is never interpreted here.
    import hashlib
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    require(path.stat().st_size == RAW_BYTES and result.hexdigest() == RAW_SHA,
            "frozen raw Table11 byte identity differs")


def source_identity(source):
    return source.get("sha256") == RAW_SHA and source.get("bytes") == RAW_BYTES


def recovery_status(primary, independent, comparison_files):
    primary_manifest = checked_manifest(primary, "artifacts.json", 1)
    independent_manifest = checked_manifest(independent, "manifest.json")
    require({"binding.json", "implementation.py", "independent-source-values.json.gz", "inputs.json",
             "old-source-values.json.gz", "protocol.json", "source-parser-comparison.json", "summary.json"}
            <= set(independent_manifest["files"]), "independent recovery omits required source/result artifacts")
    protocol = primary.json("protocol.json")
    summary, inputs = primary.json("summary.json"), primary.json("inputs.json")
    binding, own_inputs = primary.json("schema-binding.json"), independent.json("inputs.json")
    require(primary.sha("protocol.json") == RECOVERY_SHA
            and primary_manifest.get("protocol_sha256") == summary.get("protocol_sha256")
            == inputs.get("protocol_sha256") == own_inputs.get("protocol_sha256") == RECOVERY_SHA,
            "recovery protocol binding differs")
    if "protocol_sha256" in independent_manifest:
        require(independent_manifest["protocol_sha256"] == RECOVERY_SHA, "independent recovery protocol differs")
    require(source_identity(summary["source"]) and source_identity(own_inputs["source"])
            and source_identity(binding["source_workbook"])
            and binding.get("semantic_binding_complete") is True
            and binding.get("status") == "bound_before_quantitative_inspection"
            and binding.get("protocol_sha256") == RECOVERY_SHA
            and own_inputs.get("schema_binding_sha256") == primary.sha("schema-binding.json"),
            "recovery raw source/schema commitments differ")
    expected_variants = protocol["fixed_pipeline_grid"]["variant_count_if_all_fields_present"]
    variants = summary["variants"]
    require(len(variants) == summary["variant_count"] == expected_variants
            and len({v["id"] for v in variants}) == expected_variants,
            "recovery omits frozen arithmetic routes")
    tuples = {(v["activity_base"], v["pre_zscore_population"], v["pre_zscore_ddof"], v["post_ESS_step"])
              for v in variants}
    expected = set(itertools.product(protocol["fixed_pipeline_grid"]["base_routes"],
        ("complete_three_context_source", "all_finite_per_context"), (0, 1),
        ("none", "column_zscore_ddof0", "column_zscore_ddof1")))
    require(tuples == expected, "recovery Cartesian route definitions differ")
    own_summary = independent.json("summary.json")
    populations = {"complete3contextsourceintersection": "complete_three_context_source",
                   "allfinite_independently_percell": "all_finite_per_context"}
    own_tuples = {(v["base_id"], populations[v["pre_population"]], v["pre_ddof"],
                   "none" if v["post_ddof"] is None else "column_zscore_ddof" + str(v["post_ddof"]))
                  for v in own_summary["pipelines"]}
    require(len(own_summary["pipelines"]) == expected_variants and own_tuples == expected
            and independent.sha("protocol.json") == RECOVERY_SHA
            and independent.sha("binding.json") == primary.sha("schema-binding.json")
            and own_summary.get("predictive_fitting_performed") is False
            and own_summary.get("historical_pipeline_identified") is False
            and own_summary.get("source_values_or_old_targets_modified") is False,
            "independent recovery grid/protocol/status differs")
    parser = independent.json("source-parser-comparison.json")
    require(parser.get("all_source_rows_and_lexemes_exact") is True
            and parser.get("differing_tables") == []
            and parser.get("primary_export_sha256") == primary.sha("source-values.json.gz"),
            "independent raw source parsing differs")
    for variant in variants:
        if variant["estimable"]:
            require(variant["derived_values_file"] in primary.files, "primary route array missing")
    for variant in own_summary["pipelines"]:
        if variant["status"] == "estimable":
            require(variant["values_file"] in independent.files, "independent route array missing")
    require(inputs["implementation_sha256"] == primary.sha("implementation.py")
            and inputs["old_s9_reader_sha256"] == primary.sha("old-s9-reader.py"),
            "recovery implementation snapshots differ")
    hypothesis = summary["pre_value_code_hypothesis"]
    tests = hypothesis.get("fixed_1e-12_test_by_cohort", {})
    full = tests.get("full_s12_source_id_set", {})
    require(hypothesis.get("estimable") is True
            and full.get("ids") == protocol["old_inputs"]["s9_workbook"]["nonempty_unique_ids"]
            and full.get("matched") is True and full.get("maximum_absolute_error", 1) <= 1e-12
            and all(test.get("matched") is True for test in tests.values())
            and summary["flags"].get("historical_pipeline_identified") is False
            and summary["flags"].get("endpoint_replaced") is False
            and summary["flags"].get("predictive_models_fitted") is False,
            "recovery source relationship/status differs")
    comparison = matched_comparison(comparison_files, RECOVERY_SHA)
    require(comparison["primary_artifacts_manifest"]["sha256"] == primary.sha("artifacts.json")
            and comparison["independent_manifest"]["sha256"] == independent.sha("manifest.json"),
            "recovery comparison identifies different result manifests")
    return {"schema_version": 1, "kind": "context-endpoint-recovery-v1", "protocol_sha256": RECOVERY_SHA,
        "numeric_reproduction": comparison["matched"], "numerical_checks": len(comparison["checks"]),
        "source_recovery": {"s11_to_s12_sample_zscore_compatible": full["matched"],
            "complete_s12_ids": full["ids"], "maximum_absolute_error": full["maximum_absolute_error"]},
        "literal_activity_binding_verified": binding["semantic_binding_complete"],
        "upstream_ess_history_identified": summary["flags"]["historical_pipeline_identified"],
        "endpoint_replaced": summary["flags"]["endpoint_replaced"], "predictive_models_fitted": False,
        "independent_biological_validation": False,
        "primary_artifacts_sha256": primary.sha("artifacts.json"),
        "primary_summary_sha256": primary.sha("summary.json"),
        "independent_manifest_sha256": independent.sha("manifest.json"),
        "primary_runtime": inputs["runtime"], "independent_runtime": own_inputs.get("runtime"),
        "source": {"repository_path": "datasets/context-endpoint-recovery/raw/jointdata.xlsx",
                   "bytes": RAW_BYTES, "sha256": RAW_SHA},
        "immutable_previous_protocols": protocol["immutable_previous_protocols"],
        "arithmetic_grid_routes": len(variants),
        "limits": ["All original arithmetic routes and their divergences are retained.",
                   "Source-file compatibility does not identify original ESS normalization history.",
                   "Final release and old source represent one experiment, not independent biological replication."]}


def raw_binding(comparison, key, files, name):
    receipt = comparison[key]
    require(receipt.get("raw_sha256") == files.sha(name)
            and receipt.get("raw_bytes") == files.files[name]["bytes"], "raw comparison binding differs: " + key)


def activity_status(core, independent, comparison_files, figures, recovery_publication):
    artifacts = checked_manifest(core, "artifacts.json", 5)
    require(CORE_FILES <= set(artifacts["files"]), "schema-5 core omits required evidence")
    own_manifest = checked_manifest(independent, "manifest.json", 1)
    protocol, inputs, summary = core.json("protocol.json"), core.json("inputs.json"), core.json("summary.json")
    source = core.json("source-manifest.json")
    require(core.sha("protocol.json") == ACTIVITY_SHA
            and inputs.get("protocol_sha256") == summary.get("protocol_sha256")
            == source.get("protocol_sha256") == own_manifest.get("protocol_sha256") == ACTIVITY_SHA
            and protocol["required_recovery_protocol_sha256"] == RECOVERY_SHA,
            "activity protocol binding differs")
    require(inputs.get("endpoint") == summary.get("endpoint") == source.get("endpoint")
            == own_manifest.get("endpoint") == ENDPOINT,
            "activity endpoint differs from unchanged published activity")
    require(source_identity(source["source"]) and inputs.get("raw_source_sha256") == RAW_SHA
            and inputs.get("raw_source_bytes") == RAW_BYTES
            and all(source["gates"].get(gate) is True for gate in GATES)
            and summary.get("source_gate_passed") is True,
            "activity source/schema/coverage gates have not passed")
    require(inputs.get("source_manifest_sha256") == core.sha("source-manifest.json")
            and inputs.get("activity_sha256") == core.sha("activity.json.gz")
            and inputs.get("grouping_sha256") == core.sha("grouping.json") == protocol["grouping"]["sha256"]
            and inputs.get("dna_panel_sha256") == protocol["panel"]["sha256"],
            "activity source/activity/DNA/grouping byte bindings differ")
    for name, expected in source["files"].items():
        require(core.sha(name) == expected, "archived activity source receipt differs: " + name)
    require(core.json("source-binding.json") == source["schema_binding"] == summary["source_schema_binding"]
            and core.json("recovery-summary.json") == source["aggregation_status"] == summary["aggregation_status"],
            "source schema/aggregation statuses were changed")
    require(recovery_publication.get("kind") == "context-endpoint-recovery-v1"
            and recovery_publication.get("protocol_sha256") == RECOVERY_SHA
            and recovery_publication.get("numeric_reproduction") is True
            and recovery_publication.get("upstream_ess_history_identified") is False
            and recovery_publication["primary_summary_sha256"] == core.sha("recovery-summary.json")
            and recovery_publication["primary_artifacts_sha256"] == core.sha("recovery-artifacts.json"),
            "activity source recovery identifies different original arithmetic evidence")
    own_inputs, own_status = independent.json("inputs.json"), independent.json("status.json")
    require(own_status.get("status") == "complete"
            and own_inputs.get("endpoint") == own_status.get("endpoint") == ENDPOINT
            and own_inputs.get("protocol_sha256") == own_status.get("protocol_sha256") == ACTIVITY_SHA
            and own_inputs.get("manifest_sha256") == core.sha("source-manifest.json")
            and own_inputs.get("source_workbook_sha256") == RAW_SHA
            and own_inputs.get("source_workbook_bytes") == RAW_BYTES
            and own_inputs.get("grouping_sha256") == inputs["grouping_sha256"]
            and own_inputs.get("prior_metadata_panel_sha256") == inputs["dna_panel_sha256"]
            and own_inputs.get("ESS_observations_used") is False,
            "independent activity run is incomplete or identifies different sources")
    require(own_inputs.get("implementation_sha256") == independent.sha("implementation.py")
            and own_inputs.get("independent_base_sha256") == independent.sha("independent-context-base.py"),
            "independent activity implementation bytes differ")
    comparison = matched_comparison(comparison_files, ACTIVITY_SHA)
    require(comparison.get("endpoint") == ENDPOINT
            and comparison["rust_artifacts_manifest"].get("sha256") == core.sha("artifacts.json")
            and comparison["independent_manifest"].get("sha256") == independent.sha("manifest.json")
            and comparison.get("implementation_sha256") == comparison_files.sha("implementation.py")
            and comparison.get("independent_base_sha256") == comparison_files.sha("independent-context-base.py"),
            "activity comparison parent/implementation bindings differ")
    for key, files, name in (("rust_study", core, "study.json.gz"), ("rust_features", core, "features.json.gz"),
            ("rust_resampling", core, "resampling.json.gz"), ("rust_common_activity", core, "common-activity.json.gz"),
            ("independent_study", independent, "evaluation.json.gz"),
            ("independent_resampling", independent, "resampling.json.gz"),
            ("independent_raw_source_records", independent, "source-records.json.gz")):
        raw_binding(comparison, key, files, name)
    figure_manifest = checked_manifest(figures, "manifest.json", 1)
    require(figure_manifest.get("endpoint") == ENDPOINT
            and figure_manifest.get("source_artifacts_sha256") == core.sha("artifacts.json")
            and figure_manifest.get("fits_performed") is False,
            "figures identify another activity study or fitted responses")
    require(summary["sequence_count"] == own_status["elements"] == protocol["panel"]["elements"]
            and summary["measurement_count"] == own_status["observations"] == protocol["panel"]["observations"],
            "activity observation coverage differs")
    return {"schema_version": 1, "kind": "context-activity-v3", "endpoint": ENDPOINT,
        "protocol_sha256": ACTIVITY_SHA, "required_recovery_protocol_sha256": RECOVERY_SHA,
        "immutable_previous_protocols": recovery_publication["immutable_previous_protocols"],
        "numeric_reproduction": comparison["matched"], "numerical_checks": len(comparison["checks"]),
        "source_recovery": recovery_publication["source_recovery"],
        "literal_activity_binding_verified": summary["source_gate_passed"],
        "upstream_ess_history_identified": recovery_publication["upstream_ess_history_identified"],
        "independent_biological_validation": False, "source_normalization_history_proven": False,
        "software": inputs["software"],
        "source": {"repository_path": "datasets/context-endpoint-recovery/raw/jointdata.xlsx", "bytes": RAW_BYTES, "sha256": RAW_SHA},
        "counts": {key: summary[key] for key in ("sequence_count", "measurement_count", "group_count", "contexts")},
        "core_artifacts_sha256": core.sha("artifacts.json"),
        "published_activity_sha256": core.sha("activity.json.gz"),
        "independent_manifest_sha256": independent.sha("manifest.json"),
        "source_manifest_sha256": core.sha("source-manifest.json"),
        "recovery_primary_artifacts_sha256": recovery_publication["primary_artifacts_sha256"],
        "recovery_primary_summary_sha256": recovery_publication["primary_summary_sha256"],
        "limits": ["Unchanged published library-normalized log2(RNA/DNA), not absolute transcription.",
                   "Same selected experiment, elements and measured cells; no independent biological replication or new mechanism.",
                   "Numerical reproduction and literal activity binding do not prove upstream ESS or source normalization history."]}


def prepare_out(out, sources):
    require(not out.exists() and not out.is_symlink(), "output exists; preserve previous snapshots")
    for source in sources:
        require(source.resolve() != out.resolve() and source.resolve() not in out.resolve().parents,
                "publication output cannot be inside an input directory")
        files = Files(source)
        require(all(entry["bytes"] < MAX_FILE_BYTES for entry in files.files.values()),
                "input contains a >=100MB file; preserve it as compressed evidence before publication")


def copy_evidence(source, out):
    files = Files(source)
    require(all(entry["bytes"] < MAX_FILE_BYTES for entry in files.files.values()),
            "input has a >=100MB file; retain its bytes in a compressed evidence artifact before packaging")
    shutil.copytree(source, out)


def checked_addendum(files, expected_summary_sha, parent_protocol_sha):
    name = "artifacts.json" if "artifacts.json" in files.files else "manifest.json"
    document = checked_manifest(files, name)
    record = {"manifest": name, "manifest_sha256": files.sha(name)}
    if "protocol.json" in files.files:
        protocol = files.json("protocol.json")
        require(document.get("protocol_sha256") == files.sha("protocol.json")
                and protocol.get("status") == "frozen"
                and protocol.get("parent_protocol_sha256") in (RECOVERY_SHA, parent_protocol_sha),
                "addendum frozen protocol/parent binding differs")
        inputs = protocol.get("inputs", {})
        if "primary_summary" in inputs:
            require(inputs["primary_summary"]["sha256"] == expected_summary_sha,
                    "addendum refers to another original recovery summary")
        if "official_jointdata" in inputs:
            require(source_identity(inputs["official_jointdata"]), "addendum source workbook differs")
        record["protocol_sha256"] = files.sha("protocol.json")
        record["parent_protocol_sha256"] = protocol["parent_protocol_sha256"]
    return record


def history_and_addenda(out, history, addenda, status):
    record = {"history": {}, "addenda": {}}
    for category, paths in (("history", history), ("addenda", addenda)):
        for source in paths:
            relative_name(source.name)
            require(source.name not in record[category], "duplicate evidence-directory basename")
            if category == "addenda":
                files = Files(source)
                record[category][source.name] = checked_addendum(files,
                    status.get("primary_summary_sha256", status.get("recovery_primary_summary_sha256")), status["protocol_sha256"])
            else:
                historical = Files(source)
                record[category][source.name] = {
                    "role": "retained attempt or input-binding history; not an alternative successful endpoint",
                    "activity_payloads": {name: {"sha256": item["sha256"],
                        "identical_to_current_activity": item["sha256"] == status["published_activity_sha256"]
                            if status.get("published_activity_sha256") else None}
                        for name, item in historical.files.items() if Path(name).name == "activity.json.gz"}}
            copy_evidence(source, out / category / source.name)
    return record


def addendum_relations(out, retained):
    """Bind each supplied numerical addendum comparison to both included runs."""
    included = {record["manifest_sha256"] for record in retained["addenda"].values()}
    results = {}
    for name in retained["addenda"]:
        files = Files(out / "addenda" / name)
        if "comparison.json" not in files.files:
            continue
        document = files.json("comparison.json")
        comparison = matched_comparison(files, document["protocol_sha256"])
        primary = comparison["primary_artifacts_manifest"]["sha256"]
        independent = comparison["independent_manifest"]["sha256"]
        require(primary in included and independent in included and primary != independent,
                "addendum comparison does not identify both included source runs")
        results[name] = {"protocol_sha256": comparison["protocol_sha256"],
            "matched": comparison["matched"], "checks": len(comparison["checks"]),
            "primary_artifacts_sha256": primary, "independent_manifest_sha256": independent,
            "metadata_label_adapters": comparison.get("metadata_label_adapters", {})}
    return results


def finish(out, status, history, addenda):
    status["retained_evidence"] = history_and_addenda(out, history, addenda, status)
    status["addendum_numerical_reproductions"] = addendum_relations(out, status["retained_evidence"])
    shutil.copyfile(__file__, out / "packaging-implementation.py")
    shutil.copyfile(Path(__file__).with_name("package_context_study.py"), out / "packaging-helpers.py")
    (out / "README.md").write_text(readme(out, status), encoding="utf-8")
    status["files"] = Files(out).files
    (out / "publication.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    verify(out)


def recovery(args):
    paths = [args.primary, args.independent, args.comparison, *args.history, *args.addendum]
    prepare_out(args.out, paths)
    check_raw_source()
    status = recovery_status(Files(args.primary), Files(args.independent), Files(args.comparison))
    args.out.mkdir(parents=True)
    for source, target in ((args.primary, "primary"), (args.independent, "independent"), (args.comparison, "comparison")):
        if target == "independent":
            deterministic_tar(Files(source), args.out / "independent-recovery.tar.gz")
        else:
            copy_evidence(source, args.out / target)
    finish(args.out, status, args.history, args.addendum)


def activity(args):
    paths = {key: getattr(args, key) or Path(str(args.prefix) + "-" + key)
             for key in ("rust", "python", "comparison", "figures")}
    prepare_out(args.out, [*paths.values(), args.recovery, *args.history, *args.addendum])
    verify(args.recovery)
    recovered = json.loads((args.recovery / "publication.json").read_bytes())
    status = activity_status(*(Files(paths[key]) for key in ("rust", "python", "comparison", "figures")), recovered)
    args.out.mkdir(parents=True)
    for key, target in (("rust", "core"), ("python", "independent-python"), ("comparison", "comparison"), ("figures", "figures")):
        copy_evidence(paths[key], args.out / target)
    shutil.copyfile(args.recovery / "publication.json", args.out / "recovery-publication.json")
    finish(args.out, status, args.history, args.addendum)


def verify(out):
    files = Files(out)
    publication = files.json("publication.json")
    require(publication.get("files") == {k: v for k, v in files.files.items() if k != "publication.json"},
            "publication inventory, size or SHA-256 differs")
    require(all(v["bytes"] < MAX_FILE_BYTES for v in files.files.values()), "published file reaches 100MB")
    if publication.get("kind") == "context-endpoint-recovery-v1":
        status = recovery_status(Files(out / "primary"), Files(out / "independent-recovery.tar.gz", archive=True), Files(out / "comparison"))
    elif publication.get("kind") == "context-activity-v3":
        status = activity_status(Files(out / "core"), Files(out / "independent-python"), Files(out / "comparison"),
                                 Files(out / "figures"), files.json("recovery-publication.json"))
    else:
        raise ValueError("unsupported publication kind")
    retained = publication["retained_evidence"]
    for name, record in retained["history"].items():
        relative_name(name)
        historical = Files(out / "history" / name)
        expected_activity = {path: {"sha256": item["sha256"],
            "identical_to_current_activity": item["sha256"] == status.get("published_activity_sha256")
                if status.get("published_activity_sha256") else None}
            for path, item in historical.files.items() if Path(path).name == "activity.json.gz"}
        require(record["activity_payloads"] == expected_activity, "historical activity byte-identity status differs")
    for name, record in retained["addenda"].items():
        addendum = Files(out / "addenda" / name)
        require(checked_addendum(addendum,
            status.get("primary_summary_sha256", status.get("recovery_primary_summary_sha256")), status["protocol_sha256"]) == record,
            "addendum protocol/manifest/parent binding differs")
    status["retained_evidence"] = retained
    status["addendum_numerical_reproductions"] = addendum_relations(out, retained)
    require(status == {k: v for k, v in publication.items() if k != "files"},
            "published machine status differs from bound evidence")
    check_raw_source()
    print(json.dumps({"publication_verified": True, "kind": status["kind"], "files": len(publication["files"]),
                      "numeric_reproduction": status["numeric_reproduction"],
                      "literal_activity_binding_verified": status["literal_activity_binding_verified"],
                      "upstream_ess_history_identified": status["upstream_ess_history_identified"]}))


def readme(out, status):
    target = shlex.quote(out.as_posix())
    shared = f"""Numerical reproduction: **{str(status['numeric_reproduction']).lower()}**
({status['numerical_checks']} recorded checks). Literal published-activity binding:
**{str(status['literal_activity_binding_verified']).lower()}**. The source S11-to-S12
sample-zscore relationship is compatible; upstream ESS history identified:
**{str(status['upstream_ess_history_identified']).lower()}**.

These flags are separate claims. No global scientific-success flag is implied.
This is source recovery and computational reanalysis of one selected experiment,
not an independent study, new biological mechanism, unseen-cell prediction,
clinical validation or absolute transcription calibration. Original published
values and all frozen earlier results remain unchanged.

Verify byte integrity and cross-artifact commitments from the repository root
without scientific Python dependencies:

```bash
python3 scripts/package_context_activity.py --verify {target}
```

The raw source is retained separately in repository
`datasets/context-endpoint-recovery/raw/jointdata.xlsx` (SHA-256 `{RAW_SHA}`,
{RAW_BYTES:,} bytes). The verifier checks it but does not interpret it. The
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
"""
    if status["kind"] == "context-endpoint-recovery-v1":
        return "# Recovered context-source evidence\n\n" + shared + """

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
python3 scripts/audit_recovered_context.py \\
  --binding """ + target + "/primary/schema-binding.json \\\n  --source " + target + "/primary/source-values.json.gz \\\n  --out runs/recovered-context-replay\n```\n"
    software = status["software"]
    history_payloads = [entry for record in status["retained_evidence"]["history"].values()
                        for entry in record["activity_payloads"].values()]
    history_note = ("The retained earlier input-binding activity payload is byte-identical to "
                    "the current activity payload. Earlier inventory/36-route validation bindings "
                    "were superseded; that hardening did not replace measured values.\n"
                    if history_payloads and all(item["identical_to_current_activity"] for item in history_payloads) else "")
    return "# Published reporter activity: frozen v3\n\n" + shared + f"""

`core/` retains every schema-5 Rust fit, scaler, prediction, common-activity
identity, fold/group assignment, paired bootstrap and source receipt.
`independent-python/` retains the complete independent implementation and outputs.
`comparison/` retains every source/numerical check; `figures/` retains standalone
figures with numeric tables and the exact parent-result commitment.
`recovery-publication.json` binds the separate complete recovery snapshot;
verify that snapshot separately to inspect its full original arithmetic outputs.
{history_note}

Rust result code commit: `{software.get('commit')}`; working-tree dirty flag:
`{software.get('working_tree_dirty')}`. These are the values recorded at the
original execution, not claims about the current checkout.

Audit the Rust evidence and replay only the saved numerical comparison:

```bash
source scripts/research-env.sh
cargo run --release --frozen -- audit {target}/core
python3 scripts/replicate_context_activity.py \\
  --compare-rust {target}/core/study.json.gz \\
  --independent-run {target}/independent-python \\
  --out runs/context-activity-published-comparison
```

The replay output directory must be new. Numerical agreement does not establish
normalization history, biological replication or a unique regulatory mechanism.
All five models and all declared comparisons remain available regardless of
the sign of the primary C-versus-A improvement.
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", type=Path)
    commands = parser.add_subparsers(dest="command")
    recover = commands.add_parser("recovery")
    for name in ("primary", "independent", "comparison", "out"):
        recover.add_argument("--" + name, required=True, type=Path)
    recover.add_argument("--history", action="append", type=Path, default=[])
    recover.add_argument("--addendum", action="append", type=Path, default=[])
    recover.set_defaults(action=recovery)
    active = commands.add_parser("activity")
    for name in ("prefix", "recovery", "out"):
        active.add_argument("--" + name, required=True, type=Path)
    for name in ("rust", "python", "comparison", "figures"):
        active.add_argument("--" + name, type=Path)
    active.add_argument("--history", action="append", type=Path, default=[])
    active.add_argument("--addendum", action="append", type=Path, default=[])
    active.set_defaults(action=activity)
    args = parser.parse_args()
    if args.verify is not None:
        if args.command is not None:
            parser.error("--verify cannot be combined with publication creation")
        verify(args.verify)
    elif args.command:
        args.action(args)
    else:
        parser.error("provide --verify, recovery, or activity")


if __name__ == "__main__":
    main()
