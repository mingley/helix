#!/usr/bin/env python3
"""Plot a sealed exploratory complete-case variant-transfer evaluation; no fitting."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import struct

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CELLS = ["K562", "HepG2"]
MODELS = ["zero", "I", "S"]
CLASSIFICATION = "exploratory_after_v1_source_failure"
TITLE = "Exploratory numeric complete-case analysis after v1 source failure"
PRODUCER_FILES = {"evaluation.json.gz", "summary.json", "normalized-input.json",
                  "protocol.json", "binding.json", "pins.json", "implementation.py"}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def finite(value):
    return type(value) in (float, int) and math.isfinite(value)


def load(path):
    def pairs(items):
        out = {}
        for key, value in items:
            require(key not in out, "duplicate JSON key")
            out[key] = value
        return out
    with Path(path).open() as stream:
        return json.load(stream, object_pairs_hook=pairs,
                         parse_constant=lambda _: (_ for _ in ()).throw(ValueError("nonfinite JSON")))


def verified_inputs(root):
    require(root.is_dir() and not root.is_symlink(), "evaluation directory missing or symlinked")
    require(not (root / "manifest.json").is_symlink(), "symlinked evaluation manifest")
    manifest = load(root / "manifest.json")
    require(manifest.get("schema_version") == 1
            and set(manifest.get("files", {})) == PRODUCER_FILES,
            "unexpected evaluator manifest inventory")
    # Verify the complete producer inventory before parsing any metrics or responses.
    for name, expected in manifest["files"].items():
        path = root / name
        require(path.is_file() and not path.is_symlink() and sha(path) == expected,
                "evaluator artifact commitment differs: " + name)
    protocol, summary, pins = (load(root / name) for name in
                               ["protocol.json", "summary.json", "pins.json"])
    require(protocol.get("protocol_id") == "variant-transfer-v2"
            and protocol.get("status") == "frozen"
            and protocol.get("analysis_classification") == CLASSIFICATION
            and summary.get("analysis_classification") == CLASSIFICATION,
            "only the explicitly exploratory v2 complete-case analysis may be plotted")
    require(summary.get("protocol_id") == "variant-transfer-v2"
            and summary.get("protocol_sha256") == pins.get("protocol_sha256")
            == manifest.get("input_pins", {}).get("protocol_sha256") == sha(root / "protocol.json")
            and pins == manifest.get("input_pins")
            and pins.get("input_sha256") == sha(root / "normalized-input.json")
            and pins.get("binding_sha256") == sha(root / "binding.json"),
            "producer input/protocol pins differ")
    require(protocol.get("operational_missingness", {}).get("rule")
            == "exclude_variant_unless_both_contexts_ordinary_finite_numeric"
            and protocol["operational_missingness"].get("source_semantics_established") is False,
            "missingness classification differs")
    with gzip.open(root / "evaluation.json.gz", "rt") as stream:
        evaluation = json.load(stream)
    normalized = load(root / "normalized-input.json")
    require(evaluation.get("schema_version") == 1 and evaluation.get("contexts") == CELLS
            and evaluation.get("model_order") == MODELS,
            "unexpected evaluation model/context schema")
    for name in ["observed", "primary", "per_component", "variant_count", "component_count"]:
        require(summary.get(name) == evaluation.get(name), "evaluation/summary differ: " + name)
    variants = sorted(normalized["variants"], key=lambda row: row["id"])
    ids = [row["id"] for row in variants]
    require(len(ids) == len(set(ids)) == evaluation["variant_count"]
            and ids == evaluation["variant_order"], "complete variant cohort differs")
    ledger = {row["id"]: row for row in evaluation["per_variant"]}
    require(set(ledger) == set(ids) and len(ledger) == len(evaluation["per_variant"]),
            "variant ledger membership differs")
    for row in variants:
        saved = ledger[row["id"]]
        require(all(row[k] == saved[k] for k in ["component", "locus", "gc_preserving"]),
                "variant metadata differs")
        for cell in CELLS:
            observed = row["observed"][cell]
            require(finite(observed) and struct.pack("<d", observed).hex()
                    == saved["observed_float64_little_endian"][cell], "observed float binding differs")
            for model in ["I", "S"]:
                predicted = row["predictions"][model][cell]
                require(finite(predicted) and struct.pack("<d", predicted).hex()
                        == saved["prediction_float64_little_endian"][model][cell],
                        "prediction float binding differs")
    components = evaluation["per_component"]
    require(len(components) == evaluation["component_count"]
            and len({row["component"] for row in components}) == len(components),
            "complete component cohort differs")
    flattened = []
    for component in components:
        members = sorted(row["id"] for row in variants if row["component"] == component["component"])
        require(members == sorted(component["variant_ids"])
                and len(members) == component["variant_count"], "component membership differs")
        flattened.extend(members)
    require(sorted(flattened) == ids, "component ledger does not partition every variant")
    return manifest, protocol, evaluation, variants, ledger


def plot(evaluation_path, out_path):
    root, out = Path(evaluation_path), Path(out_path)
    require(not out.exists(), "preserve existing figure output")
    manifest, protocol, result, variants, ledger = verified_inputs(root)
    macro = result["observed"]["equal_component"]["losses"]["paired"]
    mse = [macro[m]["mse"] for m in MODELS]
    primary = result["primary"]
    gain, interval = primary["observed_improvement"], primary["conditional_95_percent_range"]
    require(all(finite(x) and x >= 0 for x in mse) and finite(gain)
            and len(interval) == 2 and all(finite(x) for x in interval)
            and interval[0] <= interval[1]
            and primary["metric"] == "equal_component_macroMSE_zero_minus_I",
            "invalid plotted primary metrics")
    out.mkdir(parents=True)
    # Retain every normalized variant, model, cell and corresponding producer loss.
    scopes = ["paired", *CELLS]
    loss_columns = [f"{scope}_{model}_{metric}"
                    for scope in scopes for model in MODELS for metric in ["mse", "mae"]]
    with (out / "all-variants.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["id", "component", "locus", "gc_preserving",
                         *[f"{cell}_{key}" for cell in CELLS for key in ["observed", *MODELS]],
                         *loss_columns])
        for row in variants:
            values = [value for cell in CELLS for value in
                      [row["observed"][cell], 0.0, row["predictions"]["I"][cell],
                       row["predictions"]["S"][cell]]]
            losses = ledger[row["id"]]["losses"]
            writer.writerow([row["id"], row["component"], row["locus"], row["gc_preserving"],
                             *values, *[losses[scope][model][metric] for scope in scopes
                                        for model in MODELS for metric in ["mse", "mae"]]])
    with (out / "all-components.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["component", "variant_count", "loci_json", "variant_ids_json", *loss_columns])
        for row in result["per_component"]:
            writer.writerow([row["component"], row["variant_count"], json.dumps(row["loci"]),
                             json.dumps(row["variant_ids"]),
                             *[row["losses"][scope][model][metric] for scope in scopes
                               for model in MODELS for metric in ["mse", "mae"]]])
    with (out / "macro-model-errors.csv").open("w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(["model", "equal_component_equal_cell_MSE"])
        writer.writerows(zip(MODELS, mse))

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "helix-exploratory-variant-transfer-v2"})

    def save(fig, name):
        fig.savefig(out / (name + ".svg"), metadata={"Date": None,
                    "Creator": "Helix exploratory complete-case variant transfer"})
        fig.savefig(out / (name + ".png"), dpi=180,
                    metadata={"Software": "Matplotlib " + matplotlib.__version__})
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.3), gridspec_kw={"width_ratios": [1, 1.2]})
    fig.subplots_adjust(left=.08, right=.97, top=.72, bottom=.23, wspace=.38)
    axes[0].bar(range(3), mse, color=["#82909b", "#c47b33", "#377e89"])
    axes[0].set(xticks=range(3), xticklabels=["Zero effect", "Fixed I", "Fixed S"],
                ylabel="Macro MSE · published effect units²", title="Every fixed model retained")
    top = max(mse) or 1
    axes[0].set_ylim(0, top * 1.2)
    for index, value in enumerate(mse):
        axes[0].text(index, value + top * .035, f"{value:.6g}", ha="center", fontsize=9)
    axes[1].plot(interval, [0, 0], color="#c47b33", linewidth=4)
    axes[1].scatter([gain], [0], color="#222222", s=45, zorder=3)
    axes[1].axvline(0, color="#666666", linestyle="--", linewidth=1)
    low, high = min(0, gain, *interval), max(0, gain, *interval)
    padding = (high - low) * .12 or .1
    axes[1].set(xlim=(low - padding, high + padding), ylim=(-.6, .6), yticks=[],
                xlabel="Zero MSE − fixed I MSE (positive favors I)",
                title="Primary improvement · conditional 95% range")
    axes[1].text(.5, .85, f"{gain:.6g}  [{interval[0]:.6g}, {interval[1]:.6g}]",
                 transform=axes[1].transAxes, ha="center", fontsize=9)
    fig.text(.08, .94, TITLE, weight="bold", fontsize=14)
    fig.text(.08, .85,
             f"All {len(variants):,} complete-case variants · {result['component_count']:,} components · "
             "equal-component / equal-cell primary loss", fontsize=11)
    fig.text(.08, .08, "Fixed predictions; no refitting or calibration. Range conditions on fitted models "
             "and selected complete cases.\nNumerical evaluation does not resolve the failed v1 source-semantics gate.",
             fontsize=9, color="#444444", linespacing=1.6)
    save(fig, "01-model-transfer")

    values = [0.0, *[row["observed"][cell] for row in variants for cell in CELLS],
              *[row["predictions"][model][cell] for row in variants for model in ["I", "S"]
                for cell in CELLS]]
    low, high = min(values), max(values)
    padding = (high - low) * .05 or .1
    limits = (low - padding, high + padding)
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    fig.subplots_adjust(left=.09, right=.96, top=.81, bottom=.12, hspace=.3, wspace=.25)
    for row_index, cell in enumerate(CELLS):
        for column, model in enumerate(["I", "S"]):
            axis = axes[row_index, column]
            axis.scatter([r["observed"][cell] for r in variants],
                         [r["predictions"][model][cell] for r in variants],
                         s=4, alpha=.18, color="#c47b33" if model == "I" else "#377e89",
                         rasterized=True)
            axis.plot(limits, limits, color="#555555", linestyle="--", linewidth=.9)
            axis.axhline(0, color="#999999", linewidth=.7)
            axis.axvline(0, color="#cccccc", linewidth=.7)
            axis.set(xlim=limits, ylim=limits, xlabel="Observed alternate − reference effect",
                     ylabel="Fixed predicted alternate − reference effect",
                     title=f"{cell} · model {model} · all {len(variants):,} variants")
            axis.set_aspect("equal", adjustable="box")
    fig.text(.09, .95, "Exploratory numeric complete-case analysis", fontsize=15, weight="bold")
    fig.text(.09, .91, "After v1 source failure · both cells and fixed models · every complete-case variant",
             fontsize=10)
    fig.text(.09, .865, "Common full axes retain every observation and prediction; dashed diagonal denotes equality.",
             fontsize=10)
    fig.text(.09, .045, "No favorable filtering, clipping, refitting or fitted calibration. "
             "Units are the published signed reporter effect.\nThis exploratory complete-case evaluation "
             "does not establish the meaning of excluded source values.", fontsize=9,
             color="#444444", linespacing=1.6)
    save(fig, "02-all-variant-effects")
    outputs = {p.name: {"bytes": p.stat().st_size, "sha256": sha(p)} for p in sorted(out.iterdir())}
    lock = Path(__file__).resolve().parents[1] / "requirements/research.lock"
    receipt = {"schema_version": 1, "analysis_classification": CLASSIFICATION,
               "producer_manifest_sha256": sha(root / "manifest.json"),
               "producer_files": manifest["files"], "protocol_sha256": sha(root / "protocol.json"),
               "protocol_id": protocol["protocol_id"],
               "plot_script_sha256": sha(__file__), "dependency_lock_sha256": sha(lock),
               "matplotlib_version": matplotlib.__version__, "variants": len(variants),
               "components": result["component_count"], "scatter_axes_full_range": list(limits),
               "scatter_points_per_panel": len(variants),
               "filtering": "None after producer complete-case selection; every admitted variant and component retained.",
               "outputs": outputs}
    (out / "manifest.json").write_text(json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"analysis_classification": CLASSIFICATION, "variants": len(variants),
                      "components": result["component_count"], "outputs": len(outputs)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation", required=True)
    parser.add_argument("--out", required=True)
    arguments = parser.parse_args()
    plot(arguments.evaluation, arguments.out)
