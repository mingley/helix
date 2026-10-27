#!/usr/bin/env python3
"""Standalone ESS figures and exact numeric tables; never relabel ESS as activity."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(source, output, arithmetic=None):
    root, out = Path(source).resolve(), Path(output)
    if out.exists():
        raise ValueError("preserve previous figures; output exists")
    manifest = json.loads((root / "artifacts.json").read_bytes())
    if manifest["schema_version"] != 4:
        raise ValueError("expected specificity-study schema4")
    for name, expected in manifest["files"].items():
        path = (root / name).resolve()
        if root not in path.parents or digest(path) != expected:
            raise ValueError("input path/checksum mismatch: " + name)
    load = lambda name: json.loads((root / name).read_bytes())
    summary, grouping = load("summary.json"), load("grouping.json")
    with gzip.open(root / "study.json.gz", "rt") as handle:
        study = json.load(handle)
    with gzip.open(root / "resampling.json.gz", "rt") as handle:
        resampling = json.load(handle)
    out.mkdir(parents=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "helix-context-specificity-v2"})
    files = []

    def save(fig, name):
        fig.savefig(out / (name + ".svg"), bbox_inches="tight",
                    metadata={"Date": None, "Creator": "Helix context-specificity study"})
        fig.savefig(out / (name + ".png"), bbox_inches="tight", dpi=180,
                    metadata={"Software": "Helix; Matplotlib " + matplotlib.__version__})
        plt.close(fig)
        files.extend([name + ".svg", name + ".png"])

    labels = ["G: global mean", "C: cell mean", "S: shared DNA", "A: additive", "I: cell × DNA"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    for axis, field, title in zip(axes, ["sequence_weighted_mse", "macro_group_mse"],
                                 ["Primary: equal elements", "Sensitivity: equal components"]):
        values = [row[field] for row in summary["models"]]
        axis.bar(range(5), values, color=["#8795a4", "#a2b6ca", "#8795a4", "#4c78a8", "#e99431"])
        axis.set_xticks(range(5), labels, rotation=25, ha="right")
        axis.set(title=title, ylabel="MSE (published ESS)²", ylim=(0, max(values) * 1.22))
        for j, value in enumerate(values):
            axis.text(j, value + max(values)*0.018, f"{value:.5f}", ha="center", fontsize=8)
    axes[0].axhline(summary["fixed_zero_prediction_diagnostic"]["mse"], color="#444444", ls="--", lw=1,
                    label="fixed zero diagnostic")
    axes[0].legend(fontsize=8)
    fig.suptitle("Held-out families · published score columns; normalization semantics unresolved")
    save(fig, "01-model-errors")
    with (out / "model-errors.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["model", "equal_element_mse", "equal_element_mae", "equal_component_mse"])
        writer.writerows((r["model"], r["sequence_weighted_mse"], r["sequence_weighted_mae"], r["macro_group_mse"])
                         for r in summary["models"])
    files.append("model-errors.csv")

    pairs = {}
    for row in study["context_contrast_group_errors"]:
        key = row["context_left"], row["context_right"]
        item = pairs.setdefault(key, {"n": 0, "loss_sums": np.zeros(5)})
        item["n"] += row["sequences"]
        item["loss_sums"] += row["sequences"] * np.array(row["model_losses"]["mse"])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    keys = sorted(pairs)
    offsets = np.arange(len(keys))
    for model, shift, color, label in [(3, -0.18, "#4c78a8", "A: shared slopes + cell offset"),
                                        (4, 0.18, "#e99431", "I: cell-specific slopes")]:
        axes[0].bar(offsets + shift, [pairs[k]["loss_sums"][model]/pairs[k]["n"] for k in keys],
                    width=0.34, color=color, label=label)
    axes[0].set_xticks(offsets, [f"{a} ↔ {b}" for a, b in keys], rotation=20)
    axes[0].set(ylabel="Equal-element pair-contrast MSE", title="Every cell pair retained")
    axes[0].legend(fontsize=8)
    for j, row in enumerate(resampling["comparisons"]):
        lo, hi = row["mse_interval"]
        axes[1].plot([lo, hi], [j, j], lw=3, color="#4c78a8")
        axes[1].scatter([row["observed_mse_improvement"]], [j], color="#333333", zorder=3)
    axes[1].axvline(0, color="#444444", ls="--", lw=1)
    axes[1].set_yticks(range(3), ["Context offset (diagnostic)", "Shared DNA (diagnostic)", "Cell × DNA (primary)"])
    axes[1].set(xlabel="Baseline MSE − candidate MSE", title="2,000 conditional component resamples")
    axes[1].invert_yaxis()
    fig.suptitle("Contrast prediction and paired uncertainty · contrast identity is not independent validation")
    save(fig, "02-context-contrasts")
    with (out / "context-contrasts.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["left", "right", "elements", "additive_mse", "interaction_mse"])
        writer.writerows((a, b, pairs[a,b]["n"], pairs[a,b]["loss_sums"][3]/pairs[a,b]["n"],
                          pairs[a,b]["loss_sums"][4]/pairs[a,b]["n"]) for a, b in keys)
    files.append("context-contrasts.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    folds = grouping["folds"]
    bottoms = np.zeros(len(folds))
    for field, color, label in [("train_ids", "#4c78a8", "Training"), ("test_ids", "#e99431", "Test"),
                                ("guard_ids", "#b6bec8", "Guard exclusions")]:
        values = [len(f[field]) for f in folds]
        axes[0].bar(range(len(folds)), values, bottom=bottoms, color=color, label=label)
        bottoms += values
    axes[0].set(xticks=range(5), xticklabels=[f"Fold {j+1}" for j in range(5)], ylabel="Assayed elements",
                title="All contexts travel with each element")
    axes[0].legend(fontsize=8)
    sizes = [component["size"] for component in grouping["components"]]
    axes[1].hist(sizes, bins=30, color="#4c78a8")
    axes[1].set(xlabel="Elements per connected component", ylabel="Components", title="Spatial and sequence-family support")
    fig.suptitle(f"{summary['sequence_count']:,} reference-validated cores · {summary['group_count']:,} indivisible components")
    save(fig, "03-fold-support")

    fig, axis = plt.subplots(figsize=(11, 5.2), constrained_layout=True)
    axis.axis("off")
    axis.text(0.02, 0.93, "Stated ESS formula — and a failed source check", fontsize=17, weight="bold")
    axis.text(0.02, 0.75, "Published activity → full-panel per-cell z-score → subtract each element's across-cell mean", fontsize=11)
    for y, source_values, result, label in [(0.53, "[+1, +1, +1]", "[0, 0, 0]", "Uniformly high standardized activity"),
                                            (0.37, "[−1, −1, −1]", "[0, 0, 0]", "Uniformly low standardized activity"),
                                            (0.21, "[+1, 0, −1]", "[+1, 0, −1]", "Relative cell preference remains")]:
        axis.text(0.03, y, source_values + "    →    " + result, fontsize=14, family="DejaVu Sans Mono")
        axis.text(0.60, y, label, fontsize=11)
    maximum_sum=summary["fixed_zero_prediction_diagnostic"]["maximum_source_element_sum_absolute"]
    axis.text(0.02, -0.02, "Illustrative formula inputs, not measured rows. Under this formula, three scores sum to zero.\n"
              f"Frozen workbook: maximum |sum| = {maximum_sum:.6f}, contradicting unweighted centering.\n"
              "Source values are unchanged; semantics remain unresolved. Zero under the stated formula does not mean inactive.", fontsize=10, color="#444444")
    save(fig, "04-specificity-estimand")
    arithmetic_manifest_sha=None
    if arithmetic is not None:
        audit_root=Path(arithmetic).resolve()
        arithmetic_manifest=json.loads((audit_root/"artifacts.json").read_bytes())
        for name,info in arithmetic_manifest.items():
            path=(audit_root/name).resolve()
            if audit_root not in path.parents or digest(path)!=info["sha256"]:
                raise ValueError("arithmetic artifact path/checksum mismatch: "+name)
        audit=json.loads((audit_root/"summary.json").read_bytes())
        inputs=load("inputs.json")
        if (audit["hashes"]["panel_sha256"]!=inputs["panel_sha256"]
                or audit["hashes"]["parent_protocol_sha256"]!=inputs["protocol_sha256"]):
            raise ValueError("arithmetic diagnostic is not bound to these source inputs")
        arithmetic_manifest_sha=digest(audit_root/"artifacts.json")
        cohort=audit["all55338finite_exact_source_score_IDs"]
        raw=cohort["raw_svd"]
        contexts=list(raw["weights_by_context"])
        weights=[raw["weights_by_context"][c] for c in contexts]
        fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
        axes[0].bar(contexts,weights,color=["#4c78a8","#e99431","#54a24b"])
        axes[0].axhline(1,color="#444444",ls="--",lw=1,label="unweighted formula")
        axes[0].set(ylabel="Null-vector weight (normalized mean = 1)",ylim=(0.8,1.15),title="A weighted score plane")
        for j,weight in enumerate(weights):
            axes[0].text(j,weight+0.008,f"{weight:.6f}",ha="center",fontsize=9)
        axes[0].legend(fontsize=8)
        residuals=[cohort["row_sum"]["maximum_absolute"],raw["raw_matrix_weighted_residual"]["maximum_absolute"]]
        axes[1].bar(["Unweighted sum","Fitted weighted sum"],residuals,color=["#e99431","#4c78a8"])
        axes[1].set(yscale="log",ylabel="Maximum absolute residual",title=f"All {cohort['element_count']:,} source score rows")
        for j,value in enumerate(residuals):
            axes[1].text(j,value*1.8,f"{value:.3g}",ha="center",fontsize=9)
        fig.suptitle("Post-hoc source arithmetic: compatible with extra scaling, not proof of producer history")
        fig.supxlabel("Audit weights were not applied to targets, fits or evaluation. Source semantics remain unresolved.",fontsize=9)
        save(fig,"05-source-score-plane")
    metadata = {"schema_version": 1, "input_study_manifest_sha256": digest(root / "artifacts.json"),
                "script_sha256": digest(Path(__file__)),
                "arithmetic_artifacts_sha256":arithmetic_manifest_sha,
                "software": {"matplotlib": matplotlib.__version__, "numpy": np.__version__},
                "files": {name: digest(out / name) for name in files},
                "scope": "Published ESS reanalysis. Whole-source-panel preprocessing includes future held-out elements; no raw activity or new biological mechanism claim."}
    (out / "figures.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print("Saved standalone SVG/PNG figures, exact CSV tables and source commitments to", out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--arithmetic")
    args = parser.parse_args()
    run(args.run, args.out,args.arithmetic)
