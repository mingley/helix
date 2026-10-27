#!/usr/bin/env python3
"""Standalone figures for the frozen published-activity endpoint; no fitting."""
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


MODELS = ["global_mean", "context_mean", "shared_dna", "additive_context_dna",
          "independent_context_dna"]
LABELS = ["G: global mean", "C: cell mean", "S: shared DNA", "A: additive", "I: cell × DNA"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(source, out):
    source, out = Path(source).resolve(), Path(out)
    if out.exists():
        raise FileExistsError("preserve existing figures")
    artifacts = json.loads((source / "artifacts.json").read_bytes())
    if artifacts["schema_version"] != 5:
        raise ValueError("expected the separate activity study schema5")
    for name, expected in artifacts["files"].items():
        path = (source / name).resolve()
        if source not in path.parents or sha(path) != expected:
            raise ValueError("input artifact path/hash differs")
    with gzip.open(source / "study.json.gz", "rt") as stream:
        study = json.load(stream)
    with gzip.open(source / "resampling.json.gz", "rt") as stream:
        bootstrap = json.load(stream)
    out.mkdir(parents=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "helix-context-activity-v3"})

    def save(fig, name):
        fig.savefig(out / (name + ".svg"), bbox_inches="tight",
                    metadata={"Date": None, "Creator": "Helix published-activity study"})
        fig.savefig(out / (name + ".png"), bbox_inches="tight", dpi=180,
                    metadata={"Software": "Matplotlib " + matplotlib.__version__})
        plt.close(fig)

    models = study["models"]
    values = [item["sequence_weighted_mse"] for item in models]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    colors = ["#8795a4", "#789bbd", "#8795a4", "#d6872c", "#9265a8"]
    axes[0].bar(range(5), values, color=colors)
    axes[0].set_xticks(range(5), LABELS, rotation=25, ha="right")
    axes[0].set(ylabel="MSE · published log₂(RNA/DNA)²", title="Every fixed model retained",
                ylim=(0, max(values) * 1.15))
    for index, value in enumerate(values):
        axes[0].text(index, value + max(values)*.02, f"{value:.5f}", ha="center", fontsize=8)
    for row, comparison_index, label, color in [(0, 1, "Primary: C → A", "#d6872c"),
                                               (1, 2, "Secondary: A → I", "#9265a8")]:
        comparison = bootstrap["comparisons"][comparison_index]
        axes[1].plot(comparison["mse_interval"], [row, row], lw=4, color=color)
        axes[1].scatter([comparison["observed_mse_improvement"]], [row], color="#222222", zorder=3)
    axes[1].axvline(0, color="#444444", linestyle="--", lw=1)
    axes[1].set(yticks=[0, 1], yticklabels=["Primary: C → A", "Secondary: A → I"],
                xlabel="Baseline MSE − candidate MSE", title="Conditional 95% component resamples",
                ylim=(-.5, 1.5))
    axes[1].invert_yaxis()
    fig.suptitle("Grouped held-out reporter activity · same selected panel, three measured cells")
    save(fig, "01-model-activity")
    with (out / "model-errors.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["model", "equal_element_mse", "equal_element_mae", "equal_component_mse"])
        writer.writerows((r["model"], r["sequence_weighted_mse"], r["sequence_weighted_mae"],
                          r["macro_group_mse"]) for r in models)

    grouped = {}
    for prediction in study["predictions"]:
        grouped.setdefault(prediction["sequence_id"], []).append(prediction)
    rows = []
    for identity, predictions in sorted(grouped.items()):
        if len(predictions) != 3:
            raise ValueError("activity plot requires complete paired observations")
        rows.append((identity, predictions[0]["group"],
                     float(np.mean([r["observed"] for r in predictions])),
                     float(np.mean([r["predictions"]["global_mean"] for r in predictions])),
                     float(np.mean([r["predictions"]["shared_dna"] for r in predictions]))))
    with (out / "common-activity.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["id", "component", "observed_mean_activity", "training_mean_prediction",
                         "shared_DNA_prediction"])
        writer.writerows(rows)
    array = np.asarray([row[2:] for row in rows])
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), constrained_layout=True)
    for axis, column, title, color in [(axes[0], 1, "Training mean", "#789bbd"),
                                       (axes[1], 2, "Shared DNA · whole held-out groups", "#d6872c")]:
        axis.scatter(array[:, 0], array[:, column], s=3, alpha=.12, color=color, rasterized=True)
        lo = min(array[:, 0].min(), array[:, column].min())
        hi = max(array[:, 0].max(), array[:, column].max())
        axis.plot([lo, hi], [lo, hi], color="#555555", lw=.8, linestyle="--")
        axis.set(xlabel="Observed mean across three published activity scales",
                 ylabel="Held-out predicted mean activity", title=title)
    fig.suptitle("A measured common component · library-normalized activity, not absolute expression")
    save(fig, "02-common-activity")

    pairs = {}
    for row in study["context_contrast_group_errors"]:
        key = (row["context_left"], row["context_right"])
        item = pairs.setdefault(key, {"n": 0, "loss": np.zeros(5)})
        item["n"] += row["sequences"]
        item["loss"] += row["sequences"] * np.asarray(row["model_losses"]["mse"])
    keys = sorted(pairs)
    fig, axis = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    positions = np.arange(len(keys))
    for model, offset, label, color in [(3, -.18, "A: shared slopes + cell offsets", "#d6872c"),
                                      (4, .18, "I: cell-specific slopes", "#9265a8")]:
        axis.bar(positions+offset, [pairs[k]["loss"][model]/pairs[k]["n"] for k in keys],
                 width=.34, label=label, color=color)
    axis.set_xticks(positions, [f"{a} ↔ {b}" for a, b in keys])
    axis.set(ylabel="Pair-contrast MSE · log₂(RNA/DNA)²",
             title="Secondary: cell contrasts · correlated diagnostics")
    axis.legend(fontsize=8)
    save(fig, "03-cell-contrasts")

    manifest = {"schema_version": 1, "endpoint": "published_joint_log2_rna_dna_activity",
                "source_artifacts_sha256": sha(source / "artifacts.json"),
                "implementation_sha256": sha(__file__), "matplotlib_version": matplotlib.__version__,
                "fits_performed": False, "all_source_elements_plotted": len(rows),
                "files": {p.name: {"bytes": p.stat().st_size, "sha256": sha(p)}
                          for p in sorted(out.iterdir())}}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"figures": 3, "elements": len(rows), "output": str(out)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run(args.source, args.out)
