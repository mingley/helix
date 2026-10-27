#!/usr/bin/env python3
"""Create standalone research figures and exact data tables from an audited study."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="helix-mpl-"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ORDER = ["mean", "chemistry", "coordinate", "delta3mer", "motif", "full"]
COLORS = ["#6b7785", "#4c78a8", "#9ecae9", "#f58518", "#54a24b", "#b279a2"]
TITLES = {
    "hbb_guarded_position_cv": "HBB · guarded positions",
    "primary_leave_element_out": "Nine promoters · held-out elements",
    "hbb_to_external_transfer": "Train HBB · test other promoters",
    "row_random_diagnostic": "HBB · row-random diagnostic",
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def run(directory, output):
    root, out = Path(directory).resolve(), Path(output)
    if out.exists():
        raise ValueError("figure output already exists; preserve previous artifacts")
    manifest_bytes = (root / "artifacts.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest["schema_version"] != 3:
        raise ValueError("expected integrated study artifact version3")
    for name, expected in manifest["files"].items():
        path = (root / name).resolve()
        if root not in path.parents or digest(path.read_bytes()) != expected:
            raise ValueError("input artifact path/checksum mismatch: " + name)
    load = lambda name: json.loads((root / name).read_text())
    summary, contexts = load("summary.json"), load("context.json")
    controls, forensic = load("resampling.json"), load("forensics.json")
    with gzip.open(root / "study.json.gz", "rt") as handle:
        full = json.load(handle)
    out.mkdir(parents=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
        "axes.spines.right": False, "svg.hashsalt": "helix-regulatory-study-v1",
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    created = []

    def save(figure, stem):
        figure.savefig(out / (stem + ".svg"), bbox_inches="tight",
                       metadata={"Date": None, "Creator": "Helix regulatory study"})
        figure.savefig(out / (stem + ".png"), dpi=180, bbox_inches="tight",
                       metadata={"Software": "Helix; Matplotlib " + matplotlib.__version__})
        plt.close(figure)
        created.extend([stem + ".svg", stem + ".png"])

    families = {row["family"]: row for row in summary["analyses"]}
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    csv_rows = []
    for axis, family in zip(axes.flat, TITLES):
        models = {row["model"]: row for row in families[family]["models"]}
        values = [models[name]["macro_assay_position_mse"] for name in ORDER]
        axis.bar(ORDER, values, color=COLORS)
        axis.axhline(values[0], color="#333333", ls="--", lw=1, label="training-mean baseline")
        axis.set(title=TITLES[family], ylabel="Macro position MSE (log2 effect)²", ylim=(0, max(values) * 1.22))
        axis.tick_params(axis="x", rotation=25)
        for index, value in enumerate(values):
            axis.text(index, value + max(values)*0.012, f"{value:.4f}", ha="center", va="bottom", fontsize=8)
        for name in ORDER:
            csv_rows.append([family, name, models[name]["macro_assay_position_mse"],
                             models[name]["macro_assay_position_mae"]])
    fig.suptitle("Fixed regulatory models: primary errors, with all baselines retained", fontsize=14)
    save(fig, "01-model-errors")
    with (out / "model-errors.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["family", "model", "macro_position_mse", "macro_position_mae"])
        writer.writerows(csv_rows)
    created.append("model-errors.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    contrast_rows = []
    for axis, name, title in zip(axes, ["hbb_position_block", "primary_element"],
                                 ["HBB · 21bp block resampling", "Nine promoters · element resampling"]):
        comparisons = controls["bootstrap"][name]["comparisons"]
        for index, row in enumerate(comparisons):
            center = row["observed_mean_loss_improvement"]
            lo, hi = row["percentile_95_lower"], row["percentile_95_upper"]
            axis.plot([lo, hi], [index, index], color="#4c78a8", lw=3)
            axis.scatter([center], [index], color="#222222", zorder=3)
            contrast_rows.append([name, row["comparison_id"], center, lo, hi])
        axis.axvline(0, color="#555555", ls="--", lw=1)
        axis.set_yticks(range(len(comparisons)), [r["comparison_id"].replace("_", " ") for r in comparisons])
        axis.set(title=title, xlabel="Baseline MSE − candidate MSE\nPositive favors candidate")
        axis.invert_yaxis()
    fig.suptitle("Prespecified contrasts · conditional 95% resampling ranges", fontsize=13)
    save(fig, "02-planned-contrasts")
    with (out / "planned-contrasts.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["family", "comparison", "observed_improvement", "conditional_95_lower", "conditional_95_upper"])
        writer.writerows(contrast_rows)
    created.append("planned-contrasts.csv")

    fig, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)
    context_rows = []
    for axis, pair in zip(axes.flat, contexts["pairs"]):
        left = np.array([row["left_coefficient"] for row in pair["matches"]])
        right = np.array([row["right_coefficient"] for row in pair["matches"]])
        axis.scatter(right, left, s=8, alpha=0.3, color="#4c78a8", linewidths=0)
        lower, upper = min(left.min(), right.min()), max(left.max(), right.max())
        margin = max((upper - lower) * 0.04, 0.01)
        axis.plot([lower - margin, upper + margin], [lower - margin, upper + margin], color="#555555", ls="--", lw=1)
        axis.set(xlim=(lower - margin, upper + margin), ylim=(lower - margin, upper + margin),
                 title=pair["relationship"].replace("_", " "), xlabel="Right source log2 coefficient", ylabel="Left source log2 coefficient")
        axis.text(0.04, 0.97, f"n={len(left)}; Pearson={pair['pearson']:.3f}\nMSE={pair['mse']:.4f}",
                  transform=axis.transAxes, va="top", fontsize=9)
        context_rows.append([pair["relationship"], pair["left"], pair["right"], len(left),
                             pair["mse"], pair["mae"], pair["pearson"],
                             pair["finite_context_blind_mse_lower_bound"]])
    axes.flat[-1].axis("off")
    axes.flat[-1].text(0, 0.9,
        "Exact genomic allele + identical target sequence\n\n"
        "Each comparison is separate.\nDifferences include estimator noise and context.\n"
        "A library repeat in another gene is not\nan element-specific causal noise baseline.\n\n"
        "The conditional same-input label bound\nis pair MSE / 4, not biological irreducibility.",
        va="top", fontsize=11)
    fig.suptitle("Matched source estimates across repeated libraries and contexts", fontsize=14)
    save(fig, "03-context-comparisons")
    with (out / "context-comparisons.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["relationship", "left", "right", "matched_variants", "mse", "mae", "pearson", "finite_label_mse_bound"])
        writer.writerows(context_rows)
    created.append("context-comparisons.csv")

    hbb = next(a for a in full["analyses"] if a["family"] == "hbb_guarded_position_cv")
    models = {m["model"]: m for m in hbb["models"]}
    rows = models["full"]["predictions"]
    by_position = {}
    for row in rows:
        by_position.setdefault(row["position_1based"], []).append(row)
    positions = sorted(by_position)
    observed = [np.mean([x["observed"] for x in by_position[p]]) for p in positions]
    predicted = [np.mean([x["predicted"] for x in by_position[p]]) for p in positions]
    fig, axis = plt.subplots(figsize=(12, 4), constrained_layout=True)
    axis.axvspan(0.5, 50.5, color="#bcd6f0", alpha=0.5, label="Native HBB 5′UTR annotation")
    axis.scatter([x["position_1based"] for x in rows], [x["observed"] for x in rows],
                 color="#777777", s=10, alpha=0.3, label="Source SNV coefficient")
    axis.plot(positions, observed, color="#222222", lw=1.2, label="Observed position mean")
    axis.plot(positions, predicted, color="#b279a2", lw=1.4, label="Full-model guarded OOF position mean")
    axis.axhline(0, color="#777777", lw=0.6)
    axis.set(xlim=(0.5, 187.5), xlabel="Submitted target position (plus genomic orientation)",
             ylabel="Log2 reporter-effect coefficient",
             title="HBB landscape · observed effects versus fixed-model guarded predictions")
    axis.legend(fontsize=8, loc="upper right")
    save(fig, "04-hbb-landscape")

    reconciliation = forensic["reconciliation"]
    inconsistency = reconciliation["source_estimate_inconsistency"]
    qc = inconsistency["predefined_minimum_ten_tag_qc"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    counts = [reconciliation["equality_to_original_deletion_matches"],
              reconciliation["distinct_deletion_altered_target_sequences"],
              qc["retained_deletion_rows"], qc["distinct_altered_target_sequences"]]
    axes[0].bar(["Raw rows", "Raw sequence\nidentities", "QC rows", "QC sequence\nidentities"], counts,
                color=["#f58518", "#f2b277", "#54a24b", "#9bc997"])
    for i, value in enumerate(counts): axes[0].text(i, value, str(value), ha="center", va="bottom")
    axes[0].set(ylabel="Count", title="Source deletions and whole-edited-sequence identity", ylim=(0, 70))
    raw_floor = float(inconsistency["lower_bound_mse_over_all_deletion_rows"])
    qc_floor = float(qc["lower_bound_mse_from_repeated_identity"])
    axes[1].bar(["Raw 62 rows", "QC 35 rows"], [raw_floor, qc_floor], color=["#f58518", "#54a24b"])
    axes[1].set(ylabel="Finite observed-label MSE bound", title="Duplicate-label constraint disappears after QC",
                ylim=(0, raw_floor * 1.25))
    axes[1].text(0, raw_floor, f"{raw_floor:.6f}", ha="center", va="bottom")
    axes[1].text(1, 0, "0", ha="center", va="bottom")
    fig.suptitle("HBB deletion-source forensics · not a biological noise floor", fontsize=14)
    save(fig, "05-source-identity-qc")

    files = {name: digest((out / name).read_bytes()) for name in created}
    (out / "figures.json").write_text(json.dumps({
        "schema_version": 1, "input_study_manifest_sha256": digest(manifest_bytes),
        "script_sha256": digest(Path(__file__).read_bytes()),
        "software": {"python": sys.version.split()[0], "numpy": np.__version__, "matplotlib": matplotlib.__version__},
        "files": files,
        "interpretation": "Standalone figures from committed source estimates and frozen analyses; no calibrated significance, new biological mechanism or causal context assignment is asserted."
    }, indent=2) + "\n")
    print("Saved five SVG/PNG research figures and exact numeric CSV tables to", out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run(args.run, args.out)
