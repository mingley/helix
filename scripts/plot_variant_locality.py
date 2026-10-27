#!/usr/bin/env python3
"""Plot the completed outcome-free panel scan; no model or biological fitting."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def plot(source, out):
    source, out = Path(source).resolve(), Path(out)
    if out.exists():
        raise FileExistsError("preserve existing figure outputs")
    manifest = json.loads((source / "manifest.json").read_text())
    for name, commitment in manifest["files"].items():
        path = (source / name).resolve()
        if source not in path.parents or path.is_symlink():
            raise ValueError("input artifact path differs")
        if path.stat().st_size != commitment["bytes"] or sha(path) != commitment["sha256"]:
            raise ValueError("input artifact commitment differs: " + name)
    summary = json.loads((source / "summary.json").read_text())
    histogram, identities = collections.Counter(), set()
    with (source / "per-element.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            identity = row["element_id"]
            count = int(row["exact_zero_delta_interior_SNVs"])
            opportunities = int(row["all_interior_SNVs"])
            if identity in identities or opportunities != 588 or not 0 <= count <= 196:
                raise ValueError("invalid or repeated panel record")
            if abs(float(row["exact_zero_fraction"]) - count / 588) > 1e-15:
                raise ValueError("per-element fraction differs")
            identities.add(identity)
            histogram[count] += 1
    elements = sum(histogram.values())
    zero_opportunities = sum(k * n for k, n in histogram.items())
    total_opportunities = elements * 588
    recorded_histogram = {int(k): v for k, v in summary["per_element_exact_zero_count"]["histogram"].items()}
    if (elements != 55_276 or elements != summary["panel_slots"]
            or histogram != recorded_histogram
            or zero_opportunities != summary["exact_zero_feature_delta_opportunities"]
            or total_opportunities != summary["interior_SNV_opportunities"]):
        raise ValueError("complete panel counts differ from the immutable scan")
    out.mkdir(parents=True)
    bins = list(range(max(histogram) + 1))
    with (out / "feature-equivalent-counts.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["exact_zero_feature_delta_SNVs_per_element", "panel_elements"])
        writer.writerows((k, histogram[k]) for k in bins)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "helix-variant-locality-panel-scan-v1"})
    fig, axis = plt.subplots(figsize=(10, 5.8))
    fig.subplots_adjust(left=.09, right=.98, top=.76, bottom=.23)
    axis.bar(bins, [histogram[k] for k in bins], color="#377e89", width=.86)
    median = summary["per_element_exact_zero_count"]["median"]
    axis.axvline(median, color="#b45923", linestyle="--", linewidth=1.6,
                 label=f"Median: {median:g}")
    axis.set(xlabel="Feature-equivalent interior SNVs per 200nt element (out of 588)",
             ylabel="Panel elements", xlim=(-.75, max(bins) + .75),
             ylim=(0, max(histogram.values()) * 1.14))
    axis.set_xticks(range(0, max(bins) + 1, 3))
    axis.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    axis.legend(frameon=False, loc="upper right")
    fig.text(.09, .94, "Exact feature equivalence in the fixed human DNA panel",
             fontsize=16, weight="bold", ha="left")
    fig.text(.09, .875,
             f"{zero_opportunities:,} / {total_opportunities:,} proposed interior SNVs "
             f"({100 * zero_opportunities / total_opportunities:.4f}%) leave every feature unchanged.",
             fontsize=11, ha="left")
    fig.text(.09, .815, f"All {elements:,} elements retained; no response measurements or outcome filtering.",
             fontsize=10, ha="left", color="#444444")
    fig.text(.09, .10,
             "Features: 32 reverse-complement-canonical 3mer frequencies + GC fraction.\n"
             "Positions 3–198 (one-based), three alternate bases per position; panel slots counted as deposited.\n"
             "Exact feature equality is a representation limit, not evidence of biological neutrality.",
             fontsize=9, linespacing=1.5, ha="left", va="center", color="#444444")
    fig.savefig(out / "panel-feature-equivalence.svg", metadata={
        "Date": None, "Creator": "Helix outcome-free variant-locality scan"})
    fig.savefig(out / "panel-feature-equivalence.png", dpi=180,
                metadata={"Software": "Matplotlib " + matplotlib.__version__})
    plt.close(fig)
    inputs = {name: {"bytes": (source / name).stat().st_size,
                     "sha256": sha(source / name)}
              for name in ["manifest.json", "summary.json", "per-element.csv"]}
    outputs = {path.name: {"bytes": path.stat().st_size, "sha256": sha(path)}
               for path in sorted(out.iterdir())}
    lock = Path(__file__).resolve().parents[1] / "requirements/research.lock"
    receipt = {
        "schema_version": 1,
        "purpose": "Descriptive histogram of the completed outcome-free exact DNA scan.",
        "script": {"path": "scripts/plot_variant_locality.py", "sha256": sha(__file__)},
        "dependency_lock": {"path": "requirements/research.lock", "sha256": sha(lock)},
        "matplotlib_version": matplotlib.__version__,
        "inputs": inputs,
        "complete_panel_elements": elements,
        "exact_zero_feature_delta_opportunities": zero_opportunities,
        "all_interior_SNV_opportunities": total_opportunities,
        "filtering": "None; every original panel row is included.",
        "biological_interpretation": "Feature equivalence does not imply biological neutrality.",
        "outputs": outputs,
    }
    (out / "manifest.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"elements": elements, "zero_opportunities": zero_opportunities,
                      "total_opportunities": total_opportunities, "outputs": len(outputs)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default="research/variant-locality/panel-scan-v1")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    plot(args.source, args.out)
