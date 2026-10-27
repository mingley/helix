#!/usr/bin/env python3
"""Plot a fixed synthetic centering example; never reads empirical outcomes."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
INPUT_SHA256 = "52c765e96afe867d779d0b7706dee616903af410776f85545309c1c732417e2d"
PANEL_INPUT_SHA256 = "49d55f5b6bd306c305426612410ceec14f2c769a2ad2818a3cecea1544fd0ced"
OUTPUTS = (
    "identifiability-exact.csv",
    "identifiability.svg",
    "identifiability.png",
    "manifest.json",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def center(values: list[Fraction]) -> tuple[Fraction, list[Fraction]]:
    mean = sum(values, Fraction(0)) / len(values)
    return mean, [value - mean for value in values]


def plot_normalization_controlled(inputs: Path, out: Path) -> None:
    """Two distinct panels satisfy identical column moments before centering."""
    if sha256(inputs) != PANEL_INPUT_SHA256:
        raise SystemExit("Frozen normalization-controlled inputs do not match")
    data = json.loads(inputs.read_text())
    panels = {
        name: [[Fraction(value) for value in row] for row in data[key]]
        for name, key in (("X", "panel_x"), ("Y", "panel_y"))
    }
    centered = {name: [center(row)[1] for row in rows] for name, rows in panels.items()}
    if panels["X"] == panels["Y"] or centered["X"] != centered["Y"]:
        raise SystemExit("Synthetic distinct-panel centering identity failed")
    for rows in panels.values():
        if len(rows) != 4 or any(len(row) != 3 for row in rows):
            raise SystemExit("Expected a four-row, three-context panel")
        for column in zip(*rows, strict=True):
            if sum(column, Fraction(0)) != 0 or sum(value * value for value in column) != 14:
                raise SystemExit("Synthetic exact column moments failed")
    names = (
        "normalization-controlled-exact.csv",
        "normalization-controlled.svg",
        "normalization-controlled.png",
        "normalization-controlled-manifest.json",
    )
    out.mkdir(parents=True, exist_ok=True)
    if any((out / name).exists() for name in names):
        raise SystemExit("Output artifacts exist; choose a new --out directory")
    retained = out / "normalized-panel-inputs.json"
    if retained.resolve() != inputs.resolve():
        if retained.exists():
            raise SystemExit("Output normalized-panel-inputs.json already exists")
        retained.write_bytes(inputs.read_bytes())
    with (out / names[0]).open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["panel", "element", "context", "input_exact", "column_mean_exact", "column_sum_squares_exact", "column_sample_variance_exact", "input_row_mean_exact", "centered_input_exact", "normalized_center_exact"])
        for name, rows in panels.items():
            for i, row in enumerate(rows):
                mean, outputs = center(row)
                for context, value, output in zip(data["contexts"], row, outputs, strict=True):
                    expression = "0" if output == 0 else str(output) + "/sqrt(14/3)"
                    writer.writerow([name, data["rows"][i], context, str(value), "0", "14", "14/3", str(mean), str(output), expression])
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "svg.fonttype": "none", "svg.hashsalt": PANEL_INPUT_SHA256})
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 5.5))
    fig.subplots_adjust(left=0.10, right=0.98, top=0.74, bottom=0.27, wspace=0.25)
    scale = math.sqrt(14 / 3)
    matrices = [panels["X"], panels["Y"], [[float(v) / scale for v in row] for row in centered["X"]]]
    titles = ["Distinct synthetic panel X", "Distinct synthetic panel Y", "Same standardized, centered output"]
    for index, (ax, matrix, title) in enumerate(zip(axes, matrices, titles, strict=True)):
        bound = 3 if index < 2 else 1 / scale
        ax.imshow([[float(v) for v in row] for row in matrix], cmap="RdBu_r", vmin=-bound, vmax=bound, aspect="auto")
        ax.set_xticks(range(3), data["contexts"], fontsize=10)
        ax.set_yticks(range(4), data["rows"] if index == 0 else [""] * 4)
        ax.tick_params(length=0)
        ax.set_title(title, fontsize=11, pad=12)
        for i in range(4):
            for j in range(3):
                value = matrix[i][j]
                if index < 2:
                    label = str(value)
                else:
                    exact = centered["X"][i][j]
                    label = "0" if exact == 0 else ("1/s" if exact == 1 else "−1/s")
                color = "white" if abs(float(value)) > 0.5 * bound else "#17212B"
                ax.text(j, i, label, ha="center", va="center", color=color, fontsize=16)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("Known column moments still do not recover the common component", y=0.97, fontsize=15, fontweight="bold")
    fig.text(0.5, 0.88, "Synthetic full-panel counterexample — no empirical data or S9 reconstruction", ha="center", fontsize=11, color="#4A535D")
    fig.text(0.5, 0.17, r"Both panels: column means 0; column sums of squares 14; sample SD $s=\sqrt{14/3}$.", ha="center", fontsize=12)
    fig.text(0.5, 0.075, "The same per-column standardization followed by across-context row centering gives identical outputs.\nThe distinct common row activities remain unobserved; this does not identify the source's historical pipeline.", ha="center", fontsize=10, color="#4A535D")
    fig.savefig(out / names[1], metadata={"Date": None, "Creator": f"Matplotlib {matplotlib.__version__}; Helix synthetic illustration", "Description": data["kind"]})
    fig.savefig(out / names[2], dpi=200, metadata={"Software": f"Matplotlib {matplotlib.__version__}; Helix synthetic illustration"})
    plt.close(fig)
    manifest = {
        "schema_version": 1,
        "kind": data["kind"],
        "empirical_outcomes_consumed": False,
        "frozen_inputs_sha256": PANEL_INPUT_SHA256,
        "implementation": {"path": "scripts/plot_context_identifiability.py", "sha256": sha256(Path(__file__))},
        "formal_source": {"path": data["related_formal_source"], "sha256": sha256(ROOT / data["related_formal_source"]), "theorem": "HelixContextIdentifiability.no_exact_decoder_even_with_column_moment_constraints", "scope": "Ideal normalization-compatible synthetic panels, not the source's actual processing history. Source hash is not a proof-check receipt."},
        "runtime": {"python": platform.python_version(), "matplotlib": matplotlib.__version__},
        "dependency_lock_sha256": sha256(ROOT / "requirements/research.lock"),
        "exact_checks": {"panels_distinct": True, "column_sums": ["0", "0", "0"], "column_sums_squares": ["14", "14", "14"], "sample_variance": "14/3", "centered_panels_identical": True, "x_row_means": [str(center(row)[0]) for row in panels["X"]], "y_row_means": [str(center(row)[0]) for row in panels["Y"]]},
        "artifacts": {name: {"bytes": (out / name).stat().st_size, "sha256": sha256(out / name)} for name in ("normalized-panel-inputs.json", *names[:-1])},
        "interpretation_limits": data["not_claimed"],
    }
    (out / names[-1]).write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print("Wrote normalization-controlled exact CSV, SVG, PNG and manifest to", out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        type=Path,
        default=None,
    )
    parser.add_argument("--example", choices=["row-shift", "normalization-controlled"], default="row-shift")
    parser.add_argument(
        "--out", type=Path, default=ROOT / "research/context-identifiability-figures"
    )
    args = parser.parse_args()
    if args.example == "normalization-controlled":
        plot_normalization_controlled(args.inputs or ROOT / "research/context-identifiability-figures/normalized-panel-inputs.json", args.out)
        return
    args.inputs = args.inputs or ROOT / "research/context-identifiability-figures/inputs.json"
    if sha256(args.inputs) != INPUT_SHA256:
        raise SystemExit("Frozen synthetic inputs do not match the declared hash")
    data = json.loads(args.inputs.read_text())
    base = [Fraction(value) for value in data["base_profile"]]
    shift = Fraction(data["common_shift"])
    contexts = data["contexts"]
    if len(base) != 3 or len(contexts) != 3 or shift == 0:
        raise SystemExit("Expected exactly three contexts and a nonzero shift")
    shifted = [value + shift for value in base]
    base_mean, base_centered = center(base)
    shifted_mean, shifted_centered = center(shifted)
    if base_centered != shifted_centered or base == shifted:
        raise SystemExit("Synthetic exact centering identity failed")
    if sum(base_centered, Fraction(0)) != 0:
        raise SystemExit("Synthetic centered sum is not exactly zero")

    args.out.mkdir(parents=True, exist_ok=True)
    if any((args.out / name).exists() for name in OUTPUTS):
        raise SystemExit("Output artifacts already exist; choose a new --out directory")
    retained_input = args.out / "inputs.json"
    if retained_input.resolve() != args.inputs.resolve():
        if retained_input.exists():
            raise SystemExit("Output inputs.json already exists")
        retained_input.write_bytes(args.inputs.read_bytes())

    csv_path = args.out / OUTPUTS[0]
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "profile",
                "context",
                "input_exact",
                "common_shift_from_base_exact",
                "row_mean_exact",
                "centered_exact",
                "formal_center3_exact",
            ]
        )
        for name, values, mean, centered, common in (
            ("base", base, base_mean, base_centered, Fraction(0)),
            ("shifted", shifted, shifted_mean, shifted_centered, shift),
        ):
            for context, value, output in zip(contexts, values, centered, strict=True):
                writer.writerow(
                    [name, context, str(value), str(common), str(mean), str(output), str(3 * output)]
                )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "svg.fonttype": "none",
            "svg.hashsalt": INPUT_SHA256,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.8), gridspec_kw={"wspace": 0.32})
    fig.subplots_adjust(left=0.085, right=0.98, top=0.75, bottom=0.28)
    colors = ("#2368A2", "#CB6B20")
    x = list(range(3))
    width = 0.34
    axes[0].bar([v - width / 2 for v in x], list(map(float, base)), width, color=colors[0], label="Base profile")
    axes[0].bar([v + width / 2 for v in x], list(map(float, shifted)), width, color=colors[1], label="Base + common shift 5")
    axes[0].set_title("Distinct hypothetical inputs")
    axes[0].set_ylabel("Input to centering (arbitrary units)")
    axes[0].set_ylim(0, 10.5)
    axes[0].set_yticks([0, 2, 4, 6, 8, 10])
    axes[0].legend(frameon=False, loc="upper left", fontsize=9)
    for i, (low, high) in enumerate(zip(base, shifted, strict=True)):
        axes[0].text(i - width / 2, float(low) + 0.1, str(low), ha="center", fontsize=10)
        axes[0].text(i + width / 2, float(high) + 0.1, str(high), ha="center", fontsize=10)
    axes[1].plot(x, list(map(float, base_centered)), color=colors[0], linewidth=3, marker="o", markersize=11, label="Centered base")
    axes[1].plot(x, list(map(float, shifted_centered)), color=colors[1], linewidth=1.8, linestyle="--", marker="s", markersize=6, label="Centered shifted profile")
    axes[1].axhline(0, color="#68717A", linewidth=0.8)
    axes[1].set_title("Identical centered outputs")
    axes[1].set_ylabel("Centered value (same arbitrary units)")
    axes[1].set_ylim(-1.55, 1.65)
    axes[1].set_yticks([-1, 0, 1])
    axes[1].legend(frameon=False, loc="upper left", fontsize=9)
    for i, value in enumerate(base_centered):
        axes[1].text(i, float(value) - 0.27, str(value), ha="center", fontsize=11)
    for ax in axes:
        ax.set_xticks(x, contexts)
        ax.grid(axis="y", alpha=0.18)
        ax.set_axisbelow(True)
    fig.suptitle("Ideal across-context centering removes a common shift", x=0.5, y=0.97, fontsize=16, fontweight="bold")
    fig.text(0.5, 0.89, "Synthetic operator illustration — no empirical data or raw-activity reconstruction", ha="center", fontsize=11, color="#4A535D")
    fig.text(0.5, 0.16, r"$P(x)=x-\mathrm{mean}(x)\mathbf{1}$     $P(2,4,3)=P(7,9,8)=(-1,1,0)$", ha="center", fontsize=13)
    fig.text(0.5, 0.065, "Input-stage theorem only: upstream column standardization is not recomputed.\nThis does not assert which normalization produced the S9 workbook.", ha="center", fontsize=10, color="#4A535D")
    fig.savefig(
        args.out / "identifiability.svg",
        metadata={
            "Date": None,
            "Creator": f"Matplotlib {matplotlib.__version__}; Helix synthetic identifiability illustration",
            "Description": "Fixed synthetic centering-stage inputs only. Not S9 data, source-formula validation, or empirical raw-activity reconstruction.",
        },
    )
    fig.savefig(
        args.out / "identifiability.png",
        dpi=200,
        metadata={"Software": f"Matplotlib {matplotlib.__version__}; Helix synthetic illustration"},
    )
    plt.close(fig)

    proof = ROOT / data["related_formal_source"]
    manifest = {
        "schema_version": 1,
        "kind": data["kind"],
        "empirical_outcomes_consumed": False,
        "frozen_inputs_sha256": INPUT_SHA256,
        "implementation": {"path": "scripts/plot_context_identifiability.py", "sha256": sha256(Path(__file__))},
        "formal_source": {
            "path": data["related_formal_source"],
            "sha256": sha256(proof),
            "theorem": "HelixContextIdentifiability.common_shift_invariant",
            "scope": "Exact ideal centering-stage map; center3 is three times rational centering. A source hash is not itself a Lean verification receipt.",
        },
        "runtime": {"python": platform.python_version(), "matplotlib": matplotlib.__version__},
        "dependency_lock_sha256": sha256(ROOT / "requirements/research.lock"),
        "exact_checks": {
            "profiles_distinct": True,
            "centered_profiles_identical": True,
            "centered_sum_exactly_zero": True,
            "base_mean": str(base_mean),
            "shifted_mean": str(shifted_mean),
            "centered_values": list(map(str, base_centered)),
            "formal_center3_values": list(map(lambda x: str(3 * x), base_centered)),
        },
        "artifacts": {
            name: {"bytes": (args.out / name).stat().st_size, "sha256": sha256(args.out / name)}
            for name in ("inputs.json", *OUTPUTS[:-1])
        },
        "interpretation_limits": data["not_claimed"] + [
            "The invariance is at the input to centering, not arbitrary raw shifts with a recomputed upstream normalizer.",
            "Deterministic rendering is expected with the pinned runtime, font and input bytes; alternative rendering stacks may differ.",
        ],
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print("Wrote exact CSV, SVG, PNG and provenance manifest to", args.out)


if __name__ == "__main__":
    main()
