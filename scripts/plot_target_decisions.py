"""Draw Figure 5 from the archived GAT and random-forest study summaries."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter


FORWARD = "#316985"
REVERSE = "#A5714B"
INK = "#243642"
SUBTLE = "#60717C"
GUIDE = "#C2CDD2"
SUMMARY = Path("results/p0_rechecked_20261117/primary_common_draw_summary.csv")
SOURCES = Path("manuscript/source_data")
TARGETS = (
    ("HMP2_to_PRISM", "pooled", "visit", "Franzosa pooled", 220, 220),
    ("HMP2_to_PRISM", "discovery", "visit", "PRISM discovery", 155, 155),
    ("HMP2_to_PRISM", "validation", "visit", "Dutch validation", 65, 65),
    ("PRISM_to_HMP2", "pooled", "visit", "HMP2 records", 1626, 130),
    ("PRISM_to_HMP2", "pooled", "participant", "HMP2 participants", 130, 130),
)
PAIR_KEYS = (("discovery", "discovery"), ("discovery", "validation"),
             ("validation", "discovery"), ("validation", "validation"))
BRIER_STAGES = ("record_weighted_loss", "participant_equal_record_loss",
                "participant_mean_probability")


def records(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def unique_rows(rows, columns):
    result = {}
    for row in rows:
        key = tuple(row[column] for column in columns)
        if key in result:
            raise ValueError(f"Duplicate study result for {key}")
        result[key] = row
    return result


def load_study(root):
    contrasts = [r for r in records(root / SUMMARY)
                 if r["candidate"] == "no_prior_gat" and r["reference"] == "random_forest"]
    estimates = unique_rows(contrasts, ["direction", "target_stratum", "analysis_level", "metric"])
    for direction, stratum, level, _, units, participants in TARGETS:
        for metric in ("roc_auc", "brier"):
            row = estimates[direction, stratum, level, metric]
            observed_counts = tuple(int(row[k]) for k in ("n_target_samples", "n_target_participants",
                                                           "n_source_splits", "n_bootstrap"))
            if observed_counts != (units, participants, 10, 5000):
                raise ValueError(f"Unexpected evaluation panel for {direction} {stratum} {level} {metric}")
    pairs = unique_rows(records(root / SOURCES / "pooled_auroc_components.csv"),
                        ["case_origin", "control_origin"])
    stages = unique_rows(records(root / SOURCES / "participant_brier_components.csv"), ["stage"])
    if set(pairs) != set(PAIR_KEYS) or set(stages) != {(s,) for s in BRIER_STAGES}:
        raise ValueError("Unexpected origin groups or Brier scoring rules")
    pooled = float(estimates["HMP2_to_PRISM", "pooled", "visit", "roc_auc"]["split_mean_point_delta"])
    if abs(sum(float(row["weighted_delta"]) for row in pairs.values()) - pooled) > 1e-12:
        raise ValueError("Pair contributions do not sum to the observed pooled difference")
    brier_rows = [stages[s,] for s in BRIER_STAGES]
    for level, row in (("visit", brier_rows[0]), ("participant", brier_rows[-1])):
        endpoint = float(estimates["PRISM_to_HMP2", "pooled", level, "brier"]["split_mean_point_delta"])
        if abs(endpoint - float(row["gat_minus_rf"])) > 1e-12:
            raise ValueError(f"Brier endpoint does not match {level}")
    return estimates, [pairs[key] for key in PAIR_KEYS], brier_rows


def signed(value, decimals):
    return f"{value:+.{decimals}f}".replace("-", "−")


def heading(fig, x, y, letter, title, subtitle):
    fig.text(x, y, f"({letter})", weight="bold", size=11, color=INK, va="top")
    fig.text(x+0.042, y, title, weight="bold", size=10.5, color=INK, va="top")
    fig.text(x+0.042, y-0.045, subtitle, size=8.6, color=SUBTLE, va="top")


def style_axes(ax):
    ax.spines[["left", "top", "right"]].set_visible(False)
    ax.spines["bottom"].set_color(GUIDE)
    ax.tick_params(axis="y", length=0, pad=6)
    ax.tick_params(axis="x", color=GUIDE, length=3.5, pad=3)
    ax.set_axisbelow(True)


def forest(ax, estimates, metric, labels):
    style_axes(ax)
    ax.axvline(0, color="#90A0A9", lw=1.1, zorder=0)
    ax.axhline(2.5, color=GUIDE, lw=0.65, ls=(0, (2, 3)), zorder=0)
    for y, (direction, stratum, level, _, _, _) in enumerate(TARGETS):
        r = estimates[direction, stratum, level, metric]
        value, lo, hi = [float(r[k]) for k in ("split_mean_point_delta", "ci_low", "ci_high")]
        color = FORWARD if direction == "HMP2_to_PRISM" else REVERSE
        ax.hlines(y, lo, hi, color=color, lw=1.5)
        ax.vlines([lo, hi], y-0.095, y+0.095, color=color, lw=0.85)
        ax.plot(value, y, "o", ms=5.5, mec=color, mew=1.2,
                mfc="white" if lo <= 0 <= hi else color)
    ax.set_ylim(4.55, -0.55)
    ax.set_yticks(range(5), labels)
    if metric == "roc_auc":
        ax.set_xlim(-0.19, 0.19)
        ax.set_xticks([-0.16, -0.08, 0, 0.08, 0.16])
        ax.set_xlabel("ΔAUROC  (GAT − RF)")
    else:
        ax.set_xlim(-0.060, 0.035)
        ax.set_xticks([-0.06, -0.04, -0.02, 0, 0.02])
        ax.set_xlabel("ΔBrier  (GAT − RF)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: "0" if x == 0 else f"{x:.2f}".replace("-", "−")))


def pair_contributions(ax, pair_rows):
    style_axes(ax)
    group = {"discovery": "PRISM", "validation": "Dutch"}
    labels = []
    for y, row in enumerate(pair_rows):
        value = float(row["weighted_delta"])
        labels.append(f"{group[row['case_origin']]} / {group[row['control_origin']]}")
        ax.hlines(y, 0, value, color=FORWARD, lw=2.2)
        ax.plot(value, y, "o", ms=4.8, color=FORWARD)
        ax.annotate(signed(value, 7), (value, y), xytext=(5, 0), textcoords="offset points",
                    va="center", ha="left", size=8.2, color=INK)
    ax.axvline(0, color=GUIDE, lw=0.8)
    ax.set_ylim(3.55, -0.55)
    ax.set_yticks(range(4), labels)
    ax.set_xlim(-0.001, 0.068)
    ax.set_xticks([0, 0.02, 0.04, 0.06], ["0", "0.02", "0.04", "0.06"])
    ax.set_xlabel("Contribution to pooled ΔAUROC")


def brier_steps(ax, stages):
    style_axes(ax)
    values = [float(row["gat_minus_rf"]) for row in stages]
    ax.axhline(0, lw=0.85, color="#90A0A9", zorder=0)
    ax.plot(range(3), values, lw=1.6, marker="o", ms=5.5, color=REVERSE, zorder=2)
    for i, value in enumerate(values):
        ax.annotate(signed(value, 4), (i, value), xytext=(0, 7 if i else -4.5),
                    va="bottom" if i else "top", ha="center", textcoords="offset points",
                    weight="bold", size=8.7, color=INK)
    for i, row in enumerate(stages[1:]):
        increase = float(row["change_from_previous"])
        if abs(increase - (values[i+1] - values[i])) > 1e-12:
            raise ValueError("Brier change conflicts with the drawn endpoints")
        ax.annotate(signed(increase, 4), (i+0.5, (values[i]+values[i+1])/2),
                    xytext=(1, -6), va="top", ha="center", textcoords="offset points",
                    size=8.1, color=REVERSE, bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
    ax.set_xlim(-0.37, 2.38)
    ax.set_ylim(-0.015, 0.0155)
    ax.set_yticks([-0.01, 0, 0.01], ["−0.01", "0", "+0.01"])
    ax.set_ylabel("ΔBrier  (GAT − RF)", labelpad=1)
    ax.set_xticks(range(3), ["Per record", "Per person\nmean loss", "Per person\nmean probability"])
    ax.tick_params(axis="x", labelsize=8.1)
    ax.tick_params(axis="x", length=0, pad=5)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Repository root (default: current directory).")
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    estimates, pairs, stages = load_study(args.root)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.8, "axes.labelsize": 9.0,
                         "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK,
                         "ytick.color": INK, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig = plt.figure(figsize=(7.8, 5.05), facecolor="white")
    heading(fig, 0.015, 0.975, "a", "Discrimination by target", "Higher ΔAUROC favours GAT")
    heading(fig, 0.548, 0.975, "b", "Probability error by target", "Lower ΔBrier favours GAT")
    ax_a = fig.add_axes((0.188, 0.587, 0.327, 0.280))
    ax_b = fig.add_axes((0.590, 0.587, 0.372, 0.280))
    forest(ax_a, estimates, "roc_auc", [row[3] for row in TARGETS])
    forest(ax_b, estimates, "brier", [""]*len(TARGETS))
    handles = [
        Line2D([], [], color=FORWARD, lw=2, label="HMP2 → Franzosa"),
        Line2D([], [], color=REVERSE, lw=2, label="Franzosa → HMP2"),
        Line2D([], [], color=SUBTLE, marker="o", ms=5.3, mfc="white", lw=0, label="Interval includes zero"),
    ]
    fig.legend(handles=handles, ncol=3, frameon=False, fontsize=8.4,
               loc="center", bbox_to_anchor=(0.50, 0.452), handlelength=1.6,
               columnspacing=2.1, handletextpad=0.55)
    heading(fig, 0.015, 0.400, "c", "Recruitment-group contributions", "Franzosa case / control pairing")
    heading(fig, 0.548, 0.400, "d", "Brier difference by scoring rule", "Fixed predictions on HMP2")
    ax_c = fig.add_axes((0.168, 0.104, 0.347, 0.187))
    ax_d = fig.add_axes((0.640, 0.104, 0.322, 0.187))
    pair_contributions(ax_c, pairs)
    brier_steps(ax_d, stages)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(args.output_dir / f"Figure_5.{suffix}", dpi=600, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
