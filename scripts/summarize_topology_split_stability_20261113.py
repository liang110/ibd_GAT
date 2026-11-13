from __future__ import annotations

"""Create a split-level audit for the topology intervention.

The primary estimand is defined before aggregation: within each direction and
source split, the five degree-preserving rewired graphs are averaged and are
then contrasted with the single real Spearman graph.  This avoids treating the
five rewires as independent source splits.  The output is intentionally small
and machine-readable so that the manuscript and review package use the same
numbers.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "artifacts" / "hpc_results" / "gateB_edge_ablation_20261112" / "edge_ablation_metrics.csv"
OUT = ROOT / "artifacts" / "reconciled_submission_20261113"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> None:
    df = pd.read_csv(INPUT)
    required = {"direction", "split", "condition", "roc_auc", "average_precision", "brier"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required topology columns: {sorted(missing)}")

    rows = []
    for (direction, split), grp in df.groupby(["direction", "split"], sort=True):
        real = grp.loc[grp.condition.eq("real_spearman")].iloc[0]
        identity = grp.loc[grp.condition.eq("identity")].iloc[0]
        rewired = grp.loc[grp.condition.eq("rewired")]
        if len(rewired) != 5:
            raise ValueError(f"Expected five rewired records for {direction} split {split}, got {len(rewired)}")
        rows.append({
            "direction": direction,
            "split": int(split),
            "real_auroc": float(real.roc_auc),
            "rewired_auroc_mean": float(rewired.roc_auc.mean()),
            "rewired_auroc_sd": float(rewired.roc_auc.std(ddof=1)),
            "identity_auroc": float(identity.roc_auc),
            "real_minus_rewired_auroc": float(real.roc_auc - rewired.roc_auc.mean()),
            "real_minus_identity_auroc": float(real.roc_auc - identity.roc_auc),
            "real_ap": float(real.average_precision),
            "rewired_ap_mean": float(rewired.average_precision.mean()),
            "real_minus_rewired_ap": float(real.average_precision - rewired.average_precision.mean()),
            "real_brier": float(real.brier),
            "rewired_brier_mean": float(rewired.brier.mean()),
            "real_minus_rewired_brier": float(real.brier - rewired.brier.mean()),
            "n_rewired": int(len(rewired)),
            "positive_real_minus_rewired": bool(real.roc_auc > rewired.roc_auc.mean()),
            "edge_count_real": int(real.edge_count),
            "edge_count_rewired_median": float(rewired.edge_count.median()),
            "out_degree_variance_real": float(real.out_degree_variance),
            "in_degree_variance_real": float(real.in_degree_variance),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "topology_split_stability.csv", index=False)

    summary = {}
    for direction, grp in out.groupby("direction", sort=True):
        vals = grp.real_minus_rewired_auroc.to_numpy(float)
        loo = [float(np.delete(vals, i).mean()) for i in range(len(vals))]
        summary[direction] = {
            "n_splits": int(len(vals)),
            "mean_split_equal": float(vals.mean()),
            "median_split_equal": float(np.median(vals)),
            "min_split": float(vals.min()),
            "max_split": float(vals.max()),
            "positive_splits": int((vals > 0).sum()),
            "positive_fraction": float((vals > 0).mean()),
            "leave_one_split_out_min": float(min(loo)),
            "leave_one_split_out_max": float(max(loo)),
            "rewired_records_per_split": 5,
            "estimand": "split-equal mean of real AUROC minus within-split mean of five rewired AUROCs",
        }
    (OUT / "topology_split_stability.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lines = [
        "# Topology split-stability audit",
        "",
        "The five rewired graphs are averaged within each source split before the real-minus-rewired contrast is calculated.",
        "",
        "| Direction | Mean | Median | Range | Positive splits | Leave-one-split-out mean range |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for direction, s in summary.items():
        lines.append(
            f"| {direction} | {s['mean_split_equal']:.6f} | {s['median_split_equal']:.6f} | "
            f"{s['min_split']:.6f} to {s['max_split']:.6f} | {s['positive_splits']}/{s['n_splits']} | "
            f"{s['leave_one_split_out_min']:.6f} to {s['leave_one_split_out_max']:.6f} |"
        )
    (OUT / "topology_split_stability.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT / "topology_split_stability.csv")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
