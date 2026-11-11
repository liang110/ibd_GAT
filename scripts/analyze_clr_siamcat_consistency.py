#!/usr/bin/env python
"""Compare CLR and SIAMCAT only at the pre-defined protocol-comparison layer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr


OVERLAPPING_DIRECTIONS = (
    "HMP2_to_PRISM",
    "PRISM_to_HMP2",
    "HMP2_to_NielsenHB_2014",
    "PRISM_to_NielsenHB_2014",
)
CLR_MODELS = ("logistic_l2", "elastic_net", "linear_svm")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clr-metrics", type=Path, required=True)
    parser.add_argument("--clr-predictions", type=Path, required=True)
    parser.add_argument("--siamcat-metrics", type=Path, required=True)
    parser.add_argument("--siamcat-predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    clr_metrics = pd.read_csv(args.clr_metrics)
    siamcat_metrics = pd.read_csv(args.siamcat_metrics)
    clr_predictions = pd.read_csv(args.clr_predictions)
    siamcat_predictions = pd.read_csv(args.siamcat_predictions)
    clr_metrics = clr_metrics[clr_metrics.direction.isin(OVERLAPPING_DIRECTIONS) & clr_metrics.model.isin(CLR_MODELS)]
    siamcat_metrics = siamcat_metrics[siamcat_metrics.direction.isin(OVERLAPPING_DIRECTIONS)]
    metric_rows = []
    for direction in OVERLAPPING_DIRECTIONS:
        reference = siamcat_metrics[siamcat_metrics.direction == direction]
        for model in CLR_MODELS:
            candidate = clr_metrics[(clr_metrics.direction == direction) & (clr_metrics.model == model)]
            for metric in ("roc_auc", "average_precision", "brier"):
                metric_rows.append(
                    {
                        "direction": direction,
                        "clr_model": model,
                        "siamcat_model": "siamcat_lasso",
                        "metric": metric,
                        "clr_mean": float(candidate[metric].mean()),
                        "siamcat_mean": float(reference[metric].mean()),
                        "difference_clr_minus_siamcat": float(candidate[metric].mean() - reference[metric].mean()),
                    }
                )
    metric_table = pd.DataFrame(metric_rows)
    metric_table.to_csv(args.output_dir / "clr_siamcat_direction_metrics.csv", index=False)

    correlations = []
    for direction in OVERLAPPING_DIRECTIONS:
        siam = siamcat_predictions[siamcat_predictions.direction == direction]
        for model in CLR_MODELS:
            clr = clr_predictions[(clr_predictions.direction == direction) & (clr_predictions.model == model)]
            for split in sorted(set(clr.split).intersection(siam.split)):
                left = clr[clr.split == split][["sample_id", "probability"]].rename(columns={"probability": "clr_probability"})
                right = siam[siam.split == split][["sample_id", "probability"]].rename(columns={"probability": "siamcat_probability"})
                merged = left.merge(right, on="sample_id", validate="one_to_one")
                correlations.append(
                    {
                        "direction": direction,
                        "clr_model": model,
                        "split": split,
                        "n_samples": len(merged),
                        "pearson_probability_correlation": float(pearsonr(merged.clr_probability, merged.siamcat_probability).statistic),
                        "spearman_probability_correlation": float(spearmanr(merged.clr_probability, merged.siamcat_probability).statistic),
                    }
                )
    corr = pd.DataFrame(correlations)
    corr.to_csv(args.output_dir / "clr_siamcat_prediction_correlations.csv", index=False)
    summary = corr.groupby(["direction", "clr_model"], as_index=False)[
        ["pearson_probability_correlation", "spearman_probability_correlation"]
    ].agg(["mean", "std"])
    summary.columns = ["_".join(c).strip("_") if isinstance(c, tuple) else c for c in summary.columns]
    summary.to_csv(args.output_dir / "clr_siamcat_prediction_correlation_summary.csv", index=False)
    direction_rank_rows = []
    for model in CLR_MODELS:
        merged = metric_table[(metric_table.clr_model == model) & (metric_table.metric == "roc_auc")]
        direction_rank_rows.append(
            {
                "clr_model": model,
                "n_directions": len(merged),
                "spearman_rank_correlation_of_direction_auroc": float(
                    spearmanr(merged.clr_mean.rank(), merged.siamcat_mean.rank()).statistic
                ),
            }
        )
    pd.DataFrame(direction_rank_rows).to_csv(args.output_dir / "clr_siamcat_direction_rank_consistency.csv", index=False)
    report = {
        "comparison_scope": "protocol-level consistency only; different preprocessing pipelines are not pooled into one model ranking",
        "overlapping_directions": list(OVERLAPPING_DIRECTIONS),
        "clr_models": list(CLR_MODELS),
        "target_labels_used_for_comparison_only": True,
    }
    (args.output_dir / "clr_siamcat_consistency_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
