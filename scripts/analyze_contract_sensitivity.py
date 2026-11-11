#!/usr/bin/env python
"""Summarize registered feature-contract sensitivity without selecting a winner."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def _jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return float(len(left & right) / len(union)) if union else float("nan")


def _feature_stability(features: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows, frequencies = [], []
    key_columns = ["variant_id", "feature_cap", "min_prevalence", "direction"]
    for keys, frame in features.groupby(key_columns, sort=True):
        sets = {split: set(part.feature_id) for split, part in frame.groupby("split")}
        scores = [_jaccard(sets[left], sets[right]) for left, right in itertools.combinations(sets, 2)]
        summary_rows.append(
            dict(zip(key_columns, keys))
            | {
                "n_splits": len(sets),
                "pairwise_jaccard_mean": float(np.mean(scores)),
                "pairwise_jaccard_sd": float(np.std(scores, ddof=1)),
            }
        )
        counts = frame.groupby("feature_id").split.nunique().reset_index(name="n_selected_splits")
        counts["selection_frequency"] = counts.n_selected_splits / len(sets)
        for column, value in zip(key_columns, keys):
            counts[column] = value
        frequencies.append(counts)
    return pd.DataFrame(summary_rows), pd.concat(frequencies, ignore_index=True)


def _weight_stability(weights: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["variant_id", "feature_cap", "min_prevalence", "direction", "model"]
    for values, frame in weights.groupby(keys, sort=True):
        matrix = frame.pivot(index="split", columns="feature_id", values="feature_weight").fillna(0.0)
        correlations = []
        signs = []
        for left, right in itertools.combinations(matrix.index, 2):
            first, second = matrix.loc[left], matrix.loc[right]
            corr = spearmanr(first, second).statistic
            correlations.append(float(corr))
            active = (first != 0) & (second != 0)
            if active.any():
                signs.append(float((np.sign(first[active]) == np.sign(second[active])).mean()))
        rows.append(
            dict(zip(keys, values))
            | {
                "n_splits": int(matrix.shape[0]),
                "pairwise_weight_spearman_mean": float(np.nanmean(correlations)),
                "pairwise_weight_spearman_sd": float(np.nanstd(correlations, ddof=1)),
                "shared_nonzero_sign_agreement_mean": float(np.nanmean(signs)) if signs else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--weights", type=Path, required=True)
    parser.add_argument("--main-metrics", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    metrics = pd.read_csv(args.metrics)
    main_metrics = pd.read_csv(args.main_metrics)
    features = pd.read_csv(args.features)
    weights = pd.read_csv(args.weights)
    aggregate = metrics.groupby(["variant_id", "feature_cap", "min_prevalence", "direction", "model"], as_index=False)[
        ["roc_auc", "average_precision", "brier"]
    ].agg(["mean", "std"])
    aggregate.columns = [
        "_".join(column).strip("_") if isinstance(column, tuple) else column
        for column in aggregate.columns
    ]
    registered = main_metrics.groupby(["direction", "model"], as_index=False)[
        ["roc_auc", "average_precision", "brier"]
    ].mean().rename(columns={metric: f"registered_{metric}" for metric in ("roc_auc", "average_precision", "brier")})
    aggregate = aggregate.merge(registered, on=["direction", "model"], how="left")
    for metric in ("roc_auc", "average_precision", "brier"):
        aggregate[f"delta_{metric}_vs_registered"] = aggregate[f"{metric}_mean"] - aggregate[f"registered_{metric}"]
    aggregate.to_csv(out / "contract_sensitivity_aggregate.csv", index=False)
    feature_summary, frequency = _feature_stability(features)
    feature_summary.to_csv(out / "contract_sensitivity_feature_jaccard.csv", index=False)
    frequency.to_csv(out / "contract_sensitivity_feature_selection_frequency.csv", index=False)
    weight_summary = _weight_stability(weights)
    weight_summary.to_csv(out / "contract_sensitivity_weight_stability.csv", index=False)
    ranking = aggregate.assign(
        rank=aggregate.groupby(["variant_id", "direction"])["roc_auc_mean"].rank(ascending=False, method="average")
    )
    ranking.to_csv(out / "contract_sensitivity_model_ranks.csv", index=False)
    summary = {
        "analysis": "registered_feature_contract_sensitivity",
        "variants": sorted(metrics.variant_id.unique().tolist()),
        "directions": sorted(metrics.direction.unique().tolist()),
        "models": sorted(metrics.model.unique().tolist()),
        "n_metric_rows": int(len(metrics)),
        "target_labels_used_for_model_or_contract_selection": False,
        "interpretation": "No sensitivity variant is selected as a replacement primary benchmark; all differences are reported against the registered 126@0.10 protocol.",
    }
    (out / "contract_sensitivity_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
