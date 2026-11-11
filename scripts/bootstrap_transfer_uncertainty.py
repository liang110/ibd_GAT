#!/usr/bin/env python
"""Cluster-bootstrap external transfer predictions saved by the CPU baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


def _metrics(frame: pd.DataFrame, sample_weight: np.ndarray | None = None) -> dict[str, float]:
    return {
        "roc_auc": roc_auc_score(frame.outcome_ibd, frame.probability, sample_weight=sample_weight),
        "average_precision": average_precision_score(frame.outcome_ibd, frame.probability, sample_weight=sample_weight),
        "brier": brier_score_loss(frame.outcome_ibd, frame.probability, sample_weight=sample_weight),
    }


def _resample_cluster_weights(frame: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    groups = frame.bootstrap_cluster.drop_duplicates().to_numpy()
    for _ in range(100):
        sampled = rng.choice(groups, size=len(groups), replace=True)
        counts = pd.Series(sampled).value_counts()
        weights = frame.bootstrap_cluster.map(counts).fillna(0).to_numpy(dtype=float)
        if frame.loc[weights > 0, "outcome_ibd"].nunique() == 2:
            return weights
    raise RuntimeError("Cluster bootstrap could not retain both outcomes")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260815)
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    predictions = pd.read_csv(args.predictions)
    rng = np.random.default_rng(args.seed)
    summary_rows, draws = [], []
    for keys, frame in predictions.groupby(["direction", "split", "model"], sort=True):
        observed = _metrics(frame)
        boot = []
        for replicate in range(args.n_bootstrap):
            metrics = _metrics(frame, _resample_cluster_weights(frame, rng))
            boot.append(metrics)
            draws.append({"direction": keys[0], "split": keys[1], "model": keys[2], "replicate": replicate, **metrics})
        boot_frame = pd.DataFrame(boot)
        summary_rows.append({"direction": keys[0], "split": keys[1], "model": keys[2], "n_samples": len(frame), "n_clusters": frame.bootstrap_cluster.nunique(), **{f"{metric}_observed": value for metric, value in observed.items()}, **{f"{metric}_ci_low": float(boot_frame[metric].quantile(0.025)) for metric in observed}, **{f"{metric}_ci_high": float(boot_frame[metric].quantile(0.975)) for metric in observed}})
    summary = pd.DataFrame(summary_rows)
    draw_frame = pd.DataFrame(draws)
    summary.to_csv(out / "external_transfer_cluster_bootstrap.csv", index=False)
    draw_frame.to_csv(out / "external_transfer_cluster_bootstrap_draws.csv", index=False)
    aggregate = summary.groupby(["direction", "model"])[["roc_auc_observed", "average_precision_observed", "brier_observed"]].agg(["mean", "std"]).round(6)
    aggregate.columns = [f"{metric}_{stat}" for metric, stat in aggregate.columns]
    aggregate_draws = draw_frame.groupby(
        ["direction", "model", "replicate"], as_index=False
    )[["roc_auc", "average_precision", "brier"]].mean()
    aggregate_ci_rows = []
    for (direction, model), frame in aggregate_draws.groupby(["direction", "model"]):
        row = {"direction": direction, "model": model}
        for metric in ("roc_auc", "average_precision", "brier"):
            row[f"{metric}_bootstrap_mean"] = float(frame[metric].mean())
            row[f"{metric}_ci_low"] = float(frame[metric].quantile(0.025))
            row[f"{metric}_ci_high"] = float(frame[metric].quantile(0.975))
        aggregate_ci_rows.append(row)
    aggregate_ci = pd.DataFrame(aggregate_ci_rows)
    aggregate_ci.to_csv(out / "aggregate_cluster_bootstrap.csv", index=False)
    report = {
        "n_bootstrap": args.n_bootstrap,
        "seed": args.seed,
        "groups": int(len(summary)),
        "aggregate": aggregate.reset_index().to_dict(orient="records"),
        "aggregate_cluster_bootstrap": aggregate_ci.round(6).to_dict(orient="records"),
    }
    (out / "bootstrap_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
