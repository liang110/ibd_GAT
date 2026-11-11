#!/usr/bin/env python
"""Paired target-cluster bootstrap of external model differences."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


def _metric_delta(frame: pd.DataFrame, weights: np.ndarray | None = None) -> dict[str, float]:
    y = frame.outcome_ibd
    candidate = frame.probability_candidate
    reference = frame.probability_reference
    return {
        "roc_auc_delta": float(roc_auc_score(y, candidate, sample_weight=weights) - roc_auc_score(y, reference, sample_weight=weights)),
        "average_precision_delta": float(average_precision_score(y, candidate, sample_weight=weights) - average_precision_score(y, reference, sample_weight=weights)),
        "brier_delta": float(brier_score_loss(y, candidate, sample_weight=weights) - brier_score_loss(y, reference, sample_weight=weights)),
    }


def _cluster_weights(frame: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    clusters = frame.bootstrap_cluster.drop_duplicates().to_numpy()
    for _ in range(100):
        sampled = rng.choice(clusters, size=len(clusters), replace=True)
        counts = pd.Series(sampled).value_counts()
        weights = frame.bootstrap_cluster.map(counts).fillna(0).to_numpy(dtype=float)
        if frame.loc[weights > 0, "outcome_ibd"].nunique() == 2:
            return weights
    raise RuntimeError("Bootstrap could not retain both outcomes")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--reference", default="elastic_net")
    parser.add_argument("--n-bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260816)
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    predictions = pd.read_csv(args.predictions)
    if args.reference not in set(predictions.model):
        raise ValueError(f"Reference model {args.reference} absent")
    rng = np.random.default_rng(args.seed)
    rows, draws = [], []
    group_columns = ["direction", "split"]
    for keys, group in predictions.groupby(group_columns, sort=True):
        reference = group.loc[group.model == args.reference].drop(columns="model")
        for candidate_name, candidate in group.groupby("model", sort=True):
            if candidate_name == args.reference:
                continue
            joined = candidate.drop(columns="model").merge(
                reference,
                on=["direction", "split", "sample_id", "outcome_ibd", "bootstrap_cluster", "origin"],
                suffixes=("_candidate", "_reference"),
                validate="one_to_one",
            )
            observed = _metric_delta(joined)
            bootstrap = []
            for replicate in range(args.n_bootstrap):
                delta = _metric_delta(joined, _cluster_weights(joined, rng))
                bootstrap.append(delta)
                draws.append({"direction": keys[0], "split": keys[1], "candidate": candidate_name, "reference": args.reference, "replicate": replicate, **delta})
            boot = pd.DataFrame(bootstrap)
            rows.append({
                "direction": keys[0],
                "split": keys[1],
                "candidate": candidate_name,
                "reference": args.reference,
                **{f"{metric}_observed": value for metric, value in observed.items()},
                **{f"{metric}_ci_low": float(boot[metric].quantile(0.025)) for metric in observed},
                **{f"{metric}_ci_high": float(boot[metric].quantile(0.975)) for metric in observed},
            })
    pd.DataFrame(rows).to_csv(out / "paired_model_bootstrap.csv", index=False)
    pd.DataFrame(draws).to_csv(out / "paired_model_bootstrap_draws.csv", index=False)
    report = {"reference": args.reference, "n_bootstrap": args.n_bootstrap, "comparisons": len(rows)}
    (out / "paired_model_bootstrap_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
