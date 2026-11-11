#!/usr/bin/env python
"""Write descriptive target-evaluation rank stability from saved metrics."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for (direction, fraction), frame in pd.read_csv(args.metrics).groupby(["direction", "fraction"], sort=True):
        ranks = frame.copy()
        ranks["rank"] = ranks.groupby("split").roc_auc.rank(ascending=False, method="average")
        winner_frequency = ranks.assign(is_first=ranks["rank"] == 1).groupby("model").is_first.mean()
        rows.append({
            "direction": direction,
            "fraction": fraction,
            "top_model_by_mean_auroc": frame.groupby("model").roc_auc.mean().idxmax(),
            "top_model_first_rank_frequency": float(winner_frequency.max()),
            "n_models_with_first_rank": int((winner_frequency > 0).sum()),
            "rank_sd_mean": float(ranks.groupby("model")["rank"].std(ddof=1).mean()),
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
