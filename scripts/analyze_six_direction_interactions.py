#!/usr/bin/env python
"""Formal, paired interaction analysis for the six directed transfers.

The primary inferential unit is the source perturbation split.  Model contrasts
are paired within direction and split, and exact two-sided sign-flip tests are
used for the ten registered perturbations.  A model-by-direction ANOVA is
reported as an exploratory omnibus interaction because the six cohort edges
are the complete directed design, rather than a random sample of cohorts.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


REFERENCE = "elastic_net"
PRIMARY_METRIC = "roc_auc"


def _bh(values: list[float]) -> list[float]:
    values_array = np.asarray(values, dtype=float)
    order = np.argsort(values_array)
    adjusted = np.empty(len(values_array), dtype=float)
    running = 1.0
    for rank in range(len(values_array), 0, -1):
        index = order[rank - 1]
        running = min(running, values_array[index] * len(values_array) / rank)
        adjusted[index] = running
    return adjusted.tolist()


def _exact_sign_flip_p(delta: np.ndarray) -> float:
    delta = np.asarray(delta, dtype=float)
    delta = delta[np.isfinite(delta)]
    if len(delta) == 0:
        return float("nan")
    observed = abs(float(delta.mean()))
    statistics = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(delta)):
        statistics.append(abs(float(np.mean(delta * np.asarray(signs)))))
    return float((1 + np.sum(np.asarray(statistics) >= observed - 1e-15)) / (len(statistics) + 1))


def _load_bootstrap(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path)
    expected = {
        "direction",
        "model",
        "roc_auc_ci_low",
        "roc_auc_ci_high",
        "average_precision_ci_low",
        "average_precision_ci_high",
        "brier_ci_low",
        "brier_ci_high",
    }
    missing = expected.difference(frame.columns)
    if missing:
        raise ValueError(f"Bootstrap summary missing columns: {sorted(missing)}")
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--bootstrap-summary", type=Path, required=False)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    metrics = pd.read_csv(args.metrics)
    required = {"direction", "split", "model", "roc_auc", "average_precision", "brier"}
    missing = required.difference(metrics.columns)
    if missing:
        raise ValueError(f"Metrics missing columns: {sorted(missing)}")
    metrics["split"] = metrics["split"].astype(int)
    directions = sorted(metrics.direction.unique())
    models = sorted(metrics.model.unique())

    direction_ranks = {
        direction: metrics.loc[metrics.direction == direction].pivot(
            index="split", columns="model", values=PRIMARY_METRIC
        ).rank(axis=1, ascending=False, method="average")
        for direction in directions
    }
    aggregate_rows = []
    for (direction, model), group in metrics.groupby(["direction", "model"], sort=True):
        row = {"direction": direction, "model": model, "n_source_splits": int(group.split.nunique())}
        for metric in ("roc_auc", "average_precision", "brier"):
            row[f"{metric}_mean"] = float(group[metric].mean())
            row[f"{metric}_sd"] = float(group[metric].std(ddof=1))
        model_ranks = direction_ranks[direction][model]
        row["rank_mean"] = float(model_ranks.mean())
        row["top_rank_frequency"] = float((model_ranks == 1).mean())
        aggregate_rows.append(row)
    aggregate = pd.DataFrame(aggregate_rows)

    bootstrap = _load_bootstrap(args.bootstrap_summary) if args.bootstrap_summary else pd.DataFrame()
    if not bootstrap.empty:
        aggregate = aggregate.merge(
            bootstrap[
                [
                    "direction",
                    "model",
                    "roc_auc_ci_low",
                    "roc_auc_ci_high",
                    "average_precision_ci_low",
                    "average_precision_ci_high",
                    "brier_ci_low",
                    "brier_ci_high",
                ]
            ],
            on=["direction", "model"],
            how="left",
        )
    aggregate.to_csv(args.output_dir / "six_direction_aggregate.csv", index=False)

    contrasts = []
    pvalue_index = {"roc_auc": [], "average_precision": [], "brier": []}
    for direction in directions:
        reference = metrics.loc[
            (metrics.direction == direction) & (metrics.model == REFERENCE)
        ].set_index("split")
        for model in models:
            if model == REFERENCE:
                continue
            candidate = metrics.loc[
                (metrics.direction == direction) & (metrics.model == model)
            ].set_index("split")
            for metric in ("roc_auc", "average_precision", "brier"):
                delta = (candidate[metric] - reference[metric]).dropna().to_numpy(float)
                pvalue = _exact_sign_flip_p(delta)
                contrasts.append(
                    {
                        "direction": direction,
                        "candidate": model,
                        "reference": REFERENCE,
                        "metric": metric,
                        "n_splits": len(delta),
                        "delta_mean": float(delta.mean()),
                        "delta_sd": float(delta.std(ddof=1)),
                        "p_exact_sign_flip": pvalue,
                    }
                )
                pvalue_index[metric].append(len(contrasts) - 1)
    contrast_frame = pd.DataFrame(contrasts)
    contrast_frame["q_bh_within_metric"] = np.nan
    for metric, indices in pvalue_index.items():
        contrast_frame.loc[indices, "q_bh_within_metric"] = _bh(
            contrast_frame.loc[indices, "p_exact_sign_flip"].tolist()
        )
    contrast_frame.to_csv(args.output_dir / "six_direction_model_contrasts.csv", index=False)

    rank = aggregate.pivot(index="direction", columns="model", values="rank_mean")
    rank.to_csv(args.output_dir / "six_direction_rank_means.csv")
    top_frequency = aggregate.pivot(index="direction", columns="model", values="top_rank_frequency")
    top_frequency.to_csv(args.output_dir / "six_direction_top_rank_frequency.csv")

    interaction = {}
    try:
        import statsmodels.api as sm
        from statsmodels.formula.api import ols

        model = ols("roc_auc ~ C(model) * C(direction)", data=metrics).fit()
        anova = sm.stats.anova_lm(model, typ=2)
        anova.reset_index().rename(columns={"index": "term"}).to_csv(
            args.output_dir / "six_direction_interaction_anova.csv", index=False
        )
        term = "C(model):C(direction)"
        interaction = {
            "method": "ordinary_least_squares_type_II_anova",
            "unit": "source_perturbation_split",
            "interaction_f": float(anova.loc[term, "F"]),
            "interaction_p": float(anova.loc[term, "PR(>F)"]),
            "caveat": "Exploratory omnibus test; the six cohort edges are a complete directed design, not random cohort sampling.",
        }
    except Exception as error:  # pragma: no cover - environment-dependent fallback
        interaction = {"method": "unavailable", "error": repr(error)}

    summary = {
        "directions": directions,
        "models": models,
        "reference_model": REFERENCE,
        "primary_metric": PRIMARY_METRIC,
        "n_direction_model_comparisons": int(len(contrast_frame)),
        "multiplicity": "Benjamini-Hochberg within metric across six directions and six non-reference models",
        "interaction": interaction,
        "interpretation_guardrail": "Direction-level contrasts are paired across source perturbations; target participant bootstrap CIs are reported separately and do not replace source-perturbation uncertainty.",
    }
    (args.output_dir / "six_direction_interaction_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
