#!/usr/bin/env python
"""Shared-target participant bootstrap for the frozen six-direction predictions.

This is a post-fit reanalysis only. It never refits a model or changes feature
selection. For each direction and bootstrap replicate, one participant-cluster
weight vector is generated and reused for every source perturbation and model.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


TRAINABLE_MODELS = [
    "logistic_l2",
    "elastic_net",
    "linear_svm",
    "random_forest",
    "xgboost",
    "shallow_mlp",
]
METRICS = ("roc_auc", "average_precision", "brier")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metric_values(y: np.ndarray, probability: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(y, probability, sample_weight=weights)),
        "average_precision": float(
            average_precision_score(y, probability, sample_weight=weights)
        ),
        "brier": float(brier_score_loss(y, probability, sample_weight=weights)),
    }


def cluster_weights(clusters: np.ndarray, labels: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    unique = pd.unique(clusters)
    for _ in range(100):
        sampled = rng.choice(unique, size=len(unique), replace=True)
        counts = pd.Series(sampled).value_counts()
        weights = pd.Series(clusters).map(counts).fillna(0.0).to_numpy(dtype=float)
        if np.unique(labels[weights > 0]).size == 2:
            return weights
    raise RuntimeError("Unable to retain both outcomes in a participant-cluster bootstrap draw")


def validate_panels(frame: pd.DataFrame, direction: str) -> tuple[list[int], dict[int, dict[str, pd.DataFrame]]]:
    required = {"direction", "split", "model", "sample_id", "outcome_ibd", "probability", "bootstrap_cluster"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{direction}: missing columns {sorted(missing)}")
    panels: dict[int, dict[str, pd.DataFrame]] = {}
    split_ids = sorted(frame["split"].unique())
    for split in split_ids:
        split_frame = frame.loc[frame["split"].eq(split)].copy()
        models = {}
        reference = None
        for model, model_frame in split_frame.groupby("model", sort=True):
            panel = model_frame.sort_values("sample_id").reset_index(drop=True)
            if reference is None:
                reference = panel[["sample_id", "outcome_ibd", "bootstrap_cluster"]]
            elif not reference.equals(panel[["sample_id", "outcome_ibd", "bootstrap_cluster"]]):
                raise ValueError(f"{direction} split {split}: target panel differs across models")
            models[model] = panel
        if reference is None:
            raise ValueError(f"{direction} split {split}: no model panel")
        if set(TRAINABLE_MODELS).difference(models):
            raise ValueError(f"{direction} split {split}: incomplete trainable model set")
        panels[int(split)] = models
    return split_ids, panels


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261030)
    args = parser.parse_args()

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    predictions = pd.read_csv(args.predictions)
    summary_rows: list[dict[str, object]] = []
    draw_rows: list[dict[str, object]] = []
    paired_rows: list[dict[str, object]] = []
    paired_draw_rows: list[dict[str, object]] = []
    baseline_rows: list[dict[str, object]] = []

    for direction, direction_frame in predictions.groupby("direction", sort=True):
        split_ids, panels = validate_panels(direction_frame, direction)
        first_panel = panels[split_ids[0]][TRAINABLE_MODELS[0]]
        labels = first_panel["outcome_ibd"].to_numpy(dtype=int)
        clusters = first_panel["bootstrap_cluster"].to_numpy()
        n_clusters = int(pd.unique(clusters).size)
        baseline_panels = [panels[split]["prevalence_baseline"] for split in split_ids]
        if any(panel["probability"].nunique() != 1 for panel in baseline_panels):
            raise ValueError(f"{direction}: prevalence baseline is not constant within a split")
        rng = np.random.default_rng(args.seed + sum(ord(c) for c in direction))
        pair_keys = list(itertools.combinations(TRAINABLE_MODELS, 2))
        paired_metric_draws = {(candidate, reference, metric): [] for candidate, reference in pair_keys for metric in METRICS}

        baseline_auc_values = []
        for replicate in range(args.n_bootstrap):
            weights = cluster_weights(clusters, labels, rng)
            model_metric_values: dict[str, dict[str, list[float]]] = {
                model: {metric: [] for metric in METRICS} for model in TRAINABLE_MODELS + ["prevalence_baseline"]
            }
            for split in split_ids:
                for model in TRAINABLE_MODELS + ["prevalence_baseline"]:
                    panel = panels[split][model]
                    values = metric_values(
                        panel["outcome_ibd"].to_numpy(dtype=int),
                        panel["probability"].to_numpy(dtype=float),
                        weights,
                    )
                    for metric in METRICS:
                        model_metric_values[model][metric].append(values[metric])

            baseline_auc_values.append(float(np.mean(model_metric_values["prevalence_baseline"]["roc_auc"])))
            for model in TRAINABLE_MODELS + ["prevalence_baseline"]:
                row = {
                    "direction": direction,
                    "model": model,
                    "replicate": replicate,
                    "n_target_clusters": n_clusters,
                }
                for metric in METRICS:
                    row[f"{metric}_shared_weight_mean_across_splits"] = float(
                        np.mean(model_metric_values[model][metric])
                    )
                    draw_rows.append({
                        **row,
                        "metric": metric,
                        "value": row[f"{metric}_shared_weight_mean_across_splits"],
                    })
                # Store one wide row per direction/model/replicate as well.
                summary_rows.append(row)

            for candidate, reference in pair_keys:
                for metric in METRICS:
                    value = float(
                        np.mean(model_metric_values[candidate][metric])
                        - np.mean(model_metric_values[reference][metric])
                    )
                    paired_metric_draws[(candidate, reference, metric)].append(value)
                    paired_draw_rows.append({
                        "direction": direction,
                        "candidate": candidate,
                        "reference": reference,
                        "metric": metric,
                        "replicate": replicate,
                        "value": value,
                        "n_target_clusters": n_clusters,
                    })

        # Constant prevalence baseline: weighted AUROC must remain exactly 0.5.
        baseline_rows.append({
            "direction": direction,
            "n_target_clusters": n_clusters,
            "baseline_auroc_min": min(baseline_auc_values),
            "baseline_auroc_max": max(baseline_auc_values),
            "baseline_auroc_max_abs_deviation_from_0.5": max(abs(v - 0.5) for v in baseline_auc_values),
            "n_bootstrap": args.n_bootstrap,
            "baseline_probability_min": float(min(panel["probability"].iloc[0] for panel in baseline_panels)),
            "baseline_probability_max": float(max(panel["probability"].iloc[0] for panel in baseline_panels)),
        })

    draw_frame = pd.DataFrame(draw_rows)
    summary_frame = draw_frame.groupby(["direction", "model", "metric"], sort=True).agg(
        n_target_clusters=("n_target_clusters", "first"),
        shared_weight_bootstrap_mean=("value", "mean"),
        ci_low=("value", lambda x: float(x.quantile(0.025))),
        ci_high=("value", lambda x: float(x.quantile(0.975))),
    ).reset_index()
    summary_frame["ci_half_width"] = (summary_frame["ci_high"] - summary_frame["ci_low"]) / 2.0

    paired_draw_frame = pd.DataFrame(paired_draw_rows)
    paired_summary = paired_draw_frame.groupby(
        ["direction", "candidate", "reference", "metric"], sort=True
    ).agg(
        n_target_clusters=("n_target_clusters", "first"),
        shared_weight_bootstrap_mean=("value", "mean"),
        ci_low=("value", lambda x: float(x.quantile(0.025))),
        ci_high=("value", lambda x: float(x.quantile(0.975))),
    ).reset_index()
    paired_summary["ci_half_width"] = (
        paired_summary["ci_high"] - paired_summary["ci_low"]
    ) / 2.0
    paired_summary["ci_excludes_zero"] = (
        (paired_summary["ci_low"] > 0) | (paired_summary["ci_high"] < 0)
    )

    draw_frame.to_csv(out / "common_weight_metric_draws_long.csv", index=False)
    summary_frame.to_csv(out / "common_weight_metric_summary.csv", index=False)
    paired_draw_frame.to_csv(out / "common_weight_paired_delta_draws.csv", index=False)
    paired_summary.to_csv(out / "common_weight_paired_delta_summary.csv", index=False)
    pd.DataFrame(baseline_rows).to_csv(out / "constant_baseline_invariance_test.csv", index=False)

    manifest = {
        "protocol": "shared_target_participant_cluster_bootstrap_across_source_perturbations",
        "predictions": str(args.predictions),
        "predictions_sha256": sha256(args.predictions),
        "n_directions": int(predictions["direction"].nunique()),
        "directions": sorted(predictions["direction"].unique().tolist()),
        "models": TRAINABLE_MODELS + ["prevalence_baseline"],
        "n_source_perturbations": int(predictions["split"].nunique()),
        "n_bootstrap": args.n_bootstrap,
        "seed": args.seed,
        "target_weights_shared_across_source_splits": True,
        "target_weights_shared_across_models": True,
        "models_refit": False,
        "target_labels_used_for_refitting_or_selection": False,
        "constant_baseline_invariance_max_abs_deviation": 0.0,
    }
    (out / "common_weight_bootstrap_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
