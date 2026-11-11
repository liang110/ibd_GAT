#!/usr/bin/env python
"""Complete frozen-prediction ESP target audit for the two Nielsen-target edges.

No models are refit. One participant-cluster bootstrap weight vector per
replicate is shared across all source perturbations and all model families.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss


MODELS = [
    "prevalence_baseline",
    "logistic_l2",
    "elastic_net",
    "linear_svm",
    "random_forest",
    "xgboost",
    "shallow_mlp",
]
TRAINABLE = MODELS[1:]
METRICS = ("roc_auc", "average_precision", "brier")
TARGET_DIRECTIONS = ("HMP2_to_NielsenHB_2014", "PRISM_to_NielsenHB_2014")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def weighted_auc(y: np.ndarray, score: np.ndarray, w: np.ndarray) -> float:
    order = np.argsort(score, kind="mergesort")
    s = score[order]
    yy = y[order]
    ww = w[order]
    neg = ww * (1 - yy)
    pos = ww * yy
    unique_start = np.r_[0, np.flatnonzero(np.diff(s)) + 1]
    unique_end = np.r_[unique_start[1:], len(s)]
    cum_neg = np.cumsum(neg)
    numerator = 0.0
    for start, end in zip(unique_start, unique_end):
        p = pos[start:end].sum()
        n = neg[start:end].sum()
        before = cum_neg[start] - neg[start]
        numerator += p * (before + 0.5 * n)
    total_pos = pos.sum()
    total_neg = neg.sum()
    return float(numerator / (total_pos * total_neg))


def bootstrap_weights(clusters: np.ndarray, labels: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    unique = pd.unique(clusters)
    for _ in range(100):
        sampled = rng.choice(unique, size=len(unique), replace=True)
        counts = pd.Series(sampled).value_counts()
        w = pd.Series(clusters).map(counts).fillna(0).to_numpy(dtype=float)
        if np.unique(labels[w > 0]).size == 2:
            return w
    raise RuntimeError("ESP bootstrap failed to retain both outcomes")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261030)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    predictions = pd.read_csv(args.predictions)
    manifest = pd.read_csv(args.manifest)
    manifest["country"] = manifest["country"].fillna("unknown").astype(str)
    esp_ids = set(manifest.loc[manifest["country"].eq("ESP"), "sample_id"])
    predictions = predictions.loc[
        predictions["direction"].isin(TARGET_DIRECTIONS)
        & predictions["sample_id"].isin(esp_ids)
        & predictions["model"].isin(MODELS)
    ].copy()
    if predictions.empty:
        raise ValueError("No ESP target predictions found")
    target_meta = manifest.loc[manifest["sample_id"].isin(esp_ids)].set_index("sample_id")
    target_meta = target_meta.loc[predictions["sample_id"].drop_duplicates()]
    if target_meta["outcome_ibd"].nunique() != 2:
        raise ValueError("ESP target subset does not contain both outcomes")
    rows = []
    draws = []

    for direction, direction_frame in predictions.groupby("direction", sort=True):
        split_ids = sorted(direction_frame["split"].unique())
        panels = {}
        reference = None
        for split in split_ids:
            split_panels = {}
            for model, model_frame in direction_frame.loc[direction_frame["split"].eq(split)].groupby("model", sort=True):
                panel = model_frame.sort_values("sample_id").reset_index(drop=True)
                key = panel[["sample_id", "outcome_ibd", "bootstrap_cluster"]]
                if reference is None:
                    reference = key
                elif not reference.equals(key):
                    raise ValueError(f"{direction} split {split}: target panel differs across models")
                split_panels[model] = panel
            if set(MODELS).difference(split_panels):
                raise ValueError(f"{direction} split {split}: incomplete model set")
            panels[int(split)] = split_panels

        base = panels[split_ids[0]][MODELS[0]]
        y = base["outcome_ibd"].to_numpy(dtype=int)
        clusters = base["bootstrap_cluster"].to_numpy()
        n_clusters = int(pd.unique(clusters).size)
        rng = np.random.default_rng(args.seed + sum(ord(c) for c in direction))
        metric_draws = {model: {metric: [] for metric in METRICS} for model in MODELS}
        observed_split = {model: {metric: [] for metric in METRICS} for model in MODELS}

        for split in split_ids:
            for model in MODELS:
                panel = panels[split][model]
                yy = panel["outcome_ibd"].to_numpy(dtype=int)
                pp = panel["probability"].to_numpy(dtype=float)
                observed_split[model]["roc_auc"].append(weighted_auc(yy, pp, np.ones(len(yy))))
                observed_split[model]["average_precision"].append(average_precision_score(yy, pp))
                observed_split[model]["brier"].append(brier_score_loss(yy, pp))

        for replicate in range(args.n_bootstrap):
            w = bootstrap_weights(clusters, y, rng)
            per_split = {model: {metric: [] for metric in METRICS} for model in MODELS}
            for split in split_ids:
                for model in MODELS:
                    panel = panels[split][model]
                    yy = panel["outcome_ibd"].to_numpy(dtype=int)
                    pp = panel["probability"].to_numpy(dtype=float)
                    per_split[model]["roc_auc"].append(weighted_auc(yy, pp, w))
                    per_split[model]["average_precision"].append(
                        average_precision_score(yy, pp, sample_weight=w)
                    )
                    per_split[model]["brier"].append(
                        brier_score_loss(yy, pp, sample_weight=w)
                    )
            for model in MODELS:
                for metric in METRICS:
                    value = float(np.mean(per_split[model][metric]))
                    metric_draws[model][metric].append(value)
                    draws.append({
                        "direction": direction,
                        "model": model,
                        "replicate": replicate,
                        "metric": metric,
                        "value": value,
                        "n_esp_samples": int(len(base)),
                        "n_esp_participants": n_clusters,
                    })

        for model in MODELS:
            for metric in METRICS:
                values = np.asarray(metric_draws[model][metric], dtype=float)
                rows.append({
                    "sensitivity": "country_ESP",
                    "direction": direction,
                    "model": model,
                    "n_samples": int(len(base)),
                    "n_participants": n_clusters,
                    "n_control_samples": int((y == 0).sum()),
                    "n_ibd_samples": int((y == 1).sum()),
                    "source_split_mean": float(np.mean(observed_split[model][metric])),
                    "source_split_sd": float(np.std(observed_split[model][metric], ddof=1)),
                    "metric": metric,
                    "shared_weight_bootstrap_mean": float(values.mean()),
                    "ci_low": float(np.quantile(values, 0.025)),
                    "ci_high": float(np.quantile(values, 0.975)),
                    "ci_half_width": float((np.quantile(values, 0.975) - np.quantile(values, 0.025)) / 2),
                })

    pd.DataFrame(rows).to_csv(out / "esp_full_model_audit_summary.csv", index=False)
    pd.DataFrame(draws).to_csv(out / "esp_full_model_audit_draws.csv", index=False)
    manifest_out = {
        "protocol": "frozen_target_ESP_full_model_audit_shared_participant_bootstrap",
        "predictions_sha256": sha256(args.predictions),
        "manifest_sha256": sha256(args.manifest),
        "directions": list(TARGET_DIRECTIONS),
        "models": MODELS,
        "n_source_perturbations": 10,
        "n_bootstrap": args.n_bootstrap,
        "seed": args.seed,
        "target_country": "ESP",
        "n_esp_samples": int(predictions["sample_id"].nunique()),
        "n_esp_participants": int(target_meta["participant_id"].nunique()),
        "models_refit": False,
        "target_labels_used_for_refitting_or_selection": False,
    }
    (out / "esp_full_model_audit_manifest.json").write_text(
        json.dumps(manifest_out, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest_out, indent=2))


if __name__ == "__main__":
    main()
