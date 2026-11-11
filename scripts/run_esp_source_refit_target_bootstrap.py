from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(r"D:/IBD_microbiome")
PRED = ROOT / "ibd_framework_scaffold/artifacts/runs/esp_source_refit_20261030/esp_source_refit_predictions.csv"
FULL = ROOT / "ibd_framework_scaffold/manuscript/source_data/standard_model_benchmark_predictions.csv"
OUT = ROOT / "ibd_framework_scaffold/artifacts/runs/esp_source_refit_target_bootstrap_20261031"


def weighted_metric(y: np.ndarray, p: np.ndarray, weights: np.ndarray, metric: str) -> float:
    if metric == "roc_auc":
        return float(roc_auc_score(y, p, sample_weight=weights))
    if metric == "average_precision":
        return float(average_precision_score(y, p, sample_weight=weights))
    if metric == "brier":
        return float(np.average((p - y) ** 2, weights=weights))
    raise ValueError(metric)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    esp = pd.read_csv(PRED)
    full = pd.read_csv(FULL)
    model_order = ["logistic_l2", "random_forest"]
    direction_pairs = [
        ("ESP_NielsenHB_2014_to_HMP2", "NielsenHB_2014_to_HMP2"),
        ("ESP_NielsenHB_2014_to_PRISM", "NielsenHB_2014_to_PRISM"),
    ]
    metrics = ["roc_auc", "average_precision", "brier"]
    rng = np.random.default_rng(20261031)
    n_bootstrap = 500
    summary_rows = []
    paired_rows = []
    draw_rows = []
    for esp_direction, full_direction in direction_pairs:
        for model in model_order:
            e = esp[(esp.direction == esp_direction) & (esp.model == model)].copy()
            f = full[(full.direction == full_direction) & (full.model == model)].copy()
            if e.empty or f.empty:
                raise ValueError((esp_direction, full_direction, model))
            e = e.sort_values(["split", "sample_id"])
            f = f.sort_values(["split", "sample_id"])
            if set(e.sample_id) != set(f.sample_id):
                raise ValueError("Target sample IDs do not match between ESP refit and full Nielsen source")
            ref = e.drop_duplicates("sample_id")[["sample_id", "outcome_ibd", "bootstrap_cluster"]].copy()
            clusters = ref.bootstrap_cluster.drop_duplicates().tolist()
            cluster_to_index = {cluster: i for i, cluster in enumerate(clusters)}
            # Ensure source perturbation rows align to the same target unit.
            e_by_split = {int(s): x.set_index("sample_id").loc[ref.sample_id] for s, x in e.groupby("split")}
            f_by_split = {int(s): x.set_index("sample_id").loc[ref.sample_id] for s, x in f.groupby("split")}
            split_ids = sorted(set(e_by_split) & set(f_by_split))
            if len(split_ids) != 10:
                raise ValueError(f"Expected ten source perturbations, found {split_ids}")
            y = ref.outcome_ibd.to_numpy(dtype=int)
            cluster_ids = ref.bootstrap_cluster.to_numpy()
            split_values = {metric: [] for metric in metrics}
            for s in split_ids:
                for metric in metrics:
                    split_values[metric].append(weighted_metric(y, e_by_split[s].probability.to_numpy(), np.ones(len(y)), metric))
            source_sd = {metric: float(np.std(values, ddof=1)) for metric, values in split_values.items()}
            draws = {metric: [] for metric in metrics}
            paired_draws = {metric: [] for metric in metrics}
            for b in range(n_bootstrap):
                sampled = rng.choice(np.arange(len(clusters)), size=len(clusters), replace=True)
                weights = np.bincount(sampled, minlength=len(clusters)).astype(float)
                row_weights = np.array([weights[cluster_to_index[c]] for c in cluster_ids], dtype=float)
                for metric in metrics:
                    e_values = []
                    f_values = []
                    for s in split_ids:
                        e_values.append(weighted_metric(y, e_by_split[s].probability.to_numpy(), row_weights, metric))
                        f_values.append(weighted_metric(y, f_by_split[s].probability.to_numpy(), row_weights, metric))
                    e_mean = float(np.mean(e_values))
                    f_mean = float(np.mean(f_values))
                    draws[metric].append(e_mean)
                    paired_draws[metric].append(e_mean - f_mean)
                    draw_rows.append({
                        "direction": esp_direction,
                        "model": model,
                        "metric": metric,
                        "replicate": b,
                        "target_bootstrap_mean": e_mean,
                        "paired_delta": e_mean - f_mean,
                        "n_target_clusters": len(clusters),
                    })
            for metric in metrics:
                values = np.asarray(draws[metric])
                paired = np.asarray(paired_draws[metric])
                summary_rows.append({
                    "direction": esp_direction,
                    "target_label": full_direction,
                    "model": model,
                    "metric": metric,
                    "n_source_perturbations": len(split_ids),
                    "n_target_clusters": len(clusters),
                    "source_split_mean": float(np.mean(split_values[metric])),
                    "source_split_sd": source_sd[metric],
                    "target_bootstrap_mean": float(np.mean(values)),
                    "target_ci_low": float(np.quantile(values, 0.025)),
                    "target_ci_high": float(np.quantile(values, 0.975)),
                    "estimand": "ESP-source fixed-contract refit; perturbation-averaged target metric with shared target participant bootstrap",
                })
                paired_rows.append({
                    "direction": esp_direction,
                    "target_label": full_direction,
                    "model": model,
                    "metric": metric,
                    "reference": "full_Nielsen_source",
                    "n_target_clusters": len(clusters),
                    "paired_delta_mean": float(np.mean(paired)),
                    "paired_ci_low": float(np.quantile(paired, 0.025)),
                    "paired_ci_high": float(np.quantile(paired, 0.975)),
                    "source_selection_sd": source_sd[metric],
                    "estimand": "ESP-source minus full-Nielsen-source predictions on the same locked target participants and shared bootstrap weights",
                })
    pd.DataFrame(summary_rows).to_csv(OUT / "esp_source_refit_target_bootstrap_summary.csv", index=False)
    pd.DataFrame(paired_rows).to_csv(OUT / "esp_source_refit_target_paired_summary.csv", index=False)
    pd.DataFrame(draw_rows).to_csv(OUT / "esp_source_refit_target_bootstrap_draws.csv", index=False)
    manifest = {
        "protocol": "post_fit_ESP_source_refit_target_participant_bootstrap",
        "predictions": str(PRED),
        "full_source_predictions": str(FULL),
        "n_bootstrap": n_bootstrap,
        "seed": 20261031,
        "n_source_perturbations": 10,
        "target_weights_shared_across_source_perturbations": True,
        "models_refit": False,
        "target_labels_used_for_refitting_or_selection": False,
    }
    (OUT / "esp_source_refit_target_bootstrap_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
