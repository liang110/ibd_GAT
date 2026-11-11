#!/usr/bin/env python
"""Measure known-truth recovery across deterministic semi-synthetic seeds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

from ibd_microbiome_framework.evaluation.synthetic import make_portability_benchmark
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


def evaluate_seed(background: pd.DataFrame, seed: int) -> dict:
    source_x, source_y, target_x, target_y, truth = make_portability_benchmark(background, random_state=seed)
    pipeline = CLRFeaturePipeline(min_prevalence=0.10, top_k=20, random_state=seed).fit(source_x, source_y)
    model = LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000, random_state=seed)
    model.fit(pipeline.transform(source_x), source_y)
    probability = model.predict_proba(pipeline.transform(target_x))[:, 1]
    selected = set(pipeline.selected_features_)
    disease_recovered = len(selected.intersection(truth.disease_features))
    batch_selected = len(selected.intersection(truth.batch_features))
    return {
        "seed": seed,
        "target_roc_auc": float(roc_auc_score(target_y, probability)),
        "target_average_precision": float(average_precision_score(target_y, probability)),
        "n_disease_features_recovered": disease_recovered,
        "n_batch_features_selected": batch_selected,
        "roc_auc_pass": bool(roc_auc_score(target_y, probability) >= 0.70),
        "disease_recovery_pass": bool(disease_recovered >= 3),
        "all_acceptance_pass": bool(roc_auc_score(target_y, probability) >= 0.70 and disease_recovered >= 3),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(range(20260815, 20260825)))
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    background = pd.read_csv(args.background, index_col=0)
    metrics = pd.DataFrame([evaluate_seed(background, seed) for seed in args.seeds])
    metrics.to_csv(out / "synthetic_sweep_metrics.csv", index=False)
    summary = {
        "seeds": args.seeds,
        "n_seeds": len(args.seeds),
        "metrics_mean": metrics[["target_roc_auc", "target_average_precision", "n_disease_features_recovered", "n_batch_features_selected"]].mean().round(6).to_dict(),
        "metrics_std": metrics[["target_roc_auc", "target_average_precision", "n_disease_features_recovered", "n_batch_features_selected"]].std().round(6).to_dict(),
        "roc_auc_pass_rate": float(metrics.roc_auc_pass.mean()),
        "disease_recovery_pass_rate": float(metrics.disease_recovery_pass.mean()),
        "all_acceptance_pass_rate": float(metrics.all_acceptance_pass.mean()),
    }
    (out / "synthetic_sweep_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
