#!/usr/bin/env python
"""Evaluate training-domain stable feature selection across synthetic seeds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

from ibd_microbiome_framework.evaluation.synthetic import make_portability_benchmark
from ibd_microbiome_framework.preprocessing.compositional import DomainStableCLRFeaturePipeline


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(range(20260815, 20260825)))
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    background = pd.read_csv(args.background, index_col=0)
    rows = []
    for seed in args.seeds:
        source_x, source_y, target_x, target_y, truth, source_domains, _ = make_portability_benchmark(
            background,
            random_state=seed,
            return_domains=True,
            min_truth_prevalence=0.30,
        )
        pipeline = DomainStableCLRFeaturePipeline(min_prevalence=0.10, top_k=20, random_state=seed).fit(source_x, source_y, source_domains)
        model = LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000, random_state=seed)
        model.fit(pipeline.transform(source_x), source_y)
        probability = model.predict_proba(pipeline.transform(target_x))[:, 1]
        selected = set(pipeline.selected_features_)
        disease_recovered = len(selected.intersection(truth.disease_features))
        batch_selected = len(selected.intersection(truth.batch_features))
        auc = float(roc_auc_score(target_y, probability))
        rows.append({"seed": seed, "target_roc_auc": auc, "target_average_precision": float(average_precision_score(target_y, probability)), "n_disease_features_recovered": disease_recovered, "n_batch_features_selected": batch_selected, "n_stable_candidates": pipeline.stable_feature_count_, "all_acceptance_pass": bool(auc >= 0.70 and disease_recovered >= 3)})
    metrics = pd.DataFrame(rows); metrics.to_csv(out / "domain_stable_sweep_metrics.csv", index=False)
    summary = {"benchmark_version": "v2_prevalent_truth_features", "min_truth_prevalence": 0.30, "n_seeds": len(metrics), "target_roc_auc_mean": float(metrics.target_roc_auc.mean()), "target_roc_auc_std": float(metrics.target_roc_auc.std()), "target_average_precision_mean": float(metrics.target_average_precision.mean()), "disease_features_recovered_mean": float(metrics.n_disease_features_recovered.mean()), "batch_features_selected_mean": float(metrics.n_batch_features_selected.mean()), "all_acceptance_pass_rate": float(metrics.all_acceptance_pass.mean())}
    (out / "domain_stable_sweep_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
