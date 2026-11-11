#!/usr/bin/env python
"""Run a known-truth test for the leakage-safe shared-space workflow."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

from ibd_microbiome_framework.evaluation.synthetic import make_portability_benchmark
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--random-state", type=int, default=20260815)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    background = pd.read_csv(args.background, index_col=0)
    source_x, source_y, target_x, target_y, truth = make_portability_benchmark(
        background, random_state=args.random_state
    )
    pipeline = CLRFeaturePipeline(min_prevalence=0.10, top_k=20, random_state=args.random_state)
    pipeline.fit(source_x, source_y)
    model = LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000, random_state=args.random_state)
    model.fit(pipeline.transform(source_x), source_y)
    target_probability = model.predict_proba(pipeline.transform(target_x))[:, 1]
    selected = set(pipeline.selected_features_)
    disease_hits = sorted(selected.intersection(truth.disease_features))
    batch_hits = sorted(selected.intersection(truth.batch_features))
    metrics = {
        "target_roc_auc": float(roc_auc_score(target_y, target_probability)),
        "target_average_precision": float(average_precision_score(target_y, target_probability)),
        "disease_features_injected": list(truth.disease_features),
        "batch_features_injected": list(truth.batch_features),
        "disease_features_recovered": disease_hits,
        "batch_features_selected": batch_hits,
        "n_disease_features_recovered": len(disease_hits),
        "n_batch_features_selected": len(batch_hits),
        "n_selected_features": len(selected),
        "preprocessing_audit": pipeline.audit_.__dict__,
        "acceptance": {
            "target_roc_auc_at_least_0_70": bool(roc_auc_score(target_y, target_probability) >= 0.70),
            "at_least_half_disease_features_recovered": bool(len(disease_hits) >= 3),
        },
    }
    pd.DataFrame(
        {"feature": pipeline.selected_features_, "feature_class": [
            "disease_truth" if x in truth.disease_features else "source_batch_truth" if x in truth.batch_features else "background"
            for x in pipeline.selected_features_
        ]}
    ).to_csv(output / "selected_feature_truth_audit.csv", index=False)
    (output / "synthetic_validation_summary.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
