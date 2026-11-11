#!/usr/bin/env python
"""Run clean PRISM discovery/validation bidirectional transfer on CPU."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import repeated_stratified_train_validation
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


def _metrics(y: pd.Series, probability: np.ndarray) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(y, probability)),
        "average_precision": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }


def _select_models(x_train, y_train, x_validation, y_validation, seed):
    candidates = {
        "logistic": [
            LogisticRegression(C=c, class_weight="balanced", max_iter=3000, random_state=seed)
            for c in (0.1, 1.0, 10.0)
        ],
        "random_forest": [
            RandomForestClassifier(
                n_estimators=300,
                max_depth=depth,
                min_samples_leaf=2,
                class_weight="balanced",
                n_jobs=1,
                random_state=seed,
            )
            for depth in (None, 6)
        ],
    }
    selected = {}
    for name, models in candidates.items():
        scored = []
        for model in models:
            model.fit(x_train, y_train)
            score = roc_auc_score(y_validation, model.predict_proba(x_validation)[:, 1])
            scored.append((score, model))
        selected[name] = max(scored, key=lambda item: item[0])[1]
    return selected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-dir", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-splits", type=int, default=5)
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.contract_dir / "cohort_manifest_prism.csv")
    shared = load_shared_feature_list(args.contract_dir)
    matrix = load_canonical_microbiome(
        args.legacy_root / "prism_abundance_clean.csv", manifest, shared
    )
    manifest = manifest.set_index("sample_id", drop=False).loc[matrix.index].reset_index(drop=True)
    directions = [("discovery", "validation"), ("validation", "discovery")]
    metric_rows, predictions = [], []
    for source_origin, target_origin in directions:
        source_manifest = manifest.loc[manifest.origin == source_origin].reset_index(drop=True)
        target_manifest = manifest.loc[manifest.origin == target_origin].reset_index(drop=True)
        source = matrix.loc[source_manifest.sample_id]
        target = matrix.loc[target_manifest.sample_id]
        direction = f"PRISM_{source_origin}_to_{target_origin}"
        for split_id, train_idx, validation_idx in repeated_stratified_train_validation(
            source_manifest, args.n_splits, 0.2, 20260816
        ):
            train_ids = source_manifest.iloc[train_idx].sample_id.tolist()
            validation_ids = source_manifest.iloc[validation_idx].sample_id.tolist()
            pipeline = CLRFeaturePipeline(
                min_prevalence=0.10, top_k=100, random_state=20260816 + split_id
            ).fit(source.loc[train_ids], source_manifest.iloc[train_idx].outcome_ibd)
            x_train = pipeline.transform(source.loc[train_ids])
            x_validation = pipeline.transform(source.loc[validation_ids])
            x_target = pipeline.transform(target)
            models = _select_models(
                x_train,
                source_manifest.iloc[train_idx].outcome_ibd,
                x_validation,
                source_manifest.iloc[validation_idx].outcome_ibd,
                20260816 + split_id,
            )
            for model_name, model in models.items():
                validation_probability = model.predict_proba(x_validation)[:, 1]
                target_probability = model.predict_proba(x_target)[:, 1]
                metric_rows.append({
                    "direction": direction,
                    "split": split_id,
                    "model": model_name,
                    "evaluation": "source_validation",
                    **_metrics(source_manifest.iloc[validation_idx].outcome_ibd, validation_probability),
                    "n_train": len(train_idx),
                    "n_validation": len(validation_idx),
                    "n_target": len(target_manifest),
                    "selected_features": pipeline.audit_.n_selected_features,
                })
                metric_rows.append({
                    "direction": direction,
                    "split": split_id,
                    "model": model_name,
                    "evaluation": "external_target",
                    **_metrics(target_manifest.outcome_ibd, target_probability),
                    "n_train": len(train_idx),
                    "n_validation": len(validation_idx),
                    "n_target": len(target_manifest),
                    "selected_features": pipeline.audit_.n_selected_features,
                })
                for sample_id, outcome, probability in zip(
                    target_manifest.sample_id, target_manifest.outcome_ibd, target_probability
                ):
                    predictions.append({
                        "direction": direction,
                        "split": split_id,
                        "model": model_name,
                        "sample_id": sample_id,
                        "outcome_ibd": int(outcome),
                        "probability": float(probability),
                        "bootstrap_cluster": sample_id,
                        "origin": target_origin,
                    })
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(out / "prism_origin_metrics.csv", index=False)
    pd.DataFrame(predictions).to_csv(out / "prism_origin_external_predictions.csv", index=False)
    aggregate = metrics.loc[metrics.evaluation == "external_target"].groupby(
        ["direction", "model"]
    )[["roc_auc", "average_precision", "brier"]].agg(["mean", "std"]).round(6)
    aggregate.columns = [f"{metric}_{stat}" for metric, stat in aggregate.columns]
    class_counts = manifest.groupby(["origin", "outcome_ibd"]).size().rename("n").reset_index()
    summary = {
        "n_splits": args.n_splits,
        "class_counts": class_counts.to_dict(orient="records"),
        "external_target_summary": aggregate.reset_index().to_dict(orient="records"),
    }
    (out / "prism_origin_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
