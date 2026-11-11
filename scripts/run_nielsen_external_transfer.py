#!/usr/bin/env python
"""Evaluate source-only IBD classifiers on the locked NielsenHB_2014 target."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from xgboost import XGBClassifier

from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import (
    repeated_group_train_validation,
    repeated_stratified_train_validation,
)
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


EXPERIMENT_SEED = 20260816


def _metrics(y: pd.Series, probability: np.ndarray) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(y, probability)),
        "average_precision": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }


def _candidates(seed: int) -> dict[str, list]:
    return {
        "logistic_l2": [
            LogisticRegression(C=c, class_weight="balanced", max_iter=5000, random_state=seed)
            for c in (0.1, 1.0, 10.0)
        ],
        "elastic_net": [
            LogisticRegression(
                C=c,
                penalty="elasticnet",
                solver="saga",
                l1_ratio=l1_ratio,
                class_weight="balanced",
                max_iter=5000,
                random_state=seed,
            )
            for c in (0.1, 1.0)
            for l1_ratio in (0.1, 0.5, 0.9)
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
        "xgboost": [
            XGBClassifier(
                n_estimators=400,
                max_depth=depth,
                learning_rate=learning_rate,
                subsample=0.8,
                colsample_bytree=0.8,
                objective="binary:logistic",
                eval_metric="auc",
                early_stopping_rounds=30,
                n_jobs=1,
                random_state=seed,
                tree_method="hist",
            )
            for depth in (2, 3)
            for learning_rate in (0.03, 0.1)
        ],
    }


def _fit_select(x_train, y_train, x_validation, y_validation, seed):
    selected = {}
    for name, models in _candidates(seed).items():
        scored = []
        for model in models:
            if name == "xgboost":
                model.fit(x_train, y_train, eval_set=[(x_validation, y_validation)], verbose=False)
            else:
                model.fit(x_train, y_train)
            scored.append((roc_auc_score(y_validation, model.predict_proba(x_validation)[:, 1]), model))
        selected[name] = max(scored, key=lambda item: item[0])[1]
    return selected


def _source_splits(manifest: pd.DataFrame, grouped: bool, n_splits: int):
    if grouped:
        return repeated_group_train_validation(manifest, n_splits, 0.2, EXPERIMENT_SEED)
    return repeated_stratified_train_validation(manifest, n_splits, 0.2, EXPERIMENT_SEED)


def _load_manifest(path: Path, matrix: pd.DataFrame) -> pd.DataFrame:
    manifest = pd.read_csv(path).set_index("sample_id", drop=False)
    return manifest.loc[matrix.index].reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-dir", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--nielsen-projection", type=Path, required=True)
    parser.add_argument("--nielsen-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-splits", type=int, default=10)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    shared = load_shared_feature_list(args.contract_dir)
    hmp_manifest = pd.read_csv(args.contract_dir / "cohort_manifest_hmp2.csv")
    prism_manifest = pd.read_csv(args.contract_dir / "cohort_manifest_prism.csv")
    hmp = load_canonical_microbiome(args.legacy_root / "manual_abundance.csv", hmp_manifest, shared)
    prism = load_canonical_microbiome(args.legacy_root / "prism_abundance_clean.csv", prism_manifest, shared)
    hmp_manifest = _load_manifest(args.contract_dir / "cohort_manifest_hmp2.csv", hmp)
    prism_manifest = _load_manifest(args.contract_dir / "cohort_manifest_prism.csv", prism)
    target_manifest = pd.read_csv(args.nielsen_manifest)
    target = load_canonical_microbiome(args.nielsen_projection, target_manifest, shared)
    target_manifest = _load_manifest(args.nielsen_manifest, target)

    cohorts = [
        ("HMP2_to_NielsenHB_2014", hmp, hmp_manifest, True),
        ("PRISM_to_NielsenHB_2014", prism, prism_manifest, False),
    ]
    metrics: list[dict[str, object]] = []
    predictions: list[dict[str, object]] = []
    for direction, source, source_manifest, grouped in cohorts:
        for split_id, train_idx, validation_idx in _source_splits(source_manifest, grouped, args.n_splits):
            train_ids = source_manifest.iloc[train_idx].sample_id.tolist()
            validation_ids = source_manifest.iloc[validation_idx].sample_id.tolist()
            pipeline = CLRFeaturePipeline(
                min_prevalence=0.10, top_k=100, random_state=EXPERIMENT_SEED + split_id
            ).fit(source.loc[train_ids], source_manifest.iloc[train_idx].outcome_ibd)
            x_train = pipeline.transform(source.loc[train_ids])
            x_validation = pipeline.transform(source.loc[validation_ids])
            x_target = pipeline.transform(target)
            fitted = _fit_select(
                x_train,
                source_manifest.iloc[train_idx].outcome_ibd,
                x_validation,
                source_manifest.iloc[validation_idx].outcome_ibd,
                EXPERIMENT_SEED + split_id,
            )
            for model_name, model in fitted.items():
                for evaluation, x_eval, y_eval in (
                    ("source_validation", x_validation, source_manifest.iloc[validation_idx].outcome_ibd),
                    ("locked_external_target", x_target, target_manifest.outcome_ibd),
                ):
                    probability = model.predict_proba(x_eval)[:, 1]
                    metrics.append(
                        {
                            "direction": direction,
                            "split": split_id,
                            "model": model_name,
                            "evaluation": evaluation,
                            "n_train": len(train_idx),
                            "n_validation": len(validation_idx),
                            "n_target": len(target_manifest),
                            "selected_features": pipeline.audit_.n_selected_features,
                            **_metrics(y_eval, probability),
                        }
                    )
                    if evaluation == "locked_external_target":
                        for sample_id, outcome, score, participant_id, country in zip(
                            target.index,
                            y_eval,
                            probability,
                            target_manifest.participant_id,
                            target_manifest.origin,
                        ):
                            predictions.append(
                                {
                                    "direction": direction,
                                    "split": split_id,
                                    "model": model_name,
                                    "sample_id": sample_id,
                                    "outcome_ibd": int(outcome),
                                    "probability": float(score),
                                    "bootstrap_cluster": str(participant_id),
                                    "origin": str(country),
                                }
                            )

    metric_frame = pd.DataFrame(metrics)
    metric_frame.to_csv(out / "nielsen_external_transfer_metrics.csv", index=False)
    pd.DataFrame(predictions).to_csv(out / "nielsen_external_transfer_predictions.csv", index=False)
    external = metric_frame.loc[metric_frame.evaluation.eq("locked_external_target")]
    aggregate = external.groupby(["direction", "model"])[
        ["roc_auc", "average_precision", "brier"]
    ].agg(["mean", "std"])
    aggregate.columns = [f"{metric}_{stat}" for metric, stat in aggregate.columns]
    summary = {
        "protocol": "source_only_selection_locked_external_target",
        "experiment_seed": EXPERIMENT_SEED,
        "n_splits": args.n_splits,
        "shared_features": len(shared),
        "target_labels_used_for_model_selection": False,
        "external_target_summary": aggregate.reset_index().round(6).to_dict(orient="records"),
    }
    (out / "nielsen_external_transfer_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
