#!/usr/bin/env python
"""Leakage-safe CPU baseline and fixed-target transfer benchmark."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedShuffleSplit

from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


def _metrics(y: pd.Series, p: np.ndarray) -> dict[str, float]:
    return {"roc_auc": float(roc_auc_score(y, p)), "average_precision": float(average_precision_score(y, p)), "brier": float(brier_score_loss(y, p))}


def _fit_models(x_train: pd.DataFrame, y_train: pd.Series, x_val: pd.DataFrame, y_val: pd.Series, seed: int):
    candidates = {
        "logistic": [LogisticRegression(C=c, class_weight="balanced", max_iter=3000, random_state=seed) for c in (0.1, 1.0, 10.0)],
        "random_forest": [RandomForestClassifier(n_estimators=300, max_depth=d, min_samples_leaf=2, class_weight="balanced", n_jobs=1, random_state=seed) for d in (None, 6)],
    }
    selected = {}
    for name, models in candidates.items():
        scored = []
        for model in models:
            model.fit(x_train, y_train)
            probability = model.predict_proba(x_val)[:, 1]
            scored.append((roc_auc_score(y_val, probability), model))
        selected[name] = max(scored, key=lambda pair: pair[0])[1]
    return selected


def _source_splits(manifest: pd.DataFrame, n_splits: int, seed: int, grouped: bool):
    if grouped:
        yield from ((sid, tr, va) for sid, tr, va in repeated_group_train_validation(manifest, n_splits, 0.2, seed))
        return
    strata = manifest["origin"].astype(str) + "__" + manifest["outcome_ibd"].astype(str)
    splitter = StratifiedShuffleSplit(n_splits=n_splits, test_size=0.2, random_state=seed)
    for sid, (tr, va) in enumerate(splitter.split(manifest, strata), 1):
        yield sid, tr, va


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-dir", type=Path, required=True)
    parser.add_argument("--preprocess-dir", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-splits", type=int, default=5)
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    contract = args.contract_dir.resolve()
    hmp_manifest = pd.read_csv(contract / "cohort_manifest_hmp2.csv")
    prism_manifest = pd.read_csv(contract / "cohort_manifest_prism.csv")
    shared = load_shared_feature_list(contract)
    hmp = load_canonical_microbiome(args.legacy_root / "manual_abundance.csv", hmp_manifest, shared)
    prism = load_canonical_microbiome(args.legacy_root / "prism_abundance_clean.csv", prism_manifest, shared)
    hmp_manifest = hmp_manifest.set_index("sample_id", drop=False).loc[hmp.index].reset_index(drop=True)
    prism_manifest = prism_manifest.set_index("sample_id", drop=False).loc[prism.index].reset_index(drop=True)
    cohorts = [("HMP2_to_PRISM", hmp, hmp_manifest, prism, prism_manifest, True), ("PRISM_to_HMP2", prism, prism_manifest, hmp, hmp_manifest, False)]
    rows = []
    predictions = []
    for direction, source, source_manifest, target, target_manifest, grouped in cohorts:
        for split_id, train_idx, val_idx in _source_splits(source_manifest, args.n_splits, 20260815, grouped):
            train_ids = source_manifest.iloc[train_idx]["sample_id"].tolist(); val_ids = source_manifest.iloc[val_idx]["sample_id"].tolist()
            pipeline = CLRFeaturePipeline(min_prevalence=0.10, top_k=100, random_state=20260815 + split_id)
            pipeline.fit(source.loc[train_ids], source_manifest.iloc[train_idx]["outcome_ibd"])
            x_train = pipeline.transform(source.loc[train_ids]); x_val = pipeline.transform(source.loc[val_ids]); x_target = pipeline.transform(target)
            models = _fit_models(x_train, source_manifest.iloc[train_idx]["outcome_ibd"], x_val, source_manifest.iloc[val_idx]["outcome_ibd"], 20260815 + split_id)
            for model_name, model in models.items():
                for evaluation, x_eval, y_eval in (("source_validation", x_val, source_manifest.iloc[val_idx]["outcome_ibd"]), ("external_target", x_target, target_manifest["outcome_ibd"])):
                    probability = model.predict_proba(x_eval)[:, 1]
                    rows.append({"direction": direction, "split": split_id, "model": model_name, "evaluation": evaluation, **_metrics(y_eval, probability), "n_train": len(train_idx), "n_validation": len(val_idx), "n_target": len(target), "selected_features": pipeline.audit_.n_selected_features, "target_origin_counts": json.dumps(target_manifest.get("origin", pd.Series(dtype=str)).value_counts().to_dict(), sort_keys=True)})
                    if evaluation == "external_target":
                        for sample_id, outcome, score, cluster, origin in zip(target.index, y_eval, probability, target_manifest["participant_id"], target_manifest.get("origin", pd.Series("not_available", index=target_manifest.index))):
                            predictions.append({"direction": direction, "split": split_id, "model": model_name, "sample_id": sample_id, "outcome_ibd": int(outcome), "probability": float(score), "bootstrap_cluster": str(cluster), "origin": origin})
    results = pd.DataFrame(rows)
    results.to_csv(out / "cpu_baseline_metrics.csv", index=False)
    pd.DataFrame(predictions).to_csv(out / "cpu_baseline_external_predictions.csv", index=False)
    aggregate = results[results.evaluation == "external_target"].groupby(["direction", "model"])[["roc_auc", "average_precision", "brier"]].agg(["mean", "std"]).round(6)
    aggregate.columns = [f"{metric}_{stat}" for metric, stat in aggregate.columns]
    summary = {"n_rows": len(results), "directions": [x[0] for x in cohorts], "models": ["logistic", "random_forest"], "n_splits": args.n_splits, "shared_features": len(shared), "external_target_summary": aggregate.reset_index().to_dict(orient="records")}
    (out / "cpu_baseline_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
