#!/usr/bin/env python
"""Run descriptive CD/UC-versus-control transfer without pooling PRISM origins."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation, repeated_stratified_train_validation
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


def _metrics(y: pd.Series, probability: np.ndarray) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(y, probability)),
        "average_precision": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }


def _subset(manifest: pd.DataFrame, matrix: pd.DataFrame, subtype: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    keep = manifest.outcome_ibd.eq(0) | manifest.diagnosis_subtype.eq(subtype)
    subset_manifest = manifest.loc[keep].copy().reset_index(drop=True)
    if subset_manifest.outcome_ibd.nunique() != 2:
        raise ValueError(f"{subtype} subset does not contain both outcomes")
    return subset_manifest, matrix.loc[subset_manifest.sample_id]


def _select_elastic_net(x_train, y_train, x_validation, y_validation, seed: int):
    candidates = []
    for c in (0.1, 1.0):
        for l1_ratio in (0.1, 0.5, 0.9):
            model = LogisticRegression(
                C=c,
                penalty="elasticnet",
                solver="saga",
                l1_ratio=l1_ratio,
                class_weight="balanced",
                max_iter=5000,
                random_state=seed,
            ).fit(x_train, y_train)
            candidates.append((roc_auc_score(y_validation, model.predict_proba(x_validation)[:, 1]), model, c, l1_ratio))
    return max(candidates, key=lambda item: item[0])


def _splits(manifest: pd.DataFrame, grouped: bool, n_splits: int):
    if grouped:
        return repeated_group_train_validation(manifest, n_splits, 0.2, 20260816)
    return repeated_stratified_train_validation(manifest, n_splits, 0.2, 20260816)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-dir", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--case-subtype", choices=("CD", "UC"), required=True)
    parser.add_argument("--n-splits", type=int, default=5)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    hmp_manifest = pd.read_csv(args.contract_dir / "cohort_manifest_hmp2.csv")
    prism_manifest = pd.read_csv(args.contract_dir / "cohort_manifest_prism.csv")
    shared = load_shared_feature_list(args.contract_dir)
    hmp = load_canonical_microbiome(args.legacy_root / "manual_abundance.csv", hmp_manifest, shared)
    prism = load_canonical_microbiome(args.legacy_root / "prism_abundance_clean.csv", prism_manifest, shared)
    hmp_manifest = hmp_manifest.set_index("sample_id", drop=False).loc[hmp.index].reset_index(drop=True)
    prism_manifest = prism_manifest.set_index("sample_id", drop=False).loc[prism.index].reset_index(drop=True)
    hmp_manifest, hmp = _subset(hmp_manifest, hmp, args.case_subtype)
    prism_manifest, prism = _subset(prism_manifest, prism, args.case_subtype)
    prism_origins = {
        origin: (prism.loc[group.sample_id], group.reset_index(drop=True))
        for origin, group in prism_manifest.groupby("origin", sort=True)
    }
    cohorts = []
    for origin, (origin_matrix, origin_manifest) in prism_origins.items():
        cohorts.extend(
            [
                (f"HMP2_to_PRISM_{origin}", hmp, hmp_manifest, origin_matrix, origin_manifest, True),
                (f"PRISM_{origin}_to_HMP2", origin_matrix, origin_manifest, hmp, hmp_manifest, False),
            ]
        )
    metrics, predictions = [], []
    for direction, source, source_manifest, target, target_manifest, grouped in cohorts:
        for split_id, train_idx, validation_idx in _splits(source_manifest, grouped, args.n_splits):
            train_manifest = source_manifest.iloc[train_idx]
            validation_manifest = source_manifest.iloc[validation_idx]
            pipeline = CLRFeaturePipeline(min_prevalence=0.10, top_k=100, random_state=20260816 + split_id).fit(
                source.loc[train_manifest.sample_id], train_manifest.outcome_ibd
            )
            x_train = pipeline.transform(source.loc[train_manifest.sample_id])
            x_validation = pipeline.transform(source.loc[validation_manifest.sample_id])
            x_target = pipeline.transform(target)
            _, model, selected_c, selected_l1_ratio = _select_elastic_net(
                x_train, train_manifest.outcome_ibd, x_validation, validation_manifest.outcome_ibd, 20260816 + split_id
            )
            for evaluation, y_eval, probability in (
                ("source_validation", validation_manifest.outcome_ibd, model.predict_proba(x_validation)[:, 1]),
                ("external_target", target_manifest.outcome_ibd, model.predict_proba(x_target)[:, 1]),
            ):
                metrics.append({
                    "case_subtype": args.case_subtype,
                    "direction": direction,
                    "split": split_id,
                    "model": "elastic_net",
                    "evaluation": evaluation,
                    "selected_c": selected_c,
                    "selected_l1_ratio": selected_l1_ratio,
                    "n_train_samples": len(train_manifest),
                    "n_train_participants": train_manifest.participant_id.nunique(),
                    "n_validation_samples": len(validation_manifest),
                    "n_validation_participants": validation_manifest.participant_id.nunique(),
                    "n_target_samples": len(target_manifest),
                    "n_target_participants": target_manifest.participant_id.nunique(),
                    "selected_features": pipeline.audit_.n_selected_features,
                    **_metrics(y_eval, probability),
                })
                if evaluation == "external_target":
                    for sample_id, outcome, score, cluster, origin in zip(
                        target.index, y_eval, probability, target_manifest.participant_id, target_manifest.origin
                    ):
                        predictions.append({
                            "case_subtype": args.case_subtype,
                            "direction": direction,
                            "split": split_id,
                            "model": "elastic_net",
                            "sample_id": sample_id,
                            "outcome_ibd": int(outcome),
                            "probability": float(score),
                            "bootstrap_cluster": str(cluster),
                            "origin": origin,
                        })
    metric_frame = pd.DataFrame(metrics)
    metric_frame.to_csv(out / "subtype_transfer_metrics.csv", index=False)
    pd.DataFrame(predictions).to_csv(out / "subtype_transfer_predictions.csv", index=False)
    summary = metric_frame.loc[metric_frame.evaluation.eq("external_target")].groupby("direction")[
        ["roc_auc", "average_precision", "brier"]
    ].agg(["mean", "std"]).round(6)
    summary.columns = [f"{metric}_{stat}" for metric, stat in summary.columns]
    report = {
        "case_subtype": args.case_subtype,
        "n_splits": args.n_splits,
        "shared_features": len(shared),
        "external_target_summary": summary.reset_index().to_dict(orient="records"),
    }
    (out / "subtype_transfer_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
