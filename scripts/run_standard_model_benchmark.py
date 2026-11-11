#!/usr/bin/env python
"""Run the pre-registered source-only standard model family benchmark."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import (
    repeated_group_train_validation,
    repeated_stratified_train_validation,
)
from ibd_microbiome_framework.evaluation.standard_benchmark import (
    fit_final_source_models,
    select_standard_models,
)


EXPERIMENT_SEED = 20260816


def _metrics(y: pd.Series, probability: np.ndarray) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(y, probability)),
        "average_precision": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }


def _load_manifest(path: Path, matrix: pd.DataFrame) -> pd.DataFrame:
    manifest = pd.read_csv(path).set_index("sample_id", drop=False)
    return manifest.loc[matrix.index].reset_index(drop=True)


def _source_splits(manifest: pd.DataFrame, grouped: bool, n_splits: int):
    if grouped:
        return repeated_group_train_validation(manifest, n_splits, 0.20, EXPERIMENT_SEED)
    return repeated_stratified_train_validation(manifest, n_splits, 0.20, EXPERIMENT_SEED)


def _direction_records(
    direction: str,
    source: pd.DataFrame,
    source_manifest: pd.DataFrame,
    target: pd.DataFrame,
    target_manifest: pd.DataFrame,
    grouped: bool,
    n_splits: int,
) -> tuple[list[dict], list[dict], list[dict]]:
    metrics: list[dict] = []
    predictions: list[dict] = []
    tuning: list[dict] = []
    for split_id, train_idx, validation_idx in _source_splits(source_manifest, grouped, n_splits):
        train = source.iloc[train_idx]
        train_y = source_manifest.iloc[train_idx].outcome_ibd.to_numpy()
        groups = source_manifest.iloc[train_idx].participant_id.to_numpy() if grouped else None
        selected = select_standard_models(
            train,
            train_y,
            groups=groups,
            seed=EXPERIMENT_SEED + split_id,
            n_inner_splits=3,
        )
        # Selection only sees the outer-source training partition. The final
        # model can use all source samples, but never target data or labels.
        pipeline, fitted = fit_final_source_models(
            source,
            source_manifest.outcome_ibd.to_numpy(),
            selected,
            seed=EXPERIMENT_SEED + split_id,
        )
        x_target = pipeline.transform(target)
        for family, model in fitted.items():
            selection = selected[family]
            probability = model.predict_proba(x_target)[:, 1]
            metrics.append(
                {
                    "direction": direction,
                    "split": split_id,
                    "model": family,
                    "evaluation": "locked_external_target",
                    "n_source": len(source),
                    "n_selection_train": len(train_idx),
                    "n_selection_validation": len(validation_idx),
                    "n_target": len(target),
                    "selected_features": pipeline.audit_.n_selected_features,
                    "inner_auc_mean": selection.inner_auc_mean,
                    "inner_auc_std": selection.inner_auc_std,
                    **_metrics(target_manifest.outcome_ibd, probability),
                }
            )
            tuning.append(
                {
                    "direction": direction,
                    "split": split_id,
                    "model": family,
                    "inner_auc_mean": selection.inner_auc_mean,
                    "inner_auc_std": selection.inner_auc_std,
                    "parameters_json": json.dumps(selection.selected_parameters, sort_keys=True),
                }
            )
            for sample_id, outcome, score, cluster, origin in zip(
                target.index,
                target_manifest.outcome_ibd,
                probability,
                target_manifest.participant_id,
                target_manifest.origin,
            ):
                predictions.append(
                    {
                        "direction": direction,
                        "split": split_id,
                        "model": family,
                        "sample_id": sample_id,
                        "outcome_ibd": int(outcome),
                        "probability": float(score),
                        "bootstrap_cluster": str(cluster),
                        "origin": str(origin),
                    }
                )
    return metrics, predictions, tuning


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
    prism = load_canonical_microbiome(
        args.legacy_root / "prism_abundance_clean.csv", prism_manifest, shared
    )
    hmp_manifest = _load_manifest(args.contract_dir / "cohort_manifest_hmp2.csv", hmp)
    prism_manifest = _load_manifest(args.contract_dir / "cohort_manifest_prism.csv", prism)
    nielsen_manifest = pd.read_csv(args.nielsen_manifest)
    nielsen = load_canonical_microbiome(args.nielsen_projection, nielsen_manifest, shared)
    nielsen_manifest = _load_manifest(args.nielsen_manifest, nielsen)

    directions = [
        ("HMP2_to_PRISM", hmp, hmp_manifest, prism, prism_manifest, True),
        ("PRISM_to_HMP2", prism, prism_manifest, hmp, hmp_manifest, False),
        ("HMP2_to_NielsenHB_2014", hmp, hmp_manifest, nielsen, nielsen_manifest, True),
        ("PRISM_to_NielsenHB_2014", prism, prism_manifest, nielsen, nielsen_manifest, False),
    ]
    all_metrics: list[dict] = []
    all_predictions: list[dict] = []
    all_tuning: list[dict] = []
    for values in directions:
        metrics, predictions, tuning = _direction_records(*values, n_splits=args.n_splits)
        all_metrics.extend(metrics)
        all_predictions.extend(predictions)
        all_tuning.extend(tuning)

    metric_frame = pd.DataFrame(all_metrics)
    prediction_frame = pd.DataFrame(all_predictions)
    tuning_frame = pd.DataFrame(all_tuning)
    metric_frame.to_csv(out / "standard_model_benchmark_metrics.csv", index=False)
    prediction_frame.to_csv(out / "standard_model_benchmark_predictions.csv", index=False)
    tuning_frame.to_csv(out / "standard_model_benchmark_tuning.csv", index=False)
    aggregate = metric_frame.groupby(["direction", "model"])[
        ["roc_auc", "average_precision", "brier"]
    ].agg(["mean", "std"])
    aggregate.columns = [f"{metric}_{stat}" for metric, stat in aggregate.columns]
    summary = {
        "protocol": "nested_source_only_selection_locked_external_target",
        "experiment_seed": EXPERIMENT_SEED,
        "n_splits": args.n_splits,
        "n_inner_splits": 3,
        "shared_features": len(shared),
        "target_labels_used_for_selection": False,
        "models": sorted(metric_frame.model.unique().tolist()),
        "external_target_summary": aggregate.reset_index().round(6).to_dict(orient="records"),
    }
    (out / "standard_model_benchmark_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
