#!/usr/bin/env python
"""Run pre-registered source-only feature-contract sensitivity analyses.

The 126 shared taxa remain the immutable universe.  For each source split, an
optional feature cap is selected by source-training prevalence only; no target
abundance, metadata, or labels are used to define a variant.
"""

from __future__ import annotations

import argparse
import itertools
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


def _variant_id(feature_cap: int, min_prevalence: float) -> str:
    return f"cap{feature_cap}_prev{min_prevalence:.2f}".replace(".", "p")


def _run_variant_direction(
    *,
    variant_id: str,
    feature_cap: int,
    min_prevalence: float,
    direction: str,
    source: pd.DataFrame,
    source_manifest: pd.DataFrame,
    target: pd.DataFrame,
    target_manifest: pd.DataFrame,
    grouped: bool,
    n_splits: int,
) -> tuple[list[dict], list[dict], list[dict], list[dict], list[dict]]:
    metrics: list[dict] = []
    predictions: list[dict] = []
    tuning: list[dict] = []
    features: list[dict] = []
    weights: list[dict] = []
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
            min_prevalence=min_prevalence,
            top_k=100,
            feature_cap=feature_cap,
        )
        pipeline, fitted = fit_final_source_models(
            source,
            source_manifest.outcome_ibd.to_numpy(),
            selected,
            seed=EXPERIMENT_SEED + split_id,
            min_prevalence=min_prevalence,
            top_k=100,
            feature_cap=feature_cap,
        )
        base = {
            "variant_id": variant_id,
            "feature_cap": feature_cap,
            "min_prevalence": min_prevalence,
            "direction": direction,
            "split": split_id,
            "n_source": len(source),
            "n_selection_train": len(train_idx),
            "n_selection_validation": len(validation_idx),
            "n_target": len(target),
            "n_prevalent_features": pipeline.audit_.n_prevalent_features,
            "n_selected_features": pipeline.audit_.n_selected_features,
        }
        for feature_rank, feature_id in enumerate(pipeline.selected_features_, 1):
            features.append({**base, "feature_rank": feature_rank, "feature_id": feature_id})
        x_target = pipeline.transform(target)
        for family, model in fitted.items():
            probability = model.predict_proba(x_target)[:, 1]
            selection = selected[family]
            metrics.append(
                {
                    **base,
                    "model": family,
                    "evaluation": "locked_external_target",
                    "inner_auc_mean": selection.inner_auc_mean,
                    "inner_auc_std": selection.inner_auc_std,
                    **_metrics(target_manifest.outcome_ibd, probability),
                }
            )
            tuning.append(
                {
                    **base,
                    "model": family,
                    "inner_auc_mean": selection.inner_auc_mean,
                    "inner_auc_std": selection.inner_auc_std,
                    "parameters_json": json.dumps(selection.selected_parameters, sort_keys=True),
                }
            )
            if hasattr(model, "coef_"):
                values = np.asarray(model.coef_).reshape(-1)
            elif hasattr(model, "feature_importances_"):
                values = np.asarray(model.feature_importances_).reshape(-1)
            elif family == "shallow_mlp" and hasattr(model, "coefs_"):
                # Input-layer L2 magnitude is an explicitly labeled descriptive proxy.
                values = np.linalg.norm(np.asarray(model.coefs_[0]), axis=1)
            else:
                values = None
            if values is not None:
                for feature_id, weight in zip(pipeline.selected_features_, values):
                    weights.append({**base, "model": family, "feature_id": feature_id, "feature_weight": float(weight)})
            for sample_id, outcome, score, cluster, origin in zip(
                target.index,
                target_manifest.outcome_ibd,
                probability,
                target_manifest.participant_id,
                target_manifest.origin,
            ):
                predictions.append(
                    {
                        "variant_id": variant_id,
                        "feature_cap": feature_cap,
                        "min_prevalence": min_prevalence,
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
    return metrics, predictions, tuning, features, weights


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-dir", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--nielsen-projection", type=Path, required=True)
    parser.add_argument("--nielsen-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--feature-cap", type=int, required=True)
    parser.add_argument("--min-prevalence", type=float, required=True)
    parser.add_argument("--n-splits", type=int, default=10)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    registered_variants = {(50, 0.10), (100, 0.10), (126, 0.10), (126, 0.05), (126, 0.20)}
    if (args.feature_cap, args.min_prevalence) not in registered_variants:
        raise ValueError("Registered variants are 50@0.10, 100@0.10, 126@0.10, 126@0.05, 126@0.20")

    shared = load_shared_feature_list(args.contract_dir)
    if len(shared) != 126:
        raise ValueError(f"Expected immutable 126-species universe, got {len(shared)}")
    hmp_manifest_path = args.contract_dir / "cohort_manifest_hmp2.csv"
    prism_manifest_path = args.contract_dir / "cohort_manifest_prism.csv"
    hmp_manifest_raw = pd.read_csv(hmp_manifest_path)
    prism_manifest_raw = pd.read_csv(prism_manifest_path)
    nielsen_manifest_raw = pd.read_csv(args.nielsen_manifest)
    hmp = load_canonical_microbiome(args.legacy_root / "manual_abundance.csv", hmp_manifest_raw, shared)
    prism = load_canonical_microbiome(args.legacy_root / "prism_abundance_clean.csv", prism_manifest_raw, shared)
    nielsen = load_canonical_microbiome(args.nielsen_projection, nielsen_manifest_raw, shared)
    hmp_manifest = _load_manifest(hmp_manifest_path, hmp)
    prism_manifest = _load_manifest(prism_manifest_path, prism)
    nielsen_manifest = _load_manifest(args.nielsen_manifest, nielsen)
    directions = [
        ("HMP2_to_PRISM", hmp, hmp_manifest, prism, prism_manifest, True),
        ("PRISM_to_HMP2", prism, prism_manifest, hmp, hmp_manifest, False),
        ("HMP2_to_NielsenHB_2014", hmp, hmp_manifest, nielsen, nielsen_manifest, True),
        ("PRISM_to_NielsenHB_2014", prism, prism_manifest, nielsen, nielsen_manifest, False),
        ("NielsenHB_2014_to_HMP2", nielsen, nielsen_manifest, hmp, hmp_manifest, True),
        ("NielsenHB_2014_to_PRISM", nielsen, nielsen_manifest, prism, prism_manifest, True),
    ]

    all_metrics: list[dict] = []
    all_predictions: list[dict] = []
    all_tuning: list[dict] = []
    all_features: list[dict] = []
    all_weights: list[dict] = []
    feature_cap = args.feature_cap
    min_prevalence = args.min_prevalence
    variant_id = _variant_id(feature_cap, min_prevalence)
    for values in directions:
        metrics, predictions, tuning, features, weights = _run_variant_direction(
            variant_id=variant_id,
            feature_cap=feature_cap,
            min_prevalence=min_prevalence,
            direction=values[0],
            source=values[1],
            source_manifest=values[2],
            target=values[3],
            target_manifest=values[4],
            grouped=values[5],
            n_splits=args.n_splits,
        )
        all_metrics.extend(metrics)
        all_predictions.extend(predictions)
        all_tuning.extend(tuning)
        all_features.extend(features)
        all_weights.extend(weights)
    pd.DataFrame(all_metrics).to_csv(out / "contract_sensitivity_metrics.csv", index=False)
    pd.DataFrame(all_predictions).to_csv(out / "contract_sensitivity_predictions.csv", index=False)
    pd.DataFrame(all_tuning).to_csv(out / "contract_sensitivity_tuning.csv", index=False)
    pd.DataFrame(all_features).to_csv(out / "contract_sensitivity_features.csv", index=False)
    pd.DataFrame(all_weights).to_csv(out / "contract_sensitivity_feature_weights.csv", index=False)
    summary = {
        "protocol": "registered_feature_contract_sensitivity_source_only",
        "immutable_shared_species_universe": len(shared),
        "feature_cap": feature_cap,
        "min_prevalence": min_prevalence,
        "n_splits": args.n_splits,
        "target_labels_used_for_selection": False,
        "directions": [item[0] for item in directions],
    }
    (out / "contract_sensitivity_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
