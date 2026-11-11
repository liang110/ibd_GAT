#!/usr/bin/env python
"""Refit the frozen model families on the Nielsen ESP subset.

This is a fixed-contract, source-only sensitivity audit.  It does not rebuild
the taxonomic universe: all three cohorts are projected onto the registered
126-species contract before source fitting.  HMP2 and PRISM are locked targets.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list
from ibd_microbiome_framework.harmonization.features import canonical_taxon
from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation
from ibd_microbiome_framework.evaluation.standard_benchmark import fit_final_source_models, select_standard_models


TRAINABLE_MODELS = (
    "logistic_l2",
    "elastic_net",
    "linear_svm",
    "random_forest",
    "xgboost",
    "shallow_mlp",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_contract_projection(path: Path, required: list[str], output: Path) -> None:
    """Strip non-abundance metadata rows and retain the registered contract."""
    raw = pd.read_csv(path, index_col=0, low_memory=False)
    numeric = raw.apply(pd.to_numeric, errors="coerce")
    valid = raw.notna().sum(axis=1).gt(0) & numeric.notna().all(axis=1)
    numeric = numeric.loc[valid]
    numeric.index = [canonical_taxon(value) for value in numeric.index]
    numeric = numeric.loc[numeric.index.notna()]
    numeric = numeric.groupby(level=0, sort=True).sum()
    missing = sorted(set(required).difference(numeric.index))
    if missing:
        raise ValueError(f"{path.name} lacks {len(missing)} contract taxa after metadata filtering")
    output.parent.mkdir(parents=True, exist_ok=True)
    numeric.loc[required].to_csv(output, index=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-dir", type=Path, required=True)
    parser.add_argument("--hmp2-matrix", type=Path, required=True)
    parser.add_argument("--hmp2-manifest", type=Path, required=True)
    parser.add_argument("--prism-matrix", type=Path, required=True)
    parser.add_argument("--prism-manifest", type=Path, required=True)
    parser.add_argument("--nielsen-matrix", type=Path, required=True)
    parser.add_argument("--nielsen-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-splits", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20261030)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    shared = load_shared_feature_list(args.contract_dir)
    hmp_manifest = pd.read_csv(args.hmp2_manifest)
    prism_manifest = pd.read_csv(args.prism_manifest)
    nielsen_manifest = pd.read_csv(args.nielsen_manifest)
    hmp = load_canonical_microbiome(args.hmp2_matrix, hmp_manifest, shared)
    prism_projection = out / "prism_contract_projection.csv"
    _write_contract_projection(args.prism_matrix, shared, prism_projection)
    prism = load_canonical_microbiome(prism_projection, prism_manifest, shared)
    nielsen = load_canonical_microbiome(args.nielsen_matrix, nielsen_manifest, shared)
    def _ordered_manifest(frame: pd.DataFrame, sample_ids: pd.Index) -> pd.DataFrame:
        ordered = frame.set_index("sample_id").loc[sample_ids].copy()
        ordered.insert(0, "sample_id", ordered.index.astype(str))
        return ordered.reset_index(drop=True)

    hmp_manifest = _ordered_manifest(hmp_manifest, hmp.index)
    prism_manifest = _ordered_manifest(prism_manifest, prism.index)
    nielsen_manifest = _ordered_manifest(nielsen_manifest, nielsen.index)

    esp_ids = set(nielsen_manifest.loc[nielsen_manifest["country"].eq("ESP"), "sample_id"])
    esp_manifest = nielsen_manifest.loc[nielsen_manifest["sample_id"].isin(esp_ids)].reset_index(drop=True)
    esp = nielsen.loc[esp_manifest["sample_id"].tolist()]
    if esp_manifest["outcome_ibd"].nunique() != 2:
        raise ValueError("ESP source subset must contain both outcomes")

    metrics: list[dict[str, object]] = []
    predictions: list[dict[str, object]] = []
    tuning: list[dict[str, object]] = []
    targets = (
        ("ESP_NielsenHB_2014_to_HMP2", hmp, hmp_manifest),
        ("ESP_NielsenHB_2014_to_PRISM", prism, prism_manifest),
    )
    for split_id, train_idx, validation_idx in repeated_group_train_validation(
        esp_manifest, args.n_splits, 0.20, args.seed
    ):
        selection_train = esp.iloc[train_idx]
        selection_y = esp_manifest.iloc[train_idx]["outcome_ibd"].to_numpy()
        groups = esp_manifest.iloc[train_idx]["participant_id"].to_numpy()
        selected = select_standard_models(
            selection_train,
            selection_y,
            groups=groups,
            seed=args.seed + split_id,
            n_inner_splits=3,
        )
        pipeline, fitted = fit_final_source_models(
            esp,
            esp_manifest["outcome_ibd"].to_numpy(),
            selected,
            seed=args.seed + split_id,
        )
        for direction, target, target_manifest in targets:
            transformed = pipeline.transform(target)
            for family in TRAINABLE_MODELS:
                probability = fitted[family].predict_proba(transformed)[:, 1]
                y = target_manifest["outcome_ibd"].to_numpy()
                # Metrics are deliberately kept split-specific; the target is
                # never used for selection, refitting or calibration.
                from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

                metrics.append(
                    {
                        "direction": direction,
                        "split": split_id,
                        "model": family,
                        "n_source_esp_samples": int(len(esp)),
                        "n_source_esp_participants": int(esp_manifest["participant_id"].nunique()),
                        "n_selection_train": int(len(train_idx)),
                        "n_selection_validation": int(len(validation_idx)),
                        "n_target": int(len(target)),
                        "n_target_participants": int(target_manifest["participant_id"].nunique()),
                        "selected_features": int(pipeline.audit_.n_selected_features),
                        "inner_auc_mean": float(selected[family].inner_auc_mean),
                        "inner_auc_std": float(selected[family].inner_auc_std),
                        "roc_auc": float(roc_auc_score(y, probability)),
                        "average_precision": float(average_precision_score(y, probability)),
                        "brier": float(brier_score_loss(y, probability)),
                    }
                )
                tuning.append(
                    {
                        "direction": direction,
                        "split": split_id,
                        "model": family,
                        "inner_auc_mean": float(selected[family].inner_auc_mean),
                        "inner_auc_std": float(selected[family].inner_auc_std),
                        "parameters_json": json.dumps(selected[family].selected_parameters, sort_keys=True),
                    }
                )
                for sample_id, outcome, score, cluster, origin in zip(
                    target.index,
                    target_manifest["outcome_ibd"],
                    probability,
                    target_manifest["participant_id"],
                    target_manifest["origin"],
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

    metric_frame = pd.DataFrame(metrics)
    prediction_frame = pd.DataFrame(predictions)
    tuning_frame = pd.DataFrame(tuning)
    metric_frame.to_csv(out / "esp_source_refit_metrics.csv", index=False)
    prediction_frame.to_csv(out / "esp_source_refit_predictions.csv", index=False)
    tuning_frame.to_csv(out / "esp_source_refit_tuning.csv", index=False)
    aggregate = metric_frame.groupby(["direction", "model"])[
        ["roc_auc", "average_precision", "brier"]
    ].agg(["mean", "std"])
    aggregate.columns = [f"{metric}_{stat}" for metric, stat in aggregate.columns]
    aggregate.reset_index().to_csv(out / "esp_source_refit_aggregate.csv", index=False)
    manifest = {
        "protocol": "fixed_126_contract_ESP_source_only_refit_locked_HMP2_PRISM_targets",
        "seed": args.seed,
        "n_splits": args.n_splits,
        "n_inner_splits": 3,
        "n_fixed_features": len(shared),
        "source_country": "ESP",
        "source_samples": int(len(esp)),
        "source_participants": int(esp_manifest["participant_id"].nunique()),
        "source_outcome_counts": esp_manifest["outcome_ibd"].value_counts().sort_index().to_dict(),
        "target_labels_used_for_selection_or_refit": False,
        "models_refit": list(TRAINABLE_MODELS),
        "strict_source_specific_feature_universe_estimated": False,
        "input_sha256": {
            "hmp2_matrix": sha256_file(args.hmp2_matrix),
            "prism_matrix": sha256_file(args.prism_matrix),
            "prism_contract_projection": sha256_file(prism_projection),
            "nielsen_matrix": sha256_file(args.nielsen_matrix),
            "hmp2_manifest": sha256_file(args.hmp2_manifest),
            "prism_manifest": sha256_file(args.prism_manifest),
            "nielsen_manifest": sha256_file(args.nielsen_manifest),
        },
    }
    (out / "esp_source_refit_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
