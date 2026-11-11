#!/usr/bin/env python
"""Run the native SIAMCAT source-only reference on the locked IBD contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


MODEL = "siamcat_lasso"
SEED_BASE = 20260816


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def metrics(outcomes: pd.Series, probabilities: pd.Series) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(outcomes, probabilities)),
        "average_precision": float(average_precision_score(outcomes, probabilities)),
        "brier": float(brier_score_loss(outcomes, probabilities)),
    }


def git_revision(project_root: Path) -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=project_root, text=True
    ).strip()


def run_direction(
    *,
    direction: str,
    source_name: str,
    target_name: str,
    inputs: Path,
    output: Path,
    rscript: str,
    reference_script: Path,
    n_splits: int,
    config_sha256: str,
) -> tuple[list[dict[str, object]], list[pd.DataFrame]]:
    source_matrix = inputs / f"{source_name}_shared_species.csv"
    source_metadata = inputs / f"{source_name}_source_metadata.csv"
    target_matrix = inputs / f"{target_name}_shared_species.csv"
    target_metadata = inputs / f"{target_name}_target_metadata.csv"
    target_outcomes = pd.read_csv(inputs / f"{target_name}_locked_outcomes.csv")
    expected_ids = target_outcomes.sample_id.tolist()
    metric_rows: list[dict[str, object]] = []
    prediction_frames: list[pd.DataFrame] = []
    for split_id in range(1, n_splits + 1):
        seed = SEED_BASE + split_id
        split_output = output / "models" / direction / f"split_{split_id:02d}"
        split_output.mkdir(parents=True, exist_ok=True)
        command = [
            rscript,
            "--vanilla",
            str(reference_script),
            str(source_matrix),
            str(source_metadata),
            str(target_matrix),
            str(target_metadata),
            str(split_output),
            str(seed),
        ]
        subprocess.run(command, check=True)
        prediction = pd.read_csv(split_output / "siamcat_external_predictions.csv")
        if prediction.sample_id.tolist() != expected_ids:
            raise RuntimeError(f"{direction} split {split_id}: target sample order changed")
        if prediction.probability.isna().any() or not prediction.probability.between(0, 1).all():
            raise RuntimeError(f"{direction} split {split_id}: invalid SIAMCAT probabilities")
        merged = prediction.merge(
            target_outcomes, on="sample_id", how="left", validate="one_to_one"
        )
        if merged.outcome_ibd.isna().any():
            raise RuntimeError(f"{direction} split {split_id}: outcome join failed")
        merged.insert(0, "direction", direction)
        merged.insert(1, "split", split_id)
        merged.insert(2, "model", MODEL)
        merged["outcome_ibd"] = merged.outcome_ibd.astype(int)
        metric_rows.append(
            {
                "direction": direction,
                "split": split_id,
                "model": MODEL,
                "evaluation": "locked_external_target",
                "n_source": int(pd.read_csv(source_metadata).shape[0]),
                "n_target": int(len(merged)),
                "n_target_clusters": int(merged.participant_id.nunique()),
                "source_ensemble_folds": 5,
                "source_seed": seed,
                "config_sha256": config_sha256,
                **metrics(merged.outcome_ibd, merged.probability),
            }
        )
        merged = merged.rename(columns={"participant_id": "bootstrap_cluster"})
        prediction_frames.append(
            merged.loc[
                :, [
                    "direction", "split", "model", "sample_id", "outcome_ibd",
                    "probability", "bootstrap_cluster", "origin",
                ]
            ]
        )
    return metric_rows, prediction_frames


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--reference-script", type=Path, required=True)
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--n-splits", type=int, default=10)
    args = parser.parse_args()
    inputs = args.inputs_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    contract = json.loads((inputs / "siamcat_input_manifest.json").read_text(encoding="utf-8"))
    if contract["shared_features"] != 126:
        raise ValueError("SIAMCAT reference requires the locked 126-species contract")
    config_sha256 = sha256(args.config)
    directions = [
        ("HMP2_to_PRISM", "hmp2", "prism"),
        ("PRISM_to_HMP2", "prism", "hmp2"),
        ("HMP2_to_NielsenHB_2014", "hmp2", "nielsenhb_2014"),
        ("PRISM_to_NielsenHB_2014", "prism", "nielsenhb_2014"),
    ]
    all_metrics: list[dict[str, object]] = []
    all_predictions: list[pd.DataFrame] = []
    for direction, source_name, target_name in directions:
        rows, frames = run_direction(
            direction=direction,
            source_name=source_name,
            target_name=target_name,
            inputs=inputs,
            output=output,
            rscript=args.rscript,
            reference_script=args.reference_script.resolve(),
            n_splits=args.n_splits,
            config_sha256=config_sha256,
        )
        all_metrics.extend(rows)
        all_predictions.extend(frames)
    metric_frame = pd.DataFrame(all_metrics)
    prediction_frame = pd.concat(all_predictions, ignore_index=True)
    metric_frame.to_csv(output / "siamcat_reference_metrics.csv", index=False)
    prediction_frame.to_csv(output / "siamcat_reference_predictions.csv", index=False)
    aggregate = metric_frame.groupby(["direction", "model"])[
        ["roc_auc", "average_precision", "brier"]
    ].agg(["mean", "std"])
    aggregate.columns = [f"{metric}_{stat}" for metric, stat in aggregate.columns]
    report = {
        "protocol": "native_siamcat_source_only_5_fold_ensemble_repeated_source_seeds",
        "reference_implementation": "SIAMCAT_lasso",
        "n_source_seed_perturbations": args.n_splits,
        "source_ensemble_folds_per_perturbation": 5,
        "shared_features": contract["shared_features"],
        "target_labels_available_to_r_training_process": False,
        "config_sha256": config_sha256,
        "input_manifest_sha256": sha256(inputs / "siamcat_input_manifest.json"),
        "reference_script_sha256": sha256(args.reference_script),
        "git_revision": git_revision(Path.cwd()),
        "external_target_summary": aggregate.reset_index().round(6).to_dict(orient="records"),
    }
    (output / "siamcat_reference_summary.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
