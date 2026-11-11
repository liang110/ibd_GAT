#!/usr/bin/env python
"""Run positive and negative controls for the semi-synthetic stress test."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

from ibd_microbiome_framework.evaluation.synthetic import make_portability_benchmark
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


SCENARIOS = {
    "signal_only": dict(disease_multiplier=50.0, batch_multiplier=1.0),
    "null": dict(disease_multiplier=1.0, batch_multiplier=1.0),
    "batch_only": dict(disease_multiplier=1.0, batch_multiplier=4.0),
    "signal_reversing_batch": dict(disease_multiplier=50.0, batch_multiplier=2.0),
    "signal_topology_shift": dict(
        disease_multiplier=50.0,
        batch_multiplier=1.0,
        n_topology_features=6,
        target_topology_shift=True,
    ),
}


def evaluate(background: pd.DataFrame, scenario: str, seed: int) -> dict:
    source_x, source_y, target_x, target_y, truth = make_portability_benchmark(
        background,
        random_state=seed,
        min_truth_prevalence=0.30,
        **SCENARIOS[scenario],
    )
    pipeline = CLRFeaturePipeline(min_prevalence=0.10, top_k=20, random_state=seed).fit(
        source_x, source_y
    )
    model = LogisticRegression(
        C=1.0, class_weight="balanced", max_iter=3000, random_state=seed
    ).fit(pipeline.transform(source_x), source_y)
    probability = model.predict_proba(pipeline.transform(target_x))[:, 1]
    auc = float(roc_auc_score(target_y, probability))
    selected = set(pipeline.selected_features_)
    recovered = len(selected.intersection(truth.disease_features))
    batch_selected = len(selected.intersection(truth.batch_features))
    if scenario in {"null", "batch_only"}:
        passed = 0.35 <= auc <= 0.65
        gate = "negative_auc_0.35_to_0.65"
    else:
        passed = auc >= 0.70 and recovered >= 3
        gate = "positive_auc_0.70_and_recovery_3_of_6"
    return {
        "scenario": scenario,
        "seed": seed,
        "target_roc_auc": auc,
        "target_average_precision": float(average_precision_score(target_y, probability)),
        "n_disease_features_recovered": recovered,
        "n_batch_features_selected": batch_selected,
        "acceptance_gate": gate,
        "acceptance_pass": bool(passed),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(range(20260815, 20260825)))
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    background = pd.read_csv(args.background, index_col=0)
    rows = [evaluate(background, scenario, seed) for scenario in SCENARIOS for seed in args.seeds]
    metrics = pd.DataFrame(rows)
    metrics.to_csv(out / "synthetic_control_metrics.csv", index=False)
    aggregate = metrics.groupby("scenario").agg(
        n_seeds=("seed", "size"),
        roc_auc_mean=("target_roc_auc", "mean"),
        roc_auc_std=("target_roc_auc", "std"),
        average_precision_mean=("target_average_precision", "mean"),
        disease_recovery_mean=("n_disease_features_recovered", "mean"),
        batch_selected_mean=("n_batch_features_selected", "mean"),
        acceptance_pass_rate=("acceptance_pass", "mean"),
    ).reset_index()
    aggregate.to_csv(out / "synthetic_control_summary.csv", index=False)
    report = {
        "seeds": args.seeds,
        "benchmark_version": "v2_prevalent_truth_features",
        "min_truth_prevalence": 0.30,
        "scenarios": aggregate.round(6).to_dict(orient="records"),
        "interpretation_rule": "Positive and negative controls must pass before portability stress results are interpretable.",
    }
    (out / "synthetic_control_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
