#!/usr/bin/env python
"""Audit repeated-subject and planted-clinical leakage controls."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupShuffleSplit, train_test_split
from sklearn.neighbors import KNeighborsClassifier

from ibd_microbiome_framework.evaluation.synthetic import make_repeated_subject_benchmark
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


def _fit_score(x: pd.DataFrame, y: pd.Series, train_idx: np.ndarray, test_idx: np.ndarray, seed: int) -> float:
    pipeline = CLRFeaturePipeline(min_prevalence=0.10, top_k=50, random_state=seed).fit(
        x.iloc[train_idx], y.iloc[train_idx]
    )
    model = KNeighborsClassifier(n_neighbors=1).fit(
        pipeline.transform(x.iloc[train_idx]), y.iloc[train_idx]
    )
    probability = model.predict_proba(pipeline.transform(x.iloc[test_idx]))[:, 1]
    return float(roc_auc_score(y.iloc[test_idx], probability))


def evaluate_seed(background: pd.DataFrame, seed: int) -> dict:
    x, y, subjects = make_repeated_subject_benchmark(background, random_state=seed)
    indices = np.arange(len(x))
    random_train, random_test = train_test_split(
        indices, test_size=0.25, stratify=y, random_state=seed
    )
    grouped = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=seed)
    group_train, group_test = next(grouped.split(x, y, groups=subjects))
    random_overlap = set(subjects.iloc[random_train]).intersection(subjects.iloc[random_test])
    grouped_overlap = set(subjects.iloc[group_train]).intersection(subjects.iloc[group_test])
    planted = x.copy()
    rng = np.random.default_rng(seed)
    planted["Fecal.Calprotectin"] = y.to_numpy() * 1000 + rng.normal(0, 5, len(y))
    leakage_rejected = False
    try:
        CLRFeaturePipeline(min_prevalence=0.10, top_k=50, random_state=seed).fit(planted, y)
    except ValueError as exc:
        leakage_rejected = "forbidden" in str(exc)
    return {
        "seed": seed,
        "random_split_subject_overlap": len(random_overlap),
        "group_split_subject_overlap": len(grouped_overlap),
        "random_split_roc_auc": _fit_score(x, y, random_train, random_test, seed),
        "group_split_roc_auc": _fit_score(x, y, group_train, group_test, seed),
        "clinical_leakage_rejected": leakage_rejected,
        "acceptance_pass": bool(len(random_overlap) > 0 and len(grouped_overlap) == 0 and leakage_rejected),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(range(20260815, 20260825)))
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    background = pd.read_csv(args.background, index_col=0)
    metrics = pd.DataFrame([evaluate_seed(background, seed) for seed in args.seeds])
    metrics.to_csv(out / "repeated_leakage_metrics.csv", index=False)
    summary = {
        "n_seeds": len(metrics),
        "acceptance_pass_rate": float(metrics.acceptance_pass.mean()),
        "clinical_leakage_rejection_rate": float(metrics.clinical_leakage_rejected.mean()),
        "random_split_subject_overlap_mean": float(metrics.random_split_subject_overlap.mean()),
        "group_split_subject_overlap_max": int(metrics.group_split_subject_overlap.max()),
        "random_split_roc_auc_mean": float(metrics.random_split_roc_auc.mean()),
        "group_split_roc_auc_mean": float(metrics.group_split_roc_auc.mean()),
    }
    (out / "repeated_leakage_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
