#!/usr/bin/env python
"""Summarize source-only hyperparameter selection stability across splits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


KEYS = {
    "prevalence_baseline": (),
    "logistic_l2": ("C",),
    "elastic_net": ("C", "l1_ratio"),
    "linear_svm": ("C",),
    "random_forest": ("max_depth", "min_samples_leaf"),
    "xgboost": ("max_depth", "learning_rate"),
    "shallow_mlp": ("hidden_layer_sizes", "alpha"),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tuning", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tuning = pd.read_csv(args.tuning)
    keys = []
    for row in tuning.itertuples(index=False):
        parameters = json.loads(row.parameters_json)
        chosen = {key: parameters.get(key) for key in KEYS[row.model]}
        keys.append(json.dumps(chosen, sort_keys=True))
    tuning = tuning.assign(parameter_choice=keys)
    frequency = (
        tuning.groupby(["direction", "model", "parameter_choice"], as_index=False)
        .agg(n_splits=("split", "nunique"))
    )
    totals = tuning.groupby(["direction", "model"], as_index=False).agg(n_total_splits=("split", "nunique"))
    frequency = frequency.merge(totals, on=["direction", "model"], how="left")
    frequency["selection_frequency"] = frequency.n_splits / frequency.n_total_splits
    summary = (
        frequency.sort_values(["direction", "model", "selection_frequency"], ascending=[True, True, False])
        .groupby(["direction", "model"], as_index=False)
        .first()
        .rename(columns={"parameter_choice": "most_frequent_parameter_choice", "selection_frequency": "most_frequent_choice_frequency"})
    )
    frequency.to_csv(args.output_dir / "source_tuning_selection_frequency.csv", index=False)
    summary.to_csv(args.output_dir / "source_tuning_stability_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
