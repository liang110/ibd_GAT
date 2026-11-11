#!/usr/bin/env python
"""Validate and combine matched external predictions for paired inference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


KEY_COLUMNS = ["direction", "split", "sample_id"]
MATCH_COLUMNS = ["outcome_ibd", "bootstrap_cluster", "origin"]
REQUIRED_COLUMNS = [*KEY_COLUMNS, "model", *MATCH_COLUMNS, "probability"]


def _read_predictions(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    missing = sorted(set(REQUIRED_COLUMNS).difference(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing prediction columns: {missing}")
    if frame.duplicated([*KEY_COLUMNS, "model"]).any():
        raise ValueError(f"{path} has duplicate model predictions for a target sample")
    return frame.loc[:, REQUIRED_COLUMNS].copy()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-predictions", type=Path, required=True)
    parser.add_argument("--candidate-predictions", type=Path, required=True)
    parser.add_argument("--reference-model", default="elastic_net")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    reference = _read_predictions(args.reference_predictions)
    candidate = _read_predictions(args.candidate_predictions)
    reference = reference.loc[reference.model == args.reference_model].copy()
    if reference.empty:
        raise ValueError(f"Reference model {args.reference_model!r} is absent")
    candidate_models = sorted(candidate.model.unique().tolist())
    if len(candidate_models) != 1:
        raise ValueError("Candidate predictions must contain exactly one model")
    if candidate_models[0] == args.reference_model:
        raise ValueError("Candidate and reference model names must differ")

    reference_indexed = reference.set_index(KEY_COLUMNS).sort_index()
    candidate_indexed = candidate.set_index(KEY_COLUMNS).sort_index()
    if not reference_indexed.index.equals(candidate_indexed.index):
        raise ValueError("Reference and candidate target sample/split keys do not match")
    for column in MATCH_COLUMNS:
        if not reference_indexed[column].equals(candidate_indexed[column]):
            raise ValueError(f"Reference and candidate disagree on {column}")

    paired = pd.concat([reference, candidate], ignore_index=True)
    paired.to_csv(output_dir / "paired_predictions.csv", index=False)
    report = {
        "reference_model": args.reference_model,
        "candidate_model": candidate_models[0],
        "n_predictions_per_model": len(reference),
        "n_direction_split_pairs": int(reference.groupby(["direction", "split"]).ngroups),
    }
    (output_dir / "paired_prediction_validation.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
