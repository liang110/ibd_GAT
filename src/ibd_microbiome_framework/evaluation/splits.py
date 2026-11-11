"""Deterministic split plans that preserve grouping and held-out cohorts."""

from __future__ import annotations

from typing import Iterator, Tuple

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit, StratifiedShuffleSplit


def repeated_group_train_validation(
    manifest: pd.DataFrame, n_splits: int, validation_size: float, seed: int
) -> Iterator[Tuple[int, np.ndarray, np.ndarray]]:
    """Yield source train/validation indices with both outcomes in each partition."""
    group_column = "participant_id"
    if manifest[group_column].isna().any():
        raise ValueError("participant_id is required for grouped splitting")
    y = manifest["outcome_ibd"].to_numpy()
    groups = manifest[group_column].astype(str).to_numpy()
    for split_id in range(n_splits):
        chosen = None
        for attempt in range(100):
            splitter = GroupShuffleSplit(
                n_splits=1,
                test_size=validation_size,
                random_state=seed + split_id * 1009 + attempt,
            )
            train_idx, validation_idx = next(splitter.split(np.zeros(len(y)), y, groups))
            if len(np.unique(y[train_idx])) == 2 and len(np.unique(y[validation_idx])) == 2:
                chosen = (train_idx, validation_idx)
                break
        if chosen is None:
            raise RuntimeError(f"Unable to create valid grouped split {split_id}")
        yield split_id + 1, chosen[0], chosen[1]


def prism_origin_summary(manifest: pd.DataFrame) -> dict:
    if "origin" not in manifest.columns:
        raise ValueError("PRISM origin missing")
    return manifest.groupby(["origin", "outcome_ibd"]).size().to_dict()


def repeated_stratified_train_validation(
    manifest: pd.DataFrame,
    n_splits: int,
    validation_size: float,
    seed: int,
    strata_columns: tuple[str, ...] = ("outcome_ibd",),
) -> Iterator[Tuple[int, np.ndarray, np.ndarray]]:
    """Yield deterministic source splits preserving declared categorical strata."""
    missing = sorted(set(strata_columns).difference(manifest.columns))
    if missing:
        raise ValueError(f"Manifest lacks stratification columns: {missing}")
    strata = manifest.loc[:, list(strata_columns)].astype(str).agg("__".join, axis=1)
    if strata.value_counts().min() < 2:
        raise ValueError("Every stratum requires at least two samples")
    splitter = StratifiedShuffleSplit(
        n_splits=n_splits, test_size=validation_size, random_state=seed
    )
    for split_id, (train_idx, validation_idx) in enumerate(
        splitter.split(np.zeros(len(manifest)), strata), 1
    ):
        yield split_id, train_idx, validation_idx
