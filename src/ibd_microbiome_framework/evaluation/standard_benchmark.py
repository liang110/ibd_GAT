"""Nested source-only selection for the pre-registered standard model families."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from sklearn.base import ClassifierMixin, clone
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold

from ibd_microbiome_framework.models.standard import selected_parameters, standard_model_candidates
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


@dataclass(frozen=True)
class SelectedModel:
    """A candidate selected without target-cohort information."""

    family: str
    estimator: ClassifierMixin
    inner_auc_mean: float
    inner_auc_std: float
    selected_parameters: dict


def _inner_splits(
    y: Sequence[int], groups: Sequence[object] | None, n_splits: int, seed: int
) -> Iterable[tuple[np.ndarray, np.ndarray]]:
    labels = np.asarray(y)
    if np.unique(labels).size != 2:
        raise ValueError("Source training set must contain both outcomes")
    if groups is None:
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        yield from splitter.split(np.zeros(len(labels)), labels)
        return
    group_values = np.asarray(groups).astype(str)
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    yielded = False
    for train_idx, validation_idx in splitter.split(np.zeros(len(labels)), labels, group_values):
        if np.unique(labels[train_idx]).size == 2 and np.unique(labels[validation_idx]).size == 2:
            yielded = True
            yield train_idx, validation_idx
    if not yielded:
        raise RuntimeError("No valid grouped inner folds with both IBD outcomes")


def select_standard_models(
    x_source: pd.DataFrame,
    y_source: Sequence[int],
    groups: Sequence[object] | None,
    seed: int,
    n_inner_splits: int = 3,
    include_xgboost: bool = True,
    min_prevalence: float = 0.10,
    top_k: int = 100,
    feature_cap: int | None = None,
    allow_ungrouped_early_stopping: bool = False,
) -> dict[str, SelectedModel]:
    """Select exactly one candidate per family using source-only nested folds.

    CLR transformation and MI feature selection are refit inside every inner
    training fold. This ensures the held-out inner fold cannot influence a
    candidate's preprocessing state, even indirectly.
    """
    labels = np.asarray(y_source)
    candidates = standard_model_candidates(
        seed,
        include_xgboost=include_xgboost,
        allow_ungrouped_early_stopping=allow_ungrouped_early_stopping,
    )
    selections: dict[str, SelectedModel] = {}
    inner_folds = list(_inner_splits(labels, groups, n_inner_splits, seed))
    for family, family_candidates in candidates.items():
        scored: list[tuple[float, float, ClassifierMixin]] = []
        for candidate in family_candidates:
            fold_scores: list[float] = []
            for fold_id, (train_idx, validation_idx) in enumerate(inner_folds, 1):
                pipeline = CLRFeaturePipeline(
                    min_prevalence=min_prevalence,
                    top_k=top_k,
                    random_state=seed + fold_id,
                    feature_cap=feature_cap,
                ).fit(x_source.iloc[train_idx], labels[train_idx])
                model = clone(candidate)
                model.fit(pipeline.transform(x_source.iloc[train_idx]), labels[train_idx])
                probability = model.predict_proba(pipeline.transform(x_source.iloc[validation_idx]))[:, 1]
                fold_scores.append(float(roc_auc_score(labels[validation_idx], probability)))
            scored.append((float(np.mean(fold_scores)), float(np.std(fold_scores, ddof=0)), candidate))
        # Deterministic first-candidate tie handling preserves the declared grid order.
        mean_auc, std_auc, winner = max(scored, key=lambda item: item[0])
        selections[family] = SelectedModel(
            family=family,
            estimator=clone(winner),
            inner_auc_mean=mean_auc,
            inner_auc_std=std_auc,
            selected_parameters=selected_parameters(winner),
        )
    return selections


def fit_final_source_models(
    x_source: pd.DataFrame,
    y_source: Sequence[int],
    selections: dict[str, SelectedModel],
    seed: int,
    min_prevalence: float = 0.10,
    top_k: int = 100,
    feature_cap: int | None = None,
) -> tuple[CLRFeaturePipeline, dict[str, ClassifierMixin]]:
    """Refit selected pipelines on all and only source samples before target scoring."""
    pipeline = CLRFeaturePipeline(
        min_prevalence=min_prevalence,
        top_k=top_k,
        random_state=seed,
        feature_cap=feature_cap,
    ).fit(x_source, y_source)
    transformed_source = pipeline.transform(x_source)
    fitted = {}
    for family, selection in selections.items():
        model = clone(selection.estimator)
        model.fit(transformed_source, y_source)
        fitted[family] = model
    return pipeline, fitted
