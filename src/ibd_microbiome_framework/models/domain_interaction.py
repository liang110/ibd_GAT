"""Regularized shared-effect model with centered source-domain deviations."""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression


class DomainInteractionLogisticClassifier:
    """Learn shared effects while allowing source domains to absorb deviations.

    The shared block is the only block available at unseen-target prediction.
    Centered domain interactions are fitted from source data only and therefore
    cannot encode target labels, covariates, or target feature statistics.
    """

    def __init__(self, c: float = 0.1, random_state: int = 42):
        self.c = c
        self.random_state = random_state

    def fit(
        self, x_train: pd.DataFrame, y_train: Sequence[int], domains: Sequence[object]
    ) -> "DomainInteractionLogisticClassifier":
        if len(x_train) != len(y_train) or len(x_train) != len(domains):
            raise ValueError("Features, labels, and domains must align one-to-one")
        self.feature_names_ = x_train.columns.tolist()
        self.domain_levels_ = np.sort(pd.Series(domains).astype(str).unique())
        if len(self.domain_levels_) < 2:
            raise ValueError("At least two source domains are required")
        values = x_train.to_numpy(dtype=float)
        self.domain_proportions_ = {
            domain: float((pd.Series(domains).astype(str) == domain).mean())
            for domain in self.domain_levels_
        }
        design = self._design(values, domains)
        self.model_ = LogisticRegression(
            C=self.c,
            penalty="l2",
            class_weight="balanced",
            max_iter=5000,
            random_state=self.random_state,
        ).fit(design, y_train)
        return self

    def predict_proba(
        self, x: pd.DataFrame, domains: Sequence[object] | None = None
    ) -> np.ndarray:
        if list(x.columns) != self.feature_names_:
            raise ValueError("Prediction features must match the fitted shared feature order")
        values = x.to_numpy(dtype=float)
        if domains is None:
            # Target-domain deviations are unavailable by design; use shared effects only.
            design = np.hstack(
                [values, np.zeros((len(values), values.shape[1] * len(self.domain_levels_)))]
            )
        else:
            if len(domains) != len(x):
                raise ValueError("Prediction domains must align one-to-one with rows")
            unknown = sorted(set(pd.Series(domains).astype(str)).difference(self.domain_levels_))
            if unknown:
                raise ValueError(f"Prediction has unseen source domains: {unknown}")
            design = self._design(values, domains)
        return self.model_.predict_proba(design)

    def _design(self, values: np.ndarray, domains: Sequence[object]) -> np.ndarray:
        domain_values = pd.Series(domains).astype(str).to_numpy()
        deviation_blocks = []
        for domain in self.domain_levels_:
            centered_indicator = (
                (domain_values == domain).astype(float) - self.domain_proportions_[domain]
            )
            deviation_blocks.append(values * centered_indicator[:, None])
        return np.hstack([values, *deviation_blocks])
