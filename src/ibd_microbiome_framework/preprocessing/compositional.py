"""Training-only compositional preprocessing and supervised selection."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd
from sklearn.feature_selection import mutual_info_classif
from sklearn.preprocessing import StandardScaler

from ibd_microbiome_framework.harmonization.features import is_forbidden_feature


@dataclass
class TransformAudit:
    n_input_features: int
    n_prevalent_features: int
    n_selected_features: int
    pseudo_count: float
    training_samples: int


class CLRFeaturePipeline:
    """Fit prevalence, CLR pseudo-count, scale, and MI only on training rows."""

    def __init__(
        self,
        min_prevalence: float = 0.10,
        top_k: int = 100,
        random_state: int = 42,
        feature_cap: int | None = None,
    ):
        self.min_prevalence = min_prevalence
        self.top_k = top_k
        self.random_state = random_state
        self.feature_cap = feature_cap

    def fit(self, x_train: pd.DataFrame, y_train: Sequence[int]) -> "CLRFeaturePipeline":
        if not x_train.index.is_unique or not x_train.columns.is_unique:
            raise ValueError("Training matrix requires unique sample and feature IDs")
        forbidden = [str(column) for column in x_train.columns if is_forbidden_feature(column)]
        if forbidden:
            raise ValueError(f"Clinical, label, or identifier features are forbidden: {forbidden}")
        if len(y_train) != len(x_train):
            raise ValueError("Training labels must align one-to-one with rows")
        prevalence = (x_train > 0).mean(axis=0)
        if self.feature_cap is not None:
            if self.feature_cap < 2:
                raise ValueError("feature_cap must retain at least two features")
            # Break prevalence ties by canonical feature ID.
            ranked = sorted(prevalence.index, key=lambda name: (-prevalence[name], str(name)))
            retained = ranked[: min(self.feature_cap, len(ranked))]
            prevalence = prevalence.loc[retained]
        self.prevalent_features_ = prevalence[prevalence >= self.min_prevalence].index.tolist()
        if len(self.prevalent_features_) < 2:
            raise ValueError("Prevalence filtering left fewer than two features")
        train = x_train.loc[:, self.prevalent_features_]
        positive = train.to_numpy()[train.to_numpy() > 0]
        if positive.size == 0:
            raise ValueError("No positive values after prevalence filtering")
        self.pseudo_count_ = float(np.min(positive) * 0.5)
        clr = self._clr(train)
        self.scaler_ = StandardScaler().fit(clr)
        scaled = self.scaler_.transform(clr)
        k = min(self.top_k, scaled.shape[1])
        scores = mutual_info_classif(scaled, np.asarray(y_train), random_state=self.random_state)
        order = np.argsort(scores)[::-1][:k]
        self.selected_features_ = [self.prevalent_features_[idx] for idx in order]
        self.audit_ = TransformAudit(
            n_input_features=x_train.shape[1],
            n_prevalent_features=len(self.prevalent_features_),
            n_selected_features=len(self.selected_features_),
            pseudo_count=self.pseudo_count_,
            training_samples=x_train.shape[0],
        )
        return self

    def transform(self, x: pd.DataFrame) -> pd.DataFrame:
        if not hasattr(self, "selected_features_"):
            raise RuntimeError("Pipeline must be fitted before transform")
        forbidden = [str(column) for column in x.columns if is_forbidden_feature(column)]
        if forbidden:
            raise ValueError(f"Clinical, label, or identifier features are forbidden: {forbidden}")
        missing = sorted(set(self.prevalent_features_).difference(x.columns))
        if missing:
            raise ValueError(f"Transform matrix lacks {len(missing)} fitted features")
        clr = self._clr(x.loc[:, self.prevalent_features_])
        scaled = self.scaler_.transform(clr)
        frame = pd.DataFrame(scaled, index=x.index, columns=self.prevalent_features_)
        return frame.loc[:, self.selected_features_]

    def _clr(self, x: pd.DataFrame) -> np.ndarray:
        values = x.to_numpy(dtype=float)
        logged = np.log(values + self.pseudo_count_)
        return logged - logged.mean(axis=1, keepdims=True)
