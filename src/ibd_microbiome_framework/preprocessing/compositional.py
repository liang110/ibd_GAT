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
            # The cap is calculated from the current training fold only.
            # Canonical feature IDs give deterministic tie handling.
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


class DomainStableCLRFeaturePipeline(CLRFeaturePipeline):
    """Select features with directionally stable training-domain effects."""

    def fit(
        self,
        x_train: pd.DataFrame,
        y_train: Sequence[int],
        domains: Sequence[object],
    ) -> "DomainStableCLRFeaturePipeline":
        super().fit(x_train, y_train)
        labels = np.asarray(y_train)
        domain_values = np.asarray(domains).astype(str)
        if len(domain_values) != len(x_train):
            raise ValueError("Training domains must align one-to-one with rows")
        scaled = self.scaler_.transform(self._clr(x_train.loc[:, self.prevalent_features_]))
        scaled_frame = pd.DataFrame(scaled, index=x_train.index, columns=self.prevalent_features_)
        effects = {}
        for domain in sorted(np.unique(domain_values)):
            mask = domain_values == domain
            if np.unique(labels[mask]).size != 2:
                continue
            effects[domain] = (
                scaled_frame.loc[mask & (labels == 1)].mean(axis=0)
                - scaled_frame.loc[mask & (labels == 0)].mean(axis=0)
            )
        if len(effects) < 2:
            raise ValueError("At least two training domains containing both outcomes are required")
        effect_frame = pd.DataFrame(effects).T
        same_direction = (effect_frame.gt(0).all(axis=0) | effect_frame.lt(0).all(axis=0))
        stable = effect_frame.loc[:, same_direction]
        if stable.shape[1] < 1:
            raise ValueError("No directionally stable features remain")
        ranking = stable.abs().min(axis=0).sort_values(ascending=False)
        self.selected_features_ = ranking.head(min(self.top_k, len(ranking))).index.tolist()
        self.domain_effects_ = effect_frame.loc[:, self.selected_features_]
        self.stable_feature_count_ = int(stable.shape[1])
        self.audit_.n_selected_features = len(self.selected_features_)
        return self
