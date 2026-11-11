"""Pre-registered classifier families for source-only benchmarking."""

from __future__ import annotations

from typing import Any

from sklearn.base import ClassifierMixin
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC


def standard_model_candidates(
    seed: int, include_xgboost: bool = True, allow_ungrouped_early_stopping: bool = False
) -> dict[str, list[ClassifierMixin]]:
    """Return an equal four-setting grid for each tunable model family.

    Candidate construction is deliberately independent of cohort labels and target
    data. The prevalence baseline has no hyperparameter selection by design.
    MLP automatic early stopping is disabled by default because scikit-learn's
    internal validation split is sample-based rather than participant-grouped.
    The opt-in flag is retained only for legacy comparisons and must not be used
    for participant-aware claims.
    """
    candidates: dict[str, list[ClassifierMixin]] = {
        "prevalence_baseline": [DummyClassifier(strategy="prior")],
        "logistic_l2": [
            LogisticRegression(C=value, class_weight="balanced", max_iter=5000, random_state=seed)
            for value in (0.01, 0.1, 1.0, 10.0)
        ],
        "elastic_net": [
            LogisticRegression(
                C=value,
                penalty="elasticnet",
                solver="saga",
                l1_ratio=l1_ratio,
                class_weight="balanced",
                max_iter=5000,
                random_state=seed,
            )
            for value, l1_ratio in ((0.1, 0.1), (0.1, 0.9), (1.0, 0.1), (1.0, 0.9))
        ],
        "linear_svm": [
            SVC(C=value, kernel="linear", class_weight="balanced", probability=True, random_state=seed)
            for value in (0.01, 0.1, 1.0, 10.0)
        ],
        "random_forest": [
            RandomForestClassifier(
                n_estimators=300,
                max_depth=depth,
                min_samples_leaf=min_leaf,
                class_weight="balanced",
                n_jobs=1,
                random_state=seed,
            )
            for depth, min_leaf in ((None, 1), (None, 3), (6, 1), (6, 3))
        ],
        "shallow_mlp": [
            MLPClassifier(
                hidden_layer_sizes=(width,),
                alpha=alpha,
                early_stopping=allow_ungrouped_early_stopping,
                validation_fraction=0.15,
                n_iter_no_change=30,
                max_iter=600,
                random_state=seed,
            )
            for width, alpha in ((8, 0.001), (8, 0.01), (16, 0.001), (16, 0.01))
        ],
    }
    if include_xgboost:
        try:
            from xgboost import XGBClassifier
        except ImportError as error:
            raise RuntimeError("XGBoost is required by the registered benchmark model set") from error
        candidates["xgboost"] = [
            XGBClassifier(
                n_estimators=300,
                max_depth=depth,
                learning_rate=learning_rate,
                subsample=0.8,
                colsample_bytree=0.8,
                objective="binary:logistic",
                eval_metric="logloss",
                n_jobs=1,
                random_state=seed,
                tree_method="hist",
            )
            for depth, learning_rate in ((2, 0.03), (2, 0.1), (3, 0.03), (3, 0.1))
        ]
    return candidates


def selected_parameters(model: ClassifierMixin) -> dict[str, Any]:
    """Return JSON-safe, prediction-relevant parameters for the run manifest."""
    parameters = model.get_params(deep=False)
    excluded = {"n_jobs", "random_state", "verbose"}
    return {key: value for key, value in parameters.items() if key not in excluded}
