import numpy as np
import pandas as pd

from ibd_microbiome_framework.evaluation.standard_benchmark import (
    fit_final_source_models,
    select_standard_models,
)
from ibd_microbiome_framework.models.standard import standard_model_candidates


def _source_matrix() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(12)
    labels = np.array([0, 1] * 24)
    values = rng.uniform(0.01, 1.0, size=(len(labels), 8))
    values[:, 0] += labels * 0.7
    matrix = pd.DataFrame(values, columns=[f"taxon_{index}" for index in range(8)])
    matrix.index = [f"sample_{index}" for index in range(len(matrix))]
    groups = np.array([f"participant_{index // 2}" for index in range(len(matrix))])
    return matrix, labels, groups


def test_standard_candidate_families_have_fixed_grids_without_optional_xgboost():
    candidates = standard_model_candidates(seed=1, include_xgboost=False)
    assert set(candidates) == {
        "prevalence_baseline",
        "logistic_l2",
        "elastic_net",
        "linear_svm",
        "random_forest",
        "shallow_mlp",
    }
    assert len(candidates["prevalence_baseline"]) == 1
    assert all(len(models) == 4 for name, models in candidates.items() if name != "prevalence_baseline")


def test_mlp_early_stopping_is_disabled_by_default():
    candidates = standard_model_candidates(seed=1, include_xgboost=False)
    assert all(model.get_params()["early_stopping"] is False for model in candidates["shallow_mlp"])


def test_mlp_ungrouped_early_stopping_requires_explicit_opt_in():
    candidates = standard_model_candidates(
        seed=1, include_xgboost=False, allow_ungrouped_early_stopping=True
    )
    assert all(model.get_params()["early_stopping"] is True for model in candidates["shallow_mlp"])


def test_nested_selection_and_final_refit_use_source_only_inputs():
    matrix, labels, groups = _source_matrix()
    selected = select_standard_models(
        matrix, labels, groups=groups, seed=5, n_inner_splits=3, include_xgboost=False
    )
    pipeline, fitted = fit_final_source_models(matrix, labels, selected, seed=5)
    transformed = pipeline.transform(matrix.iloc[:4])
    assert set(selected) == set(fitted)
    assert transformed.shape[1] == pipeline.audit_.n_selected_features
    assert all(0.0 <= result.inner_auc_mean <= 1.0 for result in selected.values())
    assert all(hasattr(model, "predict_proba") for model in fitted.values())
