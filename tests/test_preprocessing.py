import pandas as pd
import pytest

from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


def test_pipeline_fits_and_transforms_with_stable_columns():
    x_train = pd.DataFrame({"a": [0.7, 0.1, 0.6, 0.0], "b": [0.3, 0.9, 0.4, 1.0]})
    x_test = pd.DataFrame({"a": [0.2, 0.0], "b": [0.8, 1.0]})
    fitted = CLRFeaturePipeline(min_prevalence=0.25, top_k=2, random_state=1).fit(x_train, [0, 1, 0, 1])
    transformed = fitted.transform(x_test)
    assert list(transformed.columns) == fitted.selected_features_
    assert transformed.shape == (2, 2)


@pytest.mark.parametrize("forbidden", ["Age", "Fecal.Calprotectin", "diagnosis", "sample_id"])
def test_pipeline_rejects_planted_clinical_or_identifier_features(forbidden):
    matrix = pd.DataFrame({"taxon_a": [0.7, 0.1, 0.6, 0.0], "taxon_b": [0.3, 0.9, 0.4, 1.0], forbidden: [30, 31, 32, 33]})
    with pytest.raises(ValueError, match="forbidden"):
        CLRFeaturePipeline(min_prevalence=0.25, top_k=2, random_state=1).fit(matrix, [0, 1, 0, 1])


def test_transform_fails_when_a_fitted_feature_is_absent():
    train = pd.DataFrame({"a": [0.7, 0.1, 0.6, 0.0], "b": [0.3, 0.9, 0.4, 1.0]})
    fitted = CLRFeaturePipeline(min_prevalence=0.25, top_k=2, random_state=1).fit(train, [0, 1, 0, 1])
    with pytest.raises(ValueError, match="lacks"):
        fitted.transform(pd.DataFrame({"a": [0.2, 0.0]}))
