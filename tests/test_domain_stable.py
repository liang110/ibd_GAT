import numpy as np
import pandas as pd

from ibd_microbiome_framework.preprocessing.compositional import DomainStableCLRFeaturePipeline


def test_domain_stable_selector_rejects_reversing_feature():
    labels = np.array([0, 0, 1, 1] * 2)
    domains = np.array(["a"] * 4 + ["b"] * 4)
    stable = np.where(labels == 1, 0.8, 0.2)
    reversing = np.array([0.2, 0.2, 0.8, 0.8, 0.8, 0.8, 0.2, 0.2])
    neutral = np.full(8, 0.1)
    matrix = pd.DataFrame({"stable": stable, "reversing": reversing, "neutral": neutral})
    matrix = matrix.div(matrix.sum(axis=1), axis=0)
    fitted = DomainStableCLRFeaturePipeline(min_prevalence=0.0, top_k=2, random_state=1).fit(matrix, labels, domains)
    assert "stable" in fitted.selected_features_
    assert "reversing" not in fitted.selected_features_
