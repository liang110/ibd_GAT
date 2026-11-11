import numpy as np
import pandas as pd

from ibd_microbiome_framework.models.domain_interaction import DomainInteractionLogisticClassifier


def test_domain_interaction_predicts_unseen_domain_with_shared_block_only():
    x = pd.DataFrame({"a": [0.1, 0.2, 0.8, 0.9, 0.2, 0.1, 0.9, 0.8], "b": [0.9, 0.8, 0.2, 0.1, 0.8, 0.9, 0.1, 0.2]})
    y = [0, 0, 1, 1, 0, 0, 1, 1]
    domains = ["site_a"] * 4 + ["site_b"] * 4
    model = DomainInteractionLogisticClassifier(c=0.1, random_state=1).fit(x, y, domains)
    probability = model.predict_proba(x.iloc[:2])
    assert probability.shape == (2, 2)
    assert np.allclose(probability.sum(axis=1), 1.0)
    observed_domain_probability = model.predict_proba(x.iloc[:2], domains=["site_a", "site_a"])
    assert observed_domain_probability.shape == (2, 2)
