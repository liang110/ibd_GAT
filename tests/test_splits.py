import pandas as pd

from ibd_microbiome_framework.evaluation.splits import repeated_stratified_train_validation


def test_repeated_stratified_split_preserves_outcomes():
    manifest = pd.DataFrame({"outcome_ibd": [0, 1] * 20, "origin": ["a"] * 40})
    for _, train_idx, validation_idx in repeated_stratified_train_validation(
        manifest, n_splits=3, validation_size=0.25, seed=11
    ):
        assert set(manifest.iloc[train_idx].outcome_ibd) == {0, 1}
        assert set(manifest.iloc[validation_idx].outcome_ibd) == {0, 1}
