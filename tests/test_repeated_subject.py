import numpy as np
import pandas as pd

from ibd_microbiome_framework.evaluation.synthetic import make_repeated_subject_benchmark


def test_repeated_subject_benchmark_has_subject_level_labels():
    background = pd.DataFrame(
        np.full((30, 10), 0.1), columns=[f"taxon_{i}" for i in range(10)]
    )
    matrix, labels, subjects = make_repeated_subject_benchmark(
        background, n_subjects=20, repeats=3, random_state=9
    )
    assert matrix.shape == (60, 10)
    assert np.allclose(matrix.sum(axis=1), 1.0)
    assert subjects.value_counts().eq(3).all()
    assert labels.groupby(subjects).nunique().eq(1).all()
