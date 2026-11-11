import numpy as np
import pandas as pd

from ibd_microbiome_framework.evaluation.synthetic import make_portability_benchmark


def test_synthetic_benchmark_is_compositional_and_reproducible():
    columns = [f"species_{i}" for i in range(20)]
    background = pd.DataFrame(np.full((80, 20), 1 / 20), columns=columns)
    first = make_portability_benchmark(background, n_source=30, n_target=30, random_state=7)
    second = make_portability_benchmark(background, n_source=30, n_target=30, random_state=7)
    source, source_y, target, target_y, truth = first
    assert np.allclose(source.sum(axis=1), 1.0)
    assert np.allclose(target.sum(axis=1), 1.0)
    assert source_y.equals(second[1]) and target_y.equals(second[3])
    assert truth == second[4]
    assert set(truth.disease_features).isdisjoint(truth.batch_features)


def test_topology_shift_requires_declared_features():
    columns = [f"species_{i}" for i in range(20)]
    background = pd.DataFrame(np.full((80, 20), 1 / 20), columns=columns)
    with np.testing.assert_raises_regex(ValueError, "Topology shift"):
        make_portability_benchmark(
            background,
            n_source=30,
            n_target=30,
            random_state=7,
            target_topology_shift=True,
        )


def test_truth_features_meet_requested_background_prevalence():
    columns = [f"species_{i}" for i in range(20)]
    background = pd.DataFrame(np.full((80, 20), 1 / 20), columns=columns)
    background.loc[:59, "species_19"] = 0.0
    background = background.div(background.sum(axis=1), axis=0)
    result = make_portability_benchmark(
        background,
        n_source=30,
        n_target=30,
        random_state=7,
        min_truth_prevalence=0.30,
    )
    truth = result[4]
    selected = set(truth.disease_features).union(truth.batch_features)
    assert "species_19" not in selected
