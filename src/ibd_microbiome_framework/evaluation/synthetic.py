"""Known-truth semi-synthetic portability benchmark for compositional data."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class SyntheticTruth:
    disease_features: tuple[str, ...]
    batch_features: tuple[str, ...]
    topology_features: tuple[str, ...]
    source_samples: int
    target_samples: int
    random_state: int


def _renormalize(values: np.ndarray) -> np.ndarray:
    return values / values.sum(axis=1, keepdims=True)


def make_portability_benchmark(
    background: pd.DataFrame,
    n_source: int = 180,
    n_target: int = 180,
    n_disease_features: int = 6,
    n_batch_features: int = 6,
    n_topology_features: int = 0,
    random_state: int = 20260815,
    return_domains: bool = False,
    disease_multiplier: float = 50.0,
    batch_multiplier: float = 2.0,
    target_topology_shift: bool = False,
    min_truth_prevalence: float = 0.0,
):
    """Inject shared disease effects and source-only batch effects into real profiles.

    Labels are synthetic. A disease-associated multiplicative effect is shared across
    source and target. The source contains two observed batch environments with
    opposite label-to-batch correlations; the target keeps one of those directions.
    A feature that is only a batch shortcut therefore has no pooled source association,
    while the disease effect has replication across all environments.
    """
    n_truth_features = n_disease_features + n_batch_features + n_topology_features
    if background.shape[1] < n_truth_features:
        raise ValueError("Background has too few features for the requested truth set")
    if len(background) < n_source + n_target:
        raise ValueError("Background has too few samples for the requested benchmark")
    if (background.sum(axis=1) <= 0).any():
        raise ValueError("Background contains zero-mass samples")

    rng = np.random.default_rng(random_state)
    if disease_multiplier <= 0 or batch_multiplier <= 0:
        raise ValueError("Synthetic multipliers must be positive")
    if not 0.0 <= min_truth_prevalence <= 1.0:
        raise ValueError("min_truth_prevalence must be between zero and one")
    eligible = background.columns[
        (background > 0).mean(axis=0) >= min_truth_prevalence
    ].to_numpy()
    if len(eligible) < n_disease_features + n_batch_features:
        raise ValueError("Too few features satisfy the truth prevalence threshold")
    selected = rng.choice(
        eligible, n_disease_features + n_batch_features, replace=False
    )
    disease = tuple(selected[:n_disease_features])
    batch = tuple(selected[n_disease_features : n_disease_features + n_batch_features])
    remaining = np.setdiff1d(background.columns.to_numpy(), selected, assume_unique=False)
    topology = tuple(rng.choice(remaining, n_topology_features, replace=False))
    sampled = background.sample(n=n_source + n_target, random_state=random_state).copy()
    values = sampled.to_numpy(dtype=float)
    values = _renormalize(values)

    labels = np.tile(np.array([0, 1], dtype=int), int(np.ceil(len(sampled) / 2)))[: len(sampled)]
    rng.shuffle(labels)
    source_y, target_y = labels[:n_source], labels[n_source:]
    source_values, target_values = values[:n_source].copy(), values[n_source:].copy()
    feature_index = {name: idx for idx, name in enumerate(background.columns)}
    disease_idx = [feature_index[name] for name in disease]
    batch_idx = [feature_index[name] for name in batch]
    topology_idx = [feature_index[name] for name in topology]

    # Shared causal effects: both cohorts carry the same abundance shifts.
    source_values[np.ix_(source_y == 1, disease_idx)] *= disease_multiplier
    target_values[np.ix_(target_y == 1, disease_idx)] *= disease_multiplier
    # Batch correlation reverses between two source environments. This is the minimum
    # controlled setting in which source data can distinguish replicated signal from a
    # site-specific shortcut; a single-source environment cannot identify that contrast.
    source_environment = np.arange(n_source) % 2
    source_values[np.ix_((source_y == 1) & (source_environment == 0), batch_idx)] *= batch_multiplier
    source_values[np.ix_((source_y == 0) & (source_environment == 1), batch_idx)] *= batch_multiplier
    target_values[np.ix_(target_y == 0, batch_idx)] *= batch_multiplier
    if target_topology_shift:
        if not topology_idx:
            raise ValueError("Topology shift requires at least one topology feature")
        # Independently permuting selected target columns preserves their marginal
        # distributions while disrupting their source-derived covariance structure.
        for feature_idx in topology_idx:
            target_values[:, feature_idx] = target_values[rng.permutation(n_target), feature_idx]
    source_values = _renormalize(source_values)
    target_values = _renormalize(target_values)

    source_ids = [f"synthetic_source_{i:04d}" for i in range(n_source)]
    target_ids = [f"synthetic_target_{i:04d}" for i in range(n_target)]
    source = pd.DataFrame(source_values, index=source_ids, columns=background.columns)
    target = pd.DataFrame(target_values, index=target_ids, columns=background.columns)
    truth = SyntheticTruth(disease, batch, topology, n_source, n_target, random_state)
    result = (source, pd.Series(source_y, index=source_ids, name="outcome_ibd"), target, pd.Series(target_y, index=target_ids, name="outcome_ibd"), truth)
    if not return_domains:
        return result
    source_domains = pd.Series(source_environment.astype(str), index=source_ids, name="source_environment")
    target_domains = pd.Series("target", index=target_ids, name="target_environment")
    return (*result, source_domains, target_domains)


def make_repeated_subject_benchmark(
    background: pd.DataFrame,
    n_subjects: int = 80,
    repeats: int = 3,
    random_state: int = 20260815,
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """Create near-replicate profiles with random subject-level labels.

    There is deliberately no disease effect. A sample-wise split can recognize
    held-out repeats from the same subject, whereas a subject-group split cannot
    use that identity shortcut.
    """
    if len(background) < n_subjects:
        raise ValueError("Background has too few profiles for repeated-subject simulation")
    if repeats < 2:
        raise ValueError("At least two repeats per subject are required")
    rng = np.random.default_rng(random_state)
    bases = background.sample(n=n_subjects, random_state=random_state).to_numpy(dtype=float)
    subject_labels = np.tile(np.array([0, 1], dtype=int), int(np.ceil(n_subjects / 2)))[:n_subjects]
    rng.shuffle(subject_labels)
    rows, labels, subjects, sample_ids = [], [], [], []
    for subject_idx, (base, label) in enumerate(zip(bases, subject_labels)):
        subject_id = f"subject_{subject_idx:04d}"
        for repeat_idx in range(repeats):
            noisy = base * rng.lognormal(mean=0.0, sigma=0.01, size=base.shape[0])
            rows.append(noisy / noisy.sum())
            labels.append(label)
            subjects.append(subject_id)
            sample_ids.append(f"{subject_id}_repeat_{repeat_idx}")
    matrix = pd.DataFrame(rows, index=sample_ids, columns=background.columns)
    return (
        matrix,
        pd.Series(labels, index=sample_ids, name="outcome_ibd"),
        pd.Series(subjects, index=sample_ids, name="participant_id"),
    )
