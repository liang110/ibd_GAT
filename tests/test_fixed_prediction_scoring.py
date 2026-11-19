"""Regression tests for the public prediction-to-score interface."""

import hashlib
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


path = Path(__file__).resolve().parents[1] / "scripts" / "score_fixed_predictions.py"
spec = importlib.util.spec_from_file_location("score_fixed_predictions", path)
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)


def example():
    target = pd.DataFrame(
        {
            "direction": "A_to_B",
            "sample_id": list("abcdef"),
            "participant_id": ["S2", "S2", "S1", "S3", "S4", "S4"],
            "outcome": [1, 1, 0, 1, 0, 0],
            "origin": ["X", "X", "X", "Y", "Y", "Y"],
        }
    )
    blocks = []
    for split, name, probs in [
        ("1", "gat", [0.5, 0.6, 0.5, 0.9, 0.4, 0.6]),
        ("2", "gat", [0.8, 0.3, 0.2, 0.5, 0.6, 0.8]),
        ("1", "rf", [0.6, 0.6, 0.5, 0.8, 0.4, 0.5]),
        ("2", "rf", [0.4, 0.3, 0.5, 0.7, 0.1, 0.9]),
    ]:
        blocks.append(target.drop(columns="origin").assign(split=split, model=name, probability=probs))
    return target, pd.concat(blocks, ignore_index=True)


@pytest.mark.parametrize("level", ["visit", "participant"])
def test_paired_intervals_match_sklearn_for_repeated_participants(level):
    target, pred = example()
    meta, matrix = scorer.aligned_predictions(pred.sample(frac=1, random_state=45), target, ["gat", "rf"])
    result = scorer.score_setting("A_to_B", "pooled", level, meta, matrix, ["gat", "rf"], 60, 45)
    ids = meta.participant_id
    y = meta.outcome
    if level == "participant":
        matrix = matrix.groupby(ids, sort=True).mean()
        y = y.groupby(ids, sort=True).first()
        ids = matrix.index
    participants = pd.unique(ids)
    offset = int.from_bytes(hashlib.sha256(f"A_to_B|pooled|{level}".encode()).digest()[:4], "big")
    rng = np.random.default_rng(45 + offset)
    draws = []
    while len(draws) < 60:
        sample = rng.choice(participants, size=len(participants), replace=True)
        weight = np.array([np.count_nonzero(sample == person) for person in ids])
        if np.unique(y[weight > 0]).size < 2:
            continue
        scores = {}
        for model in ("gat", "rf"):
            scores[model] = np.mean(
                [
                    [
                        f(y, matrix[split, model], sample_weight=weight)
                        for f in (roc_auc_score, average_precision_score, brier_score_loss)
                    ]
                    for split in ("1", "2")
                ],
                axis=0,
            )
        draws.append(scores["gat"] - scores["rf"])
    actual = pd.DataFrame(result)
    np.testing.assert_allclose(actual.bootstrap_mean, np.mean(draws, axis=0), atol=1e-14)
    np.testing.assert_allclose(actual.ci_low, np.quantile(draws, 0.025, axis=0), atol=1e-14)
    np.testing.assert_allclose(actual.ci_high, np.quantile(draws, 0.975, axis=0), atol=1e-14)
    assert set(actual.target_records) == {6} and set(actual.target_participants) == {4}


@pytest.mark.parametrize(
    "problem",
    [
        "duplicate",
        "missing_record",
        "missing_fit",
        "unknown_id",
        "wrong_outcome",
        "wrong_participant",
        "invalid_probability",
        "inconsistent_person",
    ],
)
def test_invalid_target_pairing_fails(problem):
    target, pred = example()
    if problem == "duplicate":
        pred = pd.concat([pred, pred.iloc[:1]])
    elif problem == "missing_record":
        pred = pred.loc[pred.sample_id.ne("c")]  # Missing in every model and every fit.
    elif problem == "missing_fit":
        pred = pred.loc[~(pred.split.eq("2") & pred.model.eq("gat"))]
    elif problem == "unknown_id":
        pred.loc[0, "sample_id"] = "other"
    elif problem == "wrong_outcome":
        pred.loc[0, "outcome"] = 0
    elif problem == "wrong_participant":
        pred.loc[0, "participant_id"] = "other"
    elif problem == "invalid_probability":
        pred.loc[0, "probability"] = 1.1
    else:
        target.loc[0, "outcome"] = 0
    with pytest.raises(ValueError):
        scorer.aligned_predictions(pred, target, ["gat", "rf"])
