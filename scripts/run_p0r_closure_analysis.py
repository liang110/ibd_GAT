#!/usr/bin/env python
"""Close P0-R with paired uncertainty, participant sensitivity, and calibration."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


METRICS = ("roc_auc", "average_precision", "brier")
CALIBRATION_METRICS = ("brier", "calibration_intercept", "calibration_slope", "ece_10")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _metric_values(y: np.ndarray, probability: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(y, probability, sample_weight=weights)),
        "average_precision": float(average_precision_score(y, probability, sample_weight=weights)),
        "brier": float(brier_score_loss(y, probability, sample_weight=weights)),
    }


def _calibration_values(
    y: np.ndarray, probability: np.ndarray, weights: np.ndarray | None = None, n_bins: int = 10
) -> dict[str, float]:
    probability = np.clip(np.asarray(probability, dtype=float), 1e-8, 1 - 1e-8)
    y = np.asarray(y, dtype=int)
    if weights is None:
        weights = np.ones(len(y), dtype=float)
    weights = np.asarray(weights, dtype=float)
    logit = np.log(probability / (1 - probability)).reshape(-1, 1)
    try:
        model = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000)
        model.fit(logit, y, sample_weight=weights)
        intercept = float(model.intercept_[0])
        slope = float(model.coef_[0, 0])
    except Exception:
        intercept = float("nan")
        slope = float("nan")
    edges = np.linspace(0, 1, n_bins + 1)
    bins = np.clip(np.digitize(probability, edges[1:-1], right=False), 0, n_bins - 1)
    ece = 0.0
    total_weight = float(weights.sum())
    for index in range(n_bins):
        mask = bins == index
        if not mask.any():
            continue
        bin_weight = float(weights[mask].sum())
        if bin_weight <= 0:
            continue
        mean_probability = float(np.average(probability[mask], weights=weights[mask]))
        observed = float(np.average(y[mask], weights=weights[mask]))
        ece += bin_weight / total_weight * abs(mean_probability - observed)
    return {
        "brier": float(brier_score_loss(y, probability, sample_weight=weights)),
        "calibration_intercept": intercept,
        "calibration_slope": slope,
        "ece_10": float(ece),
    }


def _cluster_weights(clusters: np.ndarray, labels: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    unique = pd.unique(clusters)
    for _ in range(100):
        sampled = rng.choice(unique, size=len(unique), replace=True)
        counts = pd.Series(sampled).value_counts()
        weights = pd.Series(clusters).map(counts).fillna(0).to_numpy(dtype=float)
        if np.unique(labels[weights > 0]).size == 2:
            return weights
    raise RuntimeError("Cluster bootstrap failed to retain both outcomes")


def _sorted_panel(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.sort_values("sample_id").reset_index(drop=True)


def _direction_panels(predictions: pd.DataFrame, models: list[str]) -> dict[str, dict[int, dict[str, pd.DataFrame]]]:
    required = {"direction", "split", "model", "sample_id", "outcome_ibd", "probability", "bootstrap_cluster"}
    missing = required.difference(predictions.columns)
    if missing:
        raise ValueError(f"Prediction table missing columns: {sorted(missing)}")
    panels: dict[str, dict[int, dict[str, pd.DataFrame]]] = {}
    for direction, direction_frame in predictions.groupby("direction", sort=True):
        available = set(direction_frame.model)
        unavailable = sorted(set(models).difference(available))
        if unavailable:
            raise ValueError(f"{direction} lacks required models: {unavailable}")
        split_panels: dict[int, dict[str, pd.DataFrame]] = {}
        for split, split_frame in direction_frame.groupby("split", sort=True):
            model_panels = {model: _sorted_panel(split_frame.loc[split_frame.model == model].copy()) for model in models}
            reference = model_panels[models[0]][["sample_id", "outcome_ibd", "bootstrap_cluster"]]
            if reference.sample_id.duplicated().any():
                raise ValueError(f"{direction} split {split} has duplicate sample IDs")
            for model, panel in model_panels.items():
                candidate = panel[["sample_id", "outcome_ibd", "bootstrap_cluster"]]
                if not reference.equals(candidate):
                    raise ValueError(f"{direction} split {split} model {model} does not share the locked target panel")
            if reference.groupby("bootstrap_cluster").outcome_ibd.nunique().gt(1).any():
                raise ValueError(f"{direction} split {split} has outcome-varying bootstrap clusters")
            split_panels[int(split)] = model_panels
        baseline = split_panels[sorted(split_panels)[0]][models[0]][["sample_id", "outcome_ibd", "bootstrap_cluster"]]
        for split, model_panels in split_panels.items():
            candidate = model_panels[models[0]][["sample_id", "outcome_ibd", "bootstrap_cluster"]]
            if not baseline.equals(candidate):
                raise ValueError(f"{direction} target panel changed across source perturbations at split {split}")
        panels[direction] = split_panels
    return panels


def _target_counts(panels: dict[str, dict[int, dict[str, pd.DataFrame]]], model: str) -> pd.DataFrame:
    rows = []
    for direction, split_panels in panels.items():
        panel = split_panels[sorted(split_panels)[0]][model]
        cluster_labels = panel.groupby("bootstrap_cluster").outcome_ibd.first()
        rows.append(
            {
                "direction": direction,
                "n_target_samples": int(len(panel)),
                "n_target_clusters": int(cluster_labels.size),
                "n_ibd_samples": int(panel.outcome_ibd.sum()),
                "n_control_samples": int((1 - panel.outcome_ibd).sum()),
                "n_ibd_clusters": int(cluster_labels.sum()),
                "n_control_clusters": int((1 - cluster_labels).sum()),
                "n_samples_in_repeated_clusters": int(panel.bootstrap_cluster.duplicated(keep=False).sum()),
                "max_samples_per_cluster": int(panel.bootstrap_cluster.value_counts().max()),
            }
        )
    return pd.DataFrame(rows)


def _paired_deltas(
    panels: dict[str, dict[int, dict[str, pd.DataFrame]]], models: list[str], n_bootstrap: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries: list[dict] = []
    draws: list[dict] = []
    for direction, split_panels in panels.items():
        split_ids = sorted(split_panels)
        reference_panel = split_panels[split_ids[0]][models[0]]
        labels = reference_panel.outcome_ibd.to_numpy()
        clusters = reference_panel.bootstrap_cluster.to_numpy()
        observed: dict[tuple[str, str], dict[str, list[float]]] = {}
        for candidate, reference in itertools.combinations(models, 2):
            observed[(candidate, reference)] = {metric: [] for metric in METRICS}
            for split in split_ids:
                candidate_frame = split_panels[split][candidate]
                reference_frame = split_panels[split][reference]
                candidate_metrics = _metric_values(labels, candidate_frame.probability.to_numpy())
                reference_metrics = _metric_values(labels, reference_frame.probability.to_numpy())
                for metric in METRICS:
                    observed[(candidate, reference)][metric].append(candidate_metrics[metric] - reference_metrics[metric])
        bootstrap_values: dict[tuple[str, str], dict[str, list[float]]] = {
            pair: {metric: [] for metric in METRICS} for pair in observed
        }
        rng = np.random.default_rng(seed + sum(ord(character) for character in direction))
        for replicate in range(n_bootstrap):
            weights = _cluster_weights(clusters, labels, rng)
            for candidate, reference in observed:
                deltas = {metric: [] for metric in METRICS}
                for split in split_ids:
                    candidate_metrics = _metric_values(labels, split_panels[split][candidate].probability.to_numpy(), weights)
                    reference_metrics = _metric_values(labels, split_panels[split][reference].probability.to_numpy(), weights)
                    for metric in METRICS:
                        deltas[metric].append(candidate_metrics[metric] - reference_metrics[metric])
                draw = {"direction": direction, "candidate": candidate, "reference": reference, "replicate": replicate}
                for metric in METRICS:
                    value = float(np.mean(deltas[metric]))
                    bootstrap_values[(candidate, reference)][metric].append(value)
                    draw[metric] = value
                draws.append(draw)
        for (candidate, reference), metric_map in observed.items():
            for metric, values in metric_map.items():
                boot = np.asarray(bootstrap_values[(candidate, reference)][metric], dtype=float)
                summaries.append(
                    {
                        "direction": direction,
                        "candidate": candidate,
                        "reference": reference,
                        "metric": metric,
                        "n_source_perturbations": len(split_ids),
                        "n_target_clusters": int(pd.unique(clusters).size),
                        "source_split_mean_delta": float(np.mean(values)),
                        "source_split_sd_delta": float(np.std(values, ddof=1)),
                        "source_split_min_delta": float(np.min(values)),
                        "source_split_max_delta": float(np.max(values)),
                        "target_cluster_bootstrap_mean_delta": float(np.mean(boot)),
                        "target_cluster_ci_low": float(np.quantile(boot, 0.025)),
                        "target_cluster_ci_high": float(np.quantile(boot, 0.975)),
                        "target_cluster_ci_excludes_zero": bool((np.quantile(boot, 0.025) > 0) or (np.quantile(boot, 0.975) < 0)),
                    }
                )
    return pd.DataFrame(summaries), pd.DataFrame(draws)


def _participant_level_hmp2(
    panels: dict[str, dict[int, dict[str, pd.DataFrame]]], models: list[str], n_bootstrap: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    by_split: list[dict] = []
    summaries: list[dict] = []
    draws: list[dict] = []
    for direction, split_panels in panels.items():
        if not direction.endswith("_to_HMP2"):
            continue
        split_ids = sorted(split_panels)
        for model in models:
            sample_metric_values = {metric: [] for metric in METRICS}
            participant_metric_values = {metric: [] for metric in METRICS}
            participant_frames: dict[int, pd.DataFrame] = {}
            for split in split_ids:
                frame = split_panels[split][model]
                participant = (
                    frame.groupby("bootstrap_cluster", as_index=False)
                    .agg(outcome_ibd=("outcome_ibd", "first"), probability=("probability", "mean"))
                    .sort_values("bootstrap_cluster")
                    .reset_index(drop=True)
                )
                participant_frames[split] = participant
                sample_metrics = _metric_values(frame.outcome_ibd.to_numpy(), frame.probability.to_numpy())
                participant_metrics = _metric_values(participant.outcome_ibd.to_numpy(), participant.probability.to_numpy())
                for metric in METRICS:
                    sample_metric_values[metric].append(sample_metrics[metric])
                    participant_metric_values[metric].append(participant_metrics[metric])
                    by_split.append(
                        {
                            "direction": direction,
                            "model": model,
                            "split": split,
                            "metric": metric,
                            "sample_level_value": sample_metrics[metric],
                            "participant_level_value": participant_metrics[metric],
                            "participant_minus_sample": participant_metrics[metric] - sample_metrics[metric],
                            "n_target_samples": len(frame),
                            "n_target_participants": len(participant),
                        }
                    )
            clusters = participant_frames[split_ids[0]].bootstrap_cluster.to_numpy()
            labels = participant_frames[split_ids[0]].outcome_ibd.to_numpy()
            rng = np.random.default_rng(seed + 1000 + sum(ord(character) for character in direction + model))
            bootstrap = {metric: [] for metric in METRICS}
            for replicate in range(n_bootstrap):
                weights = _cluster_weights(clusters, labels, rng)
                draw = {"direction": direction, "model": model, "replicate": replicate}
                for metric in METRICS:
                    values = []
                    for split in split_ids:
                        frame = participant_frames[split]
                        values.append(_metric_values(frame.outcome_ibd.to_numpy(), frame.probability.to_numpy(), weights)[metric])
                    value = float(np.mean(values))
                    bootstrap[metric].append(value)
                    draw[metric] = value
                draws.append(draw)
            for metric in METRICS:
                boot = np.asarray(bootstrap[metric], dtype=float)
                summaries.append(
                    {
                        "direction": direction,
                        "model": model,
                        "metric": metric,
                        "n_source_perturbations": len(split_ids),
                        "n_target_samples": int(len(split_panels[split_ids[0]][model])),
                        "n_target_participants": int(len(clusters)),
                        "sample_level_mean": float(np.mean(sample_metric_values[metric])),
                        "participant_level_mean": float(np.mean(participant_metric_values[metric])),
                        "participant_minus_sample_mean": float(np.mean(participant_metric_values[metric]) - np.mean(sample_metric_values[metric])),
                        "participant_level_source_split_sd": float(np.std(participant_metric_values[metric], ddof=1)),
                        "participant_level_target_cluster_ci_low": float(np.quantile(boot, 0.025)),
                        "participant_level_target_cluster_ci_high": float(np.quantile(boot, 0.975)),
                    }
                )
    return pd.DataFrame(by_split), pd.DataFrame(summaries), pd.DataFrame(draws)


def _reliability_table(frame: pd.DataFrame, n_bootstrap: int, rng: np.random.Generator) -> pd.DataFrame:
    ensemble = frame.groupby(["sample_id", "outcome_ibd", "bootstrap_cluster"], as_index=False).probability.mean()
    ensemble["bin"] = np.clip(np.digitize(ensemble.probability, np.linspace(0, 1, 11)[1:-1]), 0, 9)
    rows = []
    for bin_id, bin_frame in ensemble.groupby("bin", sort=True):
        rows.append(
            {
                "bin": int(bin_id),
                "n_target_samples": int(len(bin_frame)),
                "n_target_clusters": int(bin_frame.bootstrap_cluster.nunique()),
                "mean_probability": float(bin_frame.probability.mean()),
                "observed_frequency": float(bin_frame.outcome_ibd.mean()),
            }
        )
    result = pd.DataFrame(rows).set_index("bin")
    clusters = ensemble.bootstrap_cluster.to_numpy()
    labels = ensemble.outcome_ibd.to_numpy()
    bootstrap: dict[int, list[float]] = {int(index): [] for index in result.index}
    for _ in range(n_bootstrap):
        weights = _cluster_weights(clusters, labels, rng)
        weighted = ensemble.assign(weight=weights)
        for bin_id in bootstrap:
            subset = weighted.loc[weighted.bin == bin_id]
            if subset.weight.sum() > 0:
                bootstrap[bin_id].append(float(np.average(subset.outcome_ibd, weights=subset.weight)))
    for bin_id, values in bootstrap.items():
        if values:
            result.loc[bin_id, "observed_frequency_ci_low"] = float(np.quantile(values, 0.025))
            result.loc[bin_id, "observed_frequency_ci_high"] = float(np.quantile(values, 0.975))
    return result.reset_index()


def _calibration(
    panels: dict[str, dict[int, dict[str, pd.DataFrame]]], models: list[str], n_bootstrap: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    by_split: list[dict] = []
    summary: list[dict] = []
    reliability: list[dict] = []
    for direction, split_panels in panels.items():
        split_ids = sorted(split_panels)
        base = split_panels[split_ids[0]][models[0]]
        labels = base.outcome_ibd.to_numpy()
        clusters = base.bootstrap_cluster.to_numpy()
        for model in models:
            observed = {metric: [] for metric in CALIBRATION_METRICS}
            split_frames = []
            for split in split_ids:
                frame = split_panels[split][model]
                values = _calibration_values(frame.outcome_ibd.to_numpy(), frame.probability.to_numpy())
                by_split.append({"direction": direction, "model": model, "split": split, **values})
                for metric in CALIBRATION_METRICS:
                    observed[metric].append(values[metric])
                split_frames.append(frame[["sample_id", "outcome_ibd", "bootstrap_cluster", "probability"]].rename(columns={"probability": f"probability_{split}"}))
            # Treat source perturbations as repeated predictions of the same locked target,
            # not as independent target datasets.  The ensemble probability is the mean
            # over source perturbations; bootstrap uncertainty is then applied once at
            # the participant level to this frozen ensemble.
            ensemble = split_frames[0][["sample_id", "outcome_ibd", "bootstrap_cluster"]].copy()
            probability_columns = []
            for frame in split_frames:
                ensemble = ensemble.merge(frame[["sample_id", frame.columns[-1]]], on="sample_id", validate="one_to_one")
                probability_columns.append(frame.columns[-1])
            ensemble["probability"] = ensemble[probability_columns].mean(axis=1)
            ensemble = ensemble[["sample_id", "outcome_ibd", "bootstrap_cluster", "probability"]]
            ensemble_labels = ensemble.outcome_ibd.to_numpy()
            ensemble_probability = ensemble.probability.to_numpy()
            ensemble_clusters = ensemble.bootstrap_cluster.to_numpy()
            ensemble_observed = _calibration_values(ensemble_labels, ensemble_probability)
            rng = np.random.default_rng(seed + 2000 + sum(ord(character) for character in direction + model))
            boot = {metric: [] for metric in CALIBRATION_METRICS}
            for _ in range(n_bootstrap):
                weights = _cluster_weights(ensemble_clusters, ensemble_labels, rng)
                calibrated = _calibration_values(ensemble_labels, ensemble_probability, weights)
                for metric in CALIBRATION_METRICS:
                    boot[metric].append(calibrated[metric])
            for metric in CALIBRATION_METRICS:
                values = np.asarray(observed[metric], dtype=float)
                bootstrap = np.asarray(boot[metric], dtype=float)
                finite = bootstrap[np.isfinite(bootstrap)]
                summary.append(
                    {
                        "direction": direction,
                        "model": model,
                        "metric": metric,
                        "n_source_perturbations": len(split_ids),
                        "n_target_clusters": int(pd.unique(ensemble_clusters).size),
                        "source_split_mean": float(np.nanmean(values)),
                        "source_split_sd": float(np.nanstd(values, ddof=1)),
                        "source_ensemble_observed": float(ensemble_observed[metric]),
                        "target_cluster_bootstrap_mean": float(np.nanmean(bootstrap)),
                        "target_cluster_ci_low": float(np.quantile(finite, 0.025)) if finite.size else float("nan"),
                        "target_cluster_ci_high": float(np.quantile(finite, 0.975)) if finite.size else float("nan"),
                    }
                )
            table = _reliability_table(ensemble, n_bootstrap, np.random.default_rng(seed + 3000 + sum(ord(character) for character in direction + model)))
            reliability.extend(table.assign(direction=direction, model=model).to_dict(orient="records"))
    return pd.DataFrame(by_split), pd.DataFrame(summary), pd.DataFrame(reliability)


def _cohort_contract_summary(manifests: dict[str, Path], expected_features: int) -> pd.DataFrame:
    rows = []
    metadata_fields = ("age_years", "sex", "treatment_observed", "site", "timepoint", "origin")
    for cohort, path in manifests.items():
        frame = pd.read_csv(path)
        if frame.sample_id.duplicated().any():
            raise ValueError(f"{cohort} manifest has duplicate sample IDs")
        participants = frame.groupby("participant_id").outcome_ibd.first()
        availability = {
            field: int(frame[field].notna().sum()) if field in frame.columns else 0 for field in metadata_fields
        }
        rows.append(
            {
                "cohort": cohort,
                "manifest_path": str(path),
                "n_contract_samples": int(len(frame)),
                "n_contract_participants": int(frame.participant_id.nunique()),
                "n_ibd_samples": int(frame.outcome_ibd.sum()),
                "n_control_samples": int((1 - frame.outcome_ibd).sum()),
                "n_ibd_participants": int(participants.sum()),
                "n_control_participants": int((1 - participants).sum()),
                "n_samples_in_repeated_participants": int(frame.participant_id.duplicated(keep=False).sum()),
                "n_origins": int(frame.origin.nunique()) if "origin" in frame.columns else 0,
                "origins": json.dumps(frame.origin.value_counts().to_dict(), sort_keys=True) if "origin" in frame.columns else "{}",
                "expected_species_contract": expected_features,
                "metadata_nonmissing_counts": json.dumps(availability, sort_keys=True),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--hmp2-manifest", type=Path, required=True)
    parser.add_argument("--prism-manifest", type=Path, required=True)
    parser.add_argument("--nielsen-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument("--expected-features", type=int, default=126)
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    predictions = pd.read_csv(args.predictions)
    models = ["logistic_l2", "elastic_net", "linear_svm", "random_forest", "xgboost", "shallow_mlp"]
    panels = _direction_panels(predictions, models)

    _target_counts(panels, models[0]).to_csv(output / "p0r_direction_target_counts.csv", index=False)
    _cohort_contract_summary(
        {"HMP2": args.hmp2_manifest, "PRISM": args.prism_manifest, "NielsenHB_2014": args.nielsen_manifest},
        args.expected_features,
    ).to_csv(output / "p0r_cohort_contract_summary.csv", index=False)

    paired_summary, paired_draws = _paired_deltas(panels, models, args.n_bootstrap, args.seed)
    paired_summary.to_csv(output / "p0r_paired_model_deltas.csv", index=False)
    paired_draws.to_csv(output / "p0r_paired_model_delta_draws.csv", index=False)

    participant_by_split, participant_summary, participant_draws = _participant_level_hmp2(
        panels, models, args.n_bootstrap, args.seed
    )
    participant_by_split.to_csv(output / "p0r_hmp2_participant_level_by_split.csv", index=False)
    participant_summary.to_csv(output / "p0r_hmp2_participant_level_summary.csv", index=False)
    participant_draws.to_csv(output / "p0r_hmp2_participant_level_draws.csv", index=False)

    calibration_by_split, calibration_summary, reliability = _calibration(panels, models, args.n_bootstrap, args.seed)
    calibration_by_split.to_csv(output / "p0r_direction_calibration_by_split.csv", index=False)
    calibration_summary.to_csv(output / "p0r_direction_calibration_summary.csv", index=False)
    reliability.to_csv(output / "p0r_direction_calibration_reliability_bins.csv", index=False)

    manifest = {
        "protocol": "P0-R paired-target-cluster uncertainty, HMP2 participant sensitivity, and direction-level descriptive calibration",
        "input_predictions": str(args.predictions.resolve()),
        "input_predictions_sha256": _sha256(args.predictions),
        "manifest_sha256": {name: _sha256(path) for name, path in {"hmp2": args.hmp2_manifest, "prism": args.prism_manifest, "nielsen": args.nielsen_manifest}.items()},
        "models": models,
        "n_directions": len(panels),
        "n_source_perturbations_per_direction": sorted({len(value) for value in panels.values()}),
        "n_bootstrap": args.n_bootstrap,
        "seed": args.seed,
        "target_labels_used_for_refitting_or_selection": False,
        "calibration_interpretation": "locked-target descriptive evaluation only; calibration parameters are not reused for training",
    }
    (output / "p0r_analysis_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
