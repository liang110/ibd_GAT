#!/usr/bin/env python
"""Post-fit Nielsen sensitivity analysis without target-side model selection."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


METRICS = ("roc_auc", "average_precision", "brier")


def _metrics(frame: pd.DataFrame, weights: np.ndarray | None = None) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(frame.outcome_ibd, frame.probability, sample_weight=weights)),
        "average_precision": float(
            average_precision_score(frame.outcome_ibd, frame.probability, sample_weight=weights)
        ),
        "brier": float(brier_score_loss(frame.outcome_ibd, frame.probability, sample_weight=weights)),
    }


def _cluster_weights(frame: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    clusters = frame.bootstrap_cluster.drop_duplicates().to_numpy()
    for _ in range(100):
        counts = pd.Series(rng.choice(clusters, size=len(clusters), replace=True)).value_counts()
        weights = frame.bootstrap_cluster.map(counts).fillna(0).to_numpy(dtype=float)
        if frame.loc[weights > 0, "outcome_ibd"].nunique() == 2:
            return weights
    raise RuntimeError("Cluster bootstrap could not retain both outcomes")


def _subgroup_definitions(manifest: pd.DataFrame) -> tuple[dict[str, pd.Index], list[dict[str, object]]]:
    indexed = manifest.set_index("sample_id", drop=False)
    definitions: dict[str, pd.Index] = {}
    infeasible: list[dict[str, object]] = []

    country = indexed["country"].fillna("unknown").astype(str)
    for value in sorted(country.unique()):
        sample_ids = indexed.loc[country.eq(value), "sample_id"]
        if indexed.loc[sample_ids, "outcome_ibd"].nunique() == 2:
            definitions[f"country_{value}"] = sample_ids
        else:
            infeasible.append(
                {
                    "sensitivity": f"country_{value}",
                    "reason": "one_outcome_only",
                    "n_samples": int(len(sample_ids)),
                    "outcome_counts": indexed.loc[sample_ids, "outcome_ibd"].value_counts().to_dict(),
                }
            )

    subtype = indexed["diagnosis_subtype"].fillna("Control").astype(str)
    for label in ("CD", "UC"):
        sample_ids = indexed.loc[subtype.isin(["Control", label]), "sample_id"]
        definitions[f"{label}_vs_control"] = sample_ids

    first = (
        indexed.reset_index(drop=True)
        .assign(_timepoint=lambda frame: pd.to_numeric(frame["timepoint"], errors="coerce"))
        .sort_values(["participant_id", "_timepoint", "sample_id"], na_position="last")
        .drop_duplicates("participant_id", keep="first")
    )
    definitions["first_sample_per_participant"] = first["sample_id"]
    definitions["age_complete"] = indexed.loc[indexed["age_years"].notna(), "sample_id"]
    return definitions, infeasible


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260816)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest)
    predictions = pd.read_csv(args.predictions)
    merged = predictions.merge(
        manifest,
        on=["sample_id", "outcome_ibd"],
        how="left",
        validate="many_to_one",
    )
    if merged.country.isna().any():
        raise ValueError("Predictions contain sample IDs absent from the target manifest")
    definitions, infeasible = _subgroup_definitions(manifest)
    metric_rows: list[dict[str, object]] = []
    draw_rows: list[dict[str, object]] = []
    for sensitivity, sample_ids in definitions.items():
        subset = merged.loc[merged.sample_id.isin(sample_ids)].copy()
        for group_index, (keys, frame) in enumerate(
            subset.groupby(["direction", "split", "model"], sort=True)
        ):
            if frame.outcome_ibd.nunique() != 2:
                raise ValueError(f"{sensitivity} unexpectedly lacks both outcomes")
            rng = np.random.default_rng(args.seed + group_index + 1009 * len(metric_rows))
            observed = _metrics(frame)
            boot = []
            for replicate in range(args.n_bootstrap):
                values = _metrics(frame, _cluster_weights(frame, rng))
                boot.append(values)
                draw_rows.append(
                    {
                        "sensitivity": sensitivity,
                        "direction": keys[0],
                        "split": keys[1],
                        "model": keys[2],
                        "replicate": replicate,
                        **values,
                    }
                )
            boot_frame = pd.DataFrame(boot)
            metric_rows.append(
                {
                    "sensitivity": sensitivity,
                    "direction": keys[0],
                    "split": keys[1],
                    "model": keys[2],
                    "n_samples": int(len(frame)),
                    "n_clusters": int(frame.bootstrap_cluster.nunique()),
                    "n_control": int((frame.outcome_ibd == 0).sum()),
                    "n_ibd": int((frame.outcome_ibd == 1).sum()),
                    **{f"{metric}_observed": value for metric, value in observed.items()},
                    **{f"{metric}_ci_low": float(boot_frame[metric].quantile(0.025)) for metric in METRICS},
                    **{f"{metric}_ci_high": float(boot_frame[metric].quantile(0.975)) for metric in METRICS},
                }
            )
    metrics = pd.DataFrame(metric_rows)
    draws = pd.DataFrame(draw_rows)
    metrics.to_csv(out / "nielsen_sensitivity_metrics.csv", index=False)
    draws.to_csv(out / "nielsen_sensitivity_bootstrap_draws.csv", index=False)
    aggregate_draws = draws.groupby(
        ["sensitivity", "direction", "model", "replicate"], as_index=False
    )[list(METRICS)].mean()
    aggregate_rows: list[dict[str, object]] = []
    observed = metrics.groupby(["sensitivity", "direction", "model"])[
        [f"{metric}_observed" for metric in METRICS]
    ].agg(["mean", "std"])
    observed.columns = [f"{metric}_{stat}" for metric, stat in observed.columns]
    for keys, frame in aggregate_draws.groupby(["sensitivity", "direction", "model"]):
        row = {"sensitivity": keys[0], "direction": keys[1], "model": keys[2]}
        for metric in METRICS:
            row[f"{metric}_cluster_ci_low"] = float(frame[metric].quantile(0.025))
            row[f"{metric}_cluster_ci_high"] = float(frame[metric].quantile(0.975))
        aggregate_rows.append(row)
    aggregate = observed.reset_index().merge(
        pd.DataFrame(aggregate_rows), on=["sensitivity", "direction", "model"], how="inner"
    )
    aggregate.to_csv(out / "nielsen_sensitivity_aggregate.csv", index=False)
    pd.DataFrame(infeasible).to_csv(out / "nielsen_sensitivity_infeasible.csv", index=False)
    report = {
        "protocol": "post_fit_target_sensitivity_no_target_model_selection",
        "n_bootstrap": args.n_bootstrap,
        "seed": args.seed,
        "active_sensitivities": sorted(definitions),
        "infeasible_sensitivities": infeasible,
        "aggregate": aggregate.round(6).to_dict(orient="records"),
    }
    (out / "nielsen_sensitivity_summary.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
