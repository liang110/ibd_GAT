"""Summarize sample- versus participant-unit performance from frozen predictions.

This script never refits a model. It is intended for the reviewer-requested
unit-of-evaluation audit on the frozen six-direction benchmark predictions.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


def metric_row(y: np.ndarray, p: np.ndarray) -> dict[str, float]:
    out = {"auroc": np.nan, "pr_auc": np.nan, "brier": np.nan}
    if np.unique(y).size == 2:
        out["auroc"] = float(roc_auc_score(y, p))
        out["pr_auc"] = float(average_precision_score(y, p))
    out["brier"] = float(brier_score_loss(y, p))
    return out


def participant_table(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    conflicts = int(
        frame.groupby("bootstrap_cluster")["outcome_ibd"].nunique().gt(1).sum()
    )
    grouped = (
        frame.groupby("bootstrap_cluster", as_index=False)
        .agg(outcome_ibd=("outcome_ibd", "first"), probability=("probability", "mean"))
    )
    return grouped, conflicts


def bootstrap_ci(frame: pd.DataFrame, n_boot: int, seed: int) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    clusters = frame["bootstrap_cluster"].drop_duplicates().to_numpy()
    draws = {"auroc": [], "pr_auc": [], "brier": []}
    for _ in range(n_boot):
        sampled = rng.choice(clusters, size=len(clusters), replace=True)
        draw = frame.set_index("bootstrap_cluster").loc[sampled].reset_index()
        values = metric_row(draw["outcome_ibd"].to_numpy(), draw["probability"].to_numpy())
        for key in draws:
            draws[key].append(values[key])
    out: dict[str, float] = {}
    for key, values in draws.items():
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        out[f"{key}_ci_low"] = float(np.quantile(finite, 0.025)) if finite.size else np.nan
        out[f"{key}_ci_high"] = float(np.quantile(finite, 0.975)) if finite.size else np.nan
        out[f"{key}_bootstrap_valid"] = int(finite.size)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260819)
    parser.add_argument(
        "--direction",
        action="append",
        default=None,
        help="Directed transfer to retain; repeat for multiple directions. Defaults to all directions.",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.predictions)
    required = {"direction", "split", "model", "sample_id", "outcome_ibd", "probability", "bootstrap_cluster"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing prediction columns: {sorted(missing)}")
    df["outcome_ibd"] = df["outcome_ibd"].astype(int)
    df["probability"] = df["probability"].astype(float)
    if args.direction:
        requested = sorted(set(args.direction))
        available = set(df["direction"].unique())
        unknown = sorted(set(requested).difference(available))
        if unknown:
            raise ValueError(f"Requested directions absent from predictions: {unknown}")
        df = df.loc[df["direction"].isin(requested)].copy()
    args.outdir.mkdir(parents=True, exist_ok=True)

    split_rows: list[dict[str, object]] = []
    for (direction, split, model), sub in df.groupby(["direction", "split", "model"], sort=True):
        part, conflicts = participant_table(sub)
        sample_metrics = metric_row(sub["outcome_ibd"].to_numpy(), sub["probability"].to_numpy())
        participant_metrics = metric_row(part["outcome_ibd"].to_numpy(), part["probability"].to_numpy())
        row = {"direction": direction, "split": int(split), "model": model,
               "n_samples": len(sub), "n_participants": len(part),
               "participant_label_conflicts": conflicts}
        row.update({f"sample_{k}": v for k, v in sample_metrics.items()})
        row.update({f"participant_{k}": v for k, v in participant_metrics.items()})
        split_rows.append(row)

    split_out = pd.DataFrame(split_rows)
    split_out.to_csv(args.outdir / "participant_unit_by_split.csv", index=False)

    ensemble_rows: list[dict[str, object]] = []
    for (direction, model), sub in df.groupby(["direction", "model"], sort=True):
        # Equal-weight source perturbation ensemble for each target sample.
        sample_ensemble = (
            sub.groupby(["sample_id", "bootstrap_cluster"], as_index=False)
            .agg(outcome_ibd=("outcome_ibd", "first"), probability=("probability", "mean"))
        )
        part_ensemble, conflicts = participant_table(sample_ensemble)
        sample_metrics = metric_row(sample_ensemble["outcome_ibd"].to_numpy(), sample_ensemble["probability"].to_numpy())
        participant_metrics = metric_row(part_ensemble["outcome_ibd"].to_numpy(), part_ensemble["probability"].to_numpy())
        stable_key = f"{direction}|{model}".encode("utf-8")
        stable_offset = int(hashlib.sha256(stable_key).hexdigest()[:8], 16) % 100000
        ci = bootstrap_ci(part_ensemble, args.bootstrap, args.seed + stable_offset)
        row = {"direction": direction, "model": model,
               "n_samples": len(sample_ensemble), "n_participants": len(part_ensemble),
               "participant_label_conflicts": conflicts,
               "source_split_count": int(sub["split"].nunique()),
               "source_split_sample_auroc_mean": float(split_out.query("direction == @direction and model == @model")["sample_auroc"].mean()),
               "source_split_participant_auroc_mean": float(split_out.query("direction == @direction and model == @model")["participant_auroc"].mean())}
        row.update({f"ensemble_sample_{k}": v for k, v in sample_metrics.items()})
        row.update({f"ensemble_participant_{k}": v for k, v in participant_metrics.items()})
        row.update(ci)
        ensemble_rows.append(row)

    pd.DataFrame(ensemble_rows).to_csv(args.outdir / "participant_unit_ensemble_summary.csv", index=False)
    summary = {
        "input": str(args.predictions),
        "rows": int(len(df)),
        "directions": sorted(df["direction"].unique().tolist()),
        "requested_directions": args.direction,
        "models": sorted(df["model"].unique().tolist()),
        "bootstrap_replicates": int(args.bootstrap),
        "bootstrap_unit": "target participant cluster",
        "label_conflict_rows": int(sum(r["participant_label_conflicts"] for r in ensemble_rows)),
    }
    (args.outdir / "participant_unit_sensitivity_summary.json").write_text(
        pd.Series(summary).to_json(indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
