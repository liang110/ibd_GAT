"""Compare fixed classifier predictions across target populations and scoring units."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


PREDICTION_COLUMNS = ["direction", "split", "model", "sample_id", "participant_id", "outcome", "probability"]
TARGET_COLUMNS = ["direction", "sample_id", "participant_id", "outcome", "origin"]
METRICS = ("roc_auc", "average_precision", "brier")


def read_columns(path, columns):
    frame = pd.read_csv(
        path, usecols=columns, dtype={c: str for c in columns if c != "probability"}, keep_default_na=False
    )
    keys = [c for c in columns if c not in ("probability", "origin", "outcome")]
    if frame.empty or frame[keys].apply(lambda col: col.str.strip().eq("").any()).any():
        raise ValueError(f"{path}: empty data or missing identifier")
    if not frame.outcome.isin(["0", "1"]).all():
        raise ValueError(f"{path}: expected binary outcomes coded 0/1")
    frame["outcome"] = frame.outcome.astype(int)
    return frame


def study_manifest(prepared):
    from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list

    contract = prepared / "contract"
    taxa = load_shared_feature_list(contract)
    frames = []
    for direction, target, matrix_file in [
        ("HMP2_to_PRISM", "prism", "prism_abundance_clean.csv"),
        ("PRISM_to_HMP2", "hmp2", "manual_abundance.csv"),
    ]:
        manifest = pd.read_csv(
            contract / f"cohort_manifest_{target}.csv", dtype={"sample_id": str, "participant_id": str}
        )
        eligible = load_canonical_microbiome(prepared / matrix_file, manifest, taxa).index
        rows = manifest.set_index("sample_id").loc[eligible].rename_axis("sample_id").reset_index()
        rows = rows.rename(columns={"outcome_ibd": "outcome"}).assign(direction=direction)
        if target == "hmp2":
            rows["origin"] = ""  # HMP2 has no discovery/validation subdivision in this comparison.
        elif rows.origin.isna().any() or set(rows.origin) != {"discovery", "validation"}:
            raise ValueError("Expected complete discovery and validation origins for the Franzosa panel")
        frames.append(rows[TARGET_COLUMNS])
    return pd.concat(frames, ignore_index=True)


def aligned_predictions(frame, target, models):
    """Put every fit and procedure on exactly the manifest's sample and participant keys."""
    target = target.set_index("sample_id", verify_integrity=True).sort_index()
    if target[["participant_id", "outcome", "origin"]].isna().any().any():
        raise ValueError("Missing target metadata")
    if target.groupby("participant_id")[["outcome", "origin"]].nunique().gt(1).any().any():
        raise ValueError("Participant has inconsistent outcome or recruitment origin")
    origin = target.origin
    if origin.eq("pooled").any() or (origin.eq("").any() and origin.ne("").any()):
        raise ValueError(
            "Supply an origin for every target record or leave all origins blank; 'pooled' is reserved"
        )
    values = frame.loc[frame.model.isin(models)].copy()
    if values.empty or set(values.model) != set(models):
        raise ValueError("Missing requested classifier predictions")
    if not values.probability.between(0, 1).all():
        raise ValueError("Expected a finite probability in [0, 1]")
    if values.duplicated(["split", "model", "sample_id"]).any():
        raise ValueError("Duplicate fit/model/sample prediction")
    joined = values.join(target[["participant_id", "outcome"]], on="sample_id", rsuffix="_expected")
    if not (
        joined.participant_id.eq(joined.participant_id_expected) & joined.outcome.eq(joined.outcome_expected)
    ).all():
        raise ValueError(
            "A prediction has an unknown sample or disagrees with its manifest outcome or participant"
        )
    fits = sorted(values.split.unique())
    columns = pd.MultiIndex.from_product([fits, models], names=["split", "model"])
    wide = values.pivot(index="sample_id", columns=["split", "model"], values="probability").reindex(
        index=target.index, columns=columns
    )
    if wide.isna().any().any():
        raise ValueError("A target sample, model, or source fit is missing from the paired prediction panel")
    return target, wide


def bootstrap_counts(clusters, labels, count, seed):
    """Draw participants with replacement and discard single-class draws."""
    codes, names = pd.factorize(clusters, sort=False)
    positive = np.zeros(len(names), dtype=bool)
    positive[codes] = labels
    if np.unique(positive).size != 2:
        raise ValueError("Each requested target population must contain both outcome classes")
    rng = np.random.default_rng(seed)
    draws = np.empty((count, len(names)), dtype=np.int32)
    filled = rejected = 0
    while filled < count:
        sampled = rng.integers(0, len(names), size=len(names))
        weights = np.bincount(sampled, minlength=len(names))
        if not weights[positive].any() or not weights[~positive].any():
            rejected += 1
            if rejected > 20 * count + 100:
                raise ValueError("Target panel produced too few two-class bootstrap draws")
            continue
        draws[filled] = weights
        filled += 1
    return draws, codes, rejected


class ScoreOrder:
    """Cache tied score groups; only participant weights vary across draws."""

    def __init__(self, labels, probability):
        self.order = np.argsort(-probability, kind="stable")
        scores = probability[self.order]
        self.starts = np.r_[0, np.flatnonzero(np.diff(scores)) + 1]
        self.labels = labels[self.order]
        self.squared_error = (probability - labels) ** 2

    def metrics(self, weights):
        ordered = weights[:, self.order]
        total = np.add.reduceat(ordered, self.starts, axis=1)
        pos = np.add.reduceat(ordered * self.labels, self.starts, axis=1)
        neg = total - pos
        npos = pos.sum(axis=1)
        nneg = neg.sum(axis=1)
        lower_neg = nneg[:, None] - neg.cumsum(axis=1)
        auc = (pos * (lower_neg + 0.5 * neg)).sum(axis=1) / (npos * nneg)
        cumulative = total.cumsum(axis=1)
        precision = np.divide(
            pos.cumsum(axis=1), cumulative, out=np.zeros_like(cumulative, dtype=float), where=cumulative > 0
        )
        ap = (precision * pos).sum(axis=1) / npos
        brier = np.einsum("br,r->b", weights, self.squared_error) / (npos + nneg)
        return np.column_stack((auc, ap, brier))


def score_setting(direction, stratum, level, rows, wide, models, count, seed):
    selected = rows if stratum == "pooled" else rows.loc[rows.origin.eq(stratum)]
    scores = wide.loc[selected.index]
    if level == "participant":
        scores = scores.groupby(selected.participant_id, sort=True).mean()
        labels = selected.groupby("participant_id", sort=True).outcome.first().to_numpy()
        participants = scores.index
    else:
        labels = selected.outcome.to_numpy()
        participants = selected.participant_id
    token = "|".join([direction, stratum, level]).encode("utf-8")
    setting_seed = seed + int.from_bytes(hashlib.sha256(token).digest()[:4], "big")
    draws, codes, invalid = bootstrap_counts(participants, labels, count, setting_seed)
    fits = scores.columns.get_level_values("split").unique()
    observed = np.empty((len(models), len(fits), len(METRICS)))
    means = np.zeros((len(models), count, len(METRICS)))
    for m, model in enumerate(models):
        for s, split in enumerate(fits):
            order = ScoreOrder(labels, scores[split, model].to_numpy())
            observed[m, s] = order.metrics(np.ones((1, len(labels))))[0]
            for start in range(0, count, 64):
                batch = draws[start : start + 64, codes]
                means[m, start : start + len(batch)] += order.metrics(batch) / len(fits)
    results = []
    for m, reference in enumerate(models[1:], start=1):
        delta = means[0] - means[m]
        lo, hi = np.quantile(delta, [0.025, 0.975], axis=0)
        for k, metric in enumerate(METRICS):
            obs = observed[:, :, k]
            results.append(
                dict(
                    direction=direction,
                    target_stratum=stratum,
                    scoring_unit="record" if level == "visit" else level,
                    candidate=models[0],
                    reference=reference,
                    metric=metric,
                    candidate_mean=obs[0].mean(),
                    reference_mean=obs[m].mean(),
                    candidate_source_fit_sd=obs[0].std(ddof=1) if len(fits) > 1 else np.nan,
                    reference_source_fit_sd=obs[m].std(ddof=1) if len(fits) > 1 else np.nan,
                    observed_difference=(obs[0] - obs[m]).mean(),
                    bootstrap_mean=delta[:, k].mean(),
                    ci_low=lo[k],
                    ci_high=hi[k],
                    source_fits=len(fits),
                    target_records=len(selected),
                    target_participants=selected.participant_id.nunique(),
                    scored_units=len(labels),
                    bootstrap_draws=count,
                    invalid_draws=invalid,
                    setting_seed=setting_seed,
                )
            )
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, nargs="+", required=True)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "--prepared-dir", type=Path, help="Output directory of prepare_published_profiles.py."
    )
    inputs.add_argument(
        "--target-manifest",
        type=Path,
        help="Eligible records: direction,sample_id,participant_id,outcome,origin.",
    )
    parser.add_argument("--candidate", default="no_prior_gat")
    parser.add_argument("--references", nargs="+", default=["random_forest"])
    parser.add_argument("--direction", help="Score only this transfer direction.")
    parser.add_argument(
        "--expected-fits", type=int, help="Require this many distinct source fits in every direction."
    )
    parser.add_argument("--n-bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261112)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    models = [args.candidate, *args.references]
    if (
        len(models) != len(set(models))
        or args.n_bootstrap < 1
        or args.seed < 0
        or (args.expected_fits is not None and args.expected_fits < 1)
    ):
        parser.error("Provide distinct model names, positive counts and a non-negative seed")
    if args.output.exists():
        parser.error("Choose a fresh output filename")
    frames = [read_columns(path, PREDICTION_COLUMNS) for path in args.predictions]
    frame = pd.concat(frames, ignore_index=True)
    frame["probability"] = pd.to_numeric(frame.probability, errors="raise")
    manifest = (
        study_manifest(args.prepared_dir)
        if args.prepared_dir
        else read_columns(args.target_manifest, TARGET_COLUMNS)
    )
    if args.direction:
        frame = frame.loc[frame.direction.eq(args.direction)]
        manifest = manifest.loc[manifest.direction.eq(args.direction)]
    if manifest.empty or set(frame.direction) != set(manifest.direction):
        parser.error("Prediction and manifest directions must match")
    results = []
    for direction, targets in manifest.groupby("direction", sort=True):
        rows, wide = aligned_predictions(frame.loc[frame.direction.eq(direction)], targets, models)
        if (
            args.expected_fits is not None
            and wide.columns.get_level_values("split").nunique() != args.expected_fits
        ):
            parser.error(f"{direction}: expected {args.expected_fits} distinct source fits")
        for stratum in ["pooled", *sorted(set(rows.origin) - {""})]:
            for level in ("visit", "participant"):
                results.extend(
                    score_setting(direction, stratum, level, rows, wide, models, args.n_bootstrap, args.seed)
                )
    pd.DataFrame(results).to_csv(args.output, index=False, mode="x")
    print(f"Wrote {len(results)} paired contrasts to {args.output}")


if __name__ == "__main__":
    main()
