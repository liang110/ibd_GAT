"""Score fixed Franzosa-trained GAT and RF predictions under HMP2 sampling rules."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


MODELS = ("no_prior_gat", "random_forest")
METRICS = ("roc_auc", "average_precision", "brier")
DIRECTION = "PRISM_to_HMP2"


def read_inputs(prediction_file, metadata_file, score_file):
    keys = ["direction", "split", "model", "sample_id", "participant_id", "outcome", "probability"]
    pred = pd.read_csv(prediction_file, usecols=keys)
    pred = pred.loc[pred.direction.eq(DIRECTION) & pred.model.isin(MODELS)].copy()
    if pred.isna().any().any() or pred.duplicated(["model", "split", "sample_id"]).any():
        raise ValueError("Missing values or duplicate predictions")
    if set(pred.split.unique()) != set(range(1, 11)) or not pred.probability.between(0, 1).all():
        raise ValueError("Expected ten fits with probabilities in [0, 1]")
    identity = pred[["sample_id", "participant_id", "outcome"]].drop_duplicates()
    identity = identity.set_index("sample_id", verify_integrity=True).sort_index()
    by_fit = pred.pivot(index="sample_id", columns=["model", "split"], values="probability")
    if by_fit.shape != (1626, 20) or by_fit.isna().any().any():
        raise ValueError("Incomplete paired HMP2 prediction panel")

    cols = ["External ID", "Participant ID", "site_sub_coll", "data_type", "reads_filtered", "diagnosis"]
    provider = pd.read_csv(metadata_file, usecols=cols)
    provider = provider.loc[provider.data_type.eq("metagenomics") & provider["External ID"].notna()]
    provider = provider.set_index("External ID", verify_integrity=True).loc[identity.index]
    if provider[cols[1:]].isna().any().any():
        raise ValueError("Incomplete profile annotations")
    if not np.array_equal(provider["Participant ID"], identity.participant_id):
        raise ValueError("Participant keys differ between predictions and provider metadata")
    diagnosis = provider.diagnosis.map({"CD": 1, "UC": 1, "nonIBD": 0})
    if not np.array_equal(diagnosis, identity.outcome):
        raise ValueError("Outcome mismatch with provider diagnosis")
    identity["collection"] = provider.site_sub_coll
    if identity.groupby("collection").participant_id.nunique().max() != 1:
        raise ValueError("Collection shared between participant identifiers")
    if identity.groupby("participant_id").outcome.nunique().max() != 1:
        raise ValueError("Inconsistent participant outcome")
    scored = pd.read_csv(score_file, sep="\t", header=None, usecols=[0], dtype=str).iloc[:, 0]
    deep = provider.index[provider.reads_filtered.ge(1_000_000)]
    if len(scored) != 1595 or set(scored) != set(deep) or scored.duplicated().any():
        raise ValueError("Published sample list does not match 1,595 depth-qualified target profiles")
    identity["published"] = identity.index.isin(scored)
    return identity, by_fit.loc[identity.index]


def cluster_draws(labels):
    """Reuse the archived participant-level RNG stream and sorted person IDs."""
    token = f"{DIRECTION}|pooled|participant".encode()
    seed = 20261112 + int.from_bytes(hashlib.sha256(token).digest()[:4], "big")
    rng = np.random.default_rng(seed)
    weights = []
    while len(weights) < 5000:
        counts = np.bincount(rng.integers(len(labels), size=len(labels)), minlength=len(labels))
        if len(np.unique(labels[counts > 0])) == 2:
            weights.append(counts)
    return np.array(weights)


def score_panel(labels, probability, person_rows, draws):
    """Batch weighted tie-group ranking scores; include the observed sample first."""
    order = np.argsort(-probability, kind="stable")
    y = labels[order]
    start = np.r_[0, np.flatnonzero(np.diff(probability[order]) != 0) + 1]
    sq_error = (y - probability[order]) ** 2
    ordered_person = person_rows[order]
    values = []
    for offset in range(0, len(draws), 128):
        w = draws[offset:offset + 128, ordered_person]
        total = np.add.reduceat(w, start, axis=1)
        pos = np.add.reduceat(w * y, start, axis=1)
        neg = total - pos
        cum_pos = np.cumsum(pos, axis=1)
        cum_total = np.cumsum(total, axis=1)
        npos = cum_pos[:, -1]
        nneg = cum_total[:, -1] - npos
        above = nneg[:, None] - (cum_total - cum_pos) + 0.5 * neg
        auroc = (pos * above).sum(axis=1) / (npos * nneg)
        precision = np.divide(cum_pos, cum_total, out=np.zeros(cum_pos.shape), where=cum_total > 0)
        ap = (precision * pos).sum(axis=1) / npos
        brier = (w * sq_error).sum(axis=1) / cum_total[:, -1]
        values.append(np.column_stack((auroc, ap, brier)))
    result = np.vstack(values)
    for index in [0, 1, 37]:
        sample_weight = draws[index, person_rows]
        ref = [metric(labels, probability, sample_weight=sample_weight)
               for metric in (roc_auc_score, average_precision_score, brier_score_loss)]
        np.testing.assert_allclose(result[index], ref, rtol=0, atol=2e-12)
    return result


def run(args):
    info, predictions = read_inputs(args.predictions, args.metadata, args.published_scores)
    people = info.groupby("participant_id").outcome.first()
    draws = np.vstack((np.ones(len(people), dtype=int), cluster_draws(people.to_numpy())))
    rows = []
    for panel_name, mask in (("all_eligible", np.ones(len(info), dtype=bool)), ("published_1595", info.published)):
        meta = info.loc[mask]
        if set(meta.participant_id) != set(people.index) or (len(people), people.sum()) != (130, 103):
            raise ValueError("Expected identical 103 IBD and 27 non-IBD people in both panels")
        values = predictions.loc[meta.index]
        collections = values.groupby([meta.participant_id, meta.collection], sort=True).mean()
        for unit, block in (
            ("record", values),
            ("person_profile_mean", values.groupby(meta.participant_id, sort=True).mean()),
            ("person_collection_mean", collections.groupby(level=0, sort=True).mean()),
        ):
            ids = meta.participant_id if unit == "record" else block.index
            person_rows = people.index.get_indexer(ids)
            labels = people.to_numpy()[person_rows].astype(int)
            scores = {}
            for model in MODELS:
                scores[model] = np.mean(
                    [score_panel(labels, block[model, split].to_numpy(), person_rows, draws)
                     for split in range(1, 11)], axis=0)
            diff = scores[MODELS[0]] - scores[MODELS[1]]
            lo, hi = np.quantile(diff[1:], [0.025, 0.975], axis=0)
            for i, metric in enumerate(METRICS):
                rows.append(dict(panel=panel_name, scoring=unit, profiles=len(meta),
                                 collections=len(collections), people=len(people),
                                 ibd_profiles=int(meta.outcome.sum()), ibd_people=int(people.sum()),
                                 fits=10, draws=len(draws)-1, metric=metric,
                                 gat=scores[MODELS[0]][0, i], rf=scores[MODELS[1]][0, i],
                                 observed_delta=diff[0, i], draw_mean=diff[1:, i].mean(),
                                 ci_lower=lo[i], ci_upper=hi[i]))
    out = pd.DataFrame(rows)
    if args.reference is not None:
        ref = pd.read_csv(args.reference)
        ref = ref.loc[ref.direction.eq(DIRECTION) & ref.target_stratum.eq("pooled")
                      & ref.candidate.eq(MODELS[0]) & ref.reference.eq(MODELS[1])]
        for row in out.loc[out.panel.eq("all_eligible") & out.scoring.isin(["record", "person_profile_mean"])].itertuples():
            level = "visit" if row.scoring == "record" else "participant"
            sub = ref.loc[ref.analysis_level.eq(level) & ref.metric.eq(row.metric)]
            if len(sub) != 1:
                raise ValueError("Missing unique baseline contrast")
            r = sub.iloc[0]
            np.testing.assert_allclose(row.observed_delta, r.split_mean_point_delta, atol=1e-12, rtol=0)
            if level == "participant":
                np.testing.assert_allclose([row.draw_mean, row.ci_lower, row.ci_upper],
                                           [r.common_draw_bootstrap_mean, r.ci_low, r.ci_high], atol=1e-12, rtol=0)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True, help="Primary runner CSV or frozen long ledger.")
    parser.add_argument("--metadata", type=Path, required=True, help="IBDMDB metadata release of 20 August 2018.")
    parser.add_argument("--published-scores", type=Path, required=True, help="Headerless IBDMDB dysbiosis score TSV.")
    parser.add_argument("--reference", type=Path, help="Published primary_common_draw_summary.csv for the invariant checks.")
    parser.add_argument("--output", type=Path, required=True, help="New CSV file for the six scoring settings.")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a new output file")
    table = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="") as handle:
        table.to_csv(handle, index=False)
    print(f"Wrote {len(table)} rows to {args.output}")


if __name__ == "__main__":
    main()
