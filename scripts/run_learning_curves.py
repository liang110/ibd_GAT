#!/usr/bin/env python
"""Source-only learning curves for the six registered transfer directions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation, repeated_stratified_train_validation
from ibd_microbiome_framework.evaluation.standard_benchmark import fit_final_source_models, select_standard_models

SEED = 20260818
FRACTIONS = (0.25, 0.50, 0.75, 1.00)


def _metric(y, p):
    return {"roc_auc": float(roc_auc_score(y, p)), "average_precision": float(average_precision_score(y, p)), "brier": float(brier_score_loss(y, p))}


def _manifest(path, matrix):
    m = pd.read_csv(path).set_index("sample_id", drop=False)
    return m.loc[matrix.index].reset_index(drop=True)


def _splits(m, grouped, n):
    if grouped:
        return repeated_group_train_validation(m, n, 0.20, SEED)
    return repeated_stratified_train_validation(m, n, 0.20, SEED)


def _fraction_indices(manifest, indices, fraction, grouped, seed):
    if fraction >= 0.999:
        return np.asarray(indices, dtype=int)
    sub = manifest.iloc[indices].copy()
    rng = np.random.default_rng(seed)
    if grouped:
        groups = sub.groupby("participant_id").outcome_ibd.agg(lambda s: int(s.mean() >= 0.5))
        chosen = []
        for outcome in (0, 1):
            ids = groups[groups == outcome].index.to_numpy()
            n_keep = max(2, int(round(len(ids) * fraction)))
            chosen.extend(rng.choice(ids, size=min(n_keep, len(ids)), replace=False).tolist())
        keep = sub.participant_id.isin(chosen).to_numpy()
    else:
        strata = sub.outcome_ibd.astype(str) + "__" + sub.origin.astype(str)
        keep = np.zeros(len(sub), dtype=bool)
        for level in strata.unique():
            loc = np.flatnonzero(strata.to_numpy() == level)
            n_keep = max(2, int(round(len(loc) * fraction)))
            keep[rng.choice(loc, size=min(n_keep, len(loc)), replace=False)] = True
    selected = np.asarray(indices, dtype=int)[keep]
    if np.unique(manifest.iloc[selected].outcome_ibd).size != 2:
        return np.asarray(indices, dtype=int)
    return selected


def _run_direction(direction, source, source_m, target, target_m, grouped, n_splits, fractions):
    rows, preds, features, weights = [], [], [], []
    for split, train_idx, valid_idx in _splits(source_m, grouped, n_splits):
        for fraction in fractions:
            fit_idx = _fraction_indices(source_m, train_idx, fraction, grouped, SEED + split * 100 + int(fraction * 100))
            train = source.iloc[fit_idx]
            y = source_m.iloc[fit_idx].outcome_ibd.to_numpy()
            groups = source_m.iloc[fit_idx].participant_id.to_numpy() if grouped else None
            selected = select_standard_models(train, y, groups, SEED + split, n_inner_splits=3)
            pipeline, fitted = fit_final_source_models(source.iloc[fit_idx], y, selected, SEED + split)
            x_target = pipeline.transform(target)
            base = {"direction": direction, "split": split, "fraction": fraction, "n_source_fit": len(fit_idx), "n_source_total": len(source), "n_selected_features": pipeline.audit_.n_selected_features}
            for rank, feature in enumerate(pipeline.selected_features_, 1):
                features.append({**base, "feature_rank": rank, "feature_id": feature})
            for model_name, model in fitted.items():
                p = model.predict_proba(x_target)[:, 1]
                sel = selected[model_name]
                rows.append({**base, "model": model_name, "inner_auc_mean": sel.inner_auc_mean, "inner_auc_std": sel.inner_auc_std, **_metric(target_m.outcome_ibd, p)})
                if hasattr(model, "coef_"):
                    values = np.asarray(model.coef_).reshape(-1)
                elif hasattr(model, "feature_importances_"):
                    values = np.asarray(model.feature_importances_).reshape(-1)
                elif hasattr(model, "coefs_"):
                    values = np.linalg.norm(np.asarray(model.coefs_[0]), axis=1)
                else:
                    values = None
                if values is not None:
                    for feature, value in zip(pipeline.selected_features_, values):
                        weights.append({**base, "model": model_name, "feature_id": feature, "feature_weight": float(value)})
                for sid, outcome, score, cluster, origin in zip(target.index, target_m.outcome_ibd, p, target_m.participant_id, target_m.origin):
                    preds.append({"direction": direction, "split": split, "fraction": fraction, "model": model_name, "sample_id": sid, "outcome_ibd": int(outcome), "probability": float(score), "bootstrap_cluster": str(cluster), "origin": str(origin)})
    return rows, preds, features, weights


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract-dir", type=Path, required=True)
    ap.add_argument("--legacy-root", type=Path, required=True)
    ap.add_argument("--nielsen-projection", type=Path, required=True)
    ap.add_argument("--nielsen-manifest", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--n-splits", type=int, default=10)
    args = ap.parse_args(); args.output_dir.mkdir(parents=True, exist_ok=True)
    shared = load_shared_feature_list(args.contract_dir)
    hm = pd.read_csv(args.contract_dir / "cohort_manifest_hmp2.csv"); pm = pd.read_csv(args.contract_dir / "cohort_manifest_prism.csv"); nm = pd.read_csv(args.nielsen_manifest)
    h = load_canonical_microbiome(args.legacy_root / "manual_abundance.csv", hm, shared); p = load_canonical_microbiome(args.legacy_root / "prism_abundance_clean.csv", pm, shared); n = load_canonical_microbiome(args.nielsen_projection, nm, shared)
    hm = _manifest(args.contract_dir / "cohort_manifest_hmp2.csv", h); pm = _manifest(args.contract_dir / "cohort_manifest_prism.csv", p); nm = _manifest(args.nielsen_manifest, n)
    directions = [("HMP2_to_PRISM", h, hm, p, pm, True), ("PRISM_to_HMP2", p, pm, h, hm, False), ("HMP2_to_NielsenHB_2014", h, hm, n, nm, True), ("PRISM_to_NielsenHB_2014", p, pm, n, nm, False), ("NielsenHB_2014_to_HMP2", n, nm, h, hm, True), ("NielsenHB_2014_to_PRISM", n, nm, p, pm, True)]
    all_rows=[]; all_preds=[]; all_features=[]; all_weights=[]
    for values in directions:
        rows,preds,features,weights = _run_direction(*values, n_splits=args.n_splits, fractions=FRACTIONS)
        all_rows.extend(rows); all_preds.extend(preds); all_features.extend(features); all_weights.extend(weights)
    pd.DataFrame(all_rows).to_csv(args.output_dir / "learning_curve_metrics.csv", index=False)
    pd.DataFrame(all_preds).to_csv(args.output_dir / "learning_curve_predictions.csv", index=False)
    pd.DataFrame(all_features).to_csv(args.output_dir / "learning_curve_features.csv", index=False)
    pd.DataFrame(all_weights).to_csv(args.output_dir / "learning_curve_weights.csv", index=False)
    (args.output_dir / "learning_curve_summary.json").write_text(json.dumps({"fractions": FRACTIONS, "n_splits": args.n_splits, "directions": [x[0] for x in directions], "source_only": True, "target_labels_used_for_selection": False}, indent=2))
    print("learning curve complete")

if __name__ == "__main__":
    main()
