#!/usr/bin/env python
"""Analyze locked-target calibration and source-perturbation stability."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss
from scipy.stats import spearmanr


def _calibration(frame: pd.DataFrame, n_bins: int = 10) -> dict[str, float]:
    y = frame.outcome_ibd.to_numpy(float); p = frame.probability.to_numpy(float)
    eps = 1e-8
    x = np.clip(p, eps, 1 - eps)
    X = np.column_stack([np.ones_like(x), np.log(x / (1 - x))])
    try:
        fit = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000).fit(X[:, 1:].reshape(-1, 1), y)
        intercept = float(fit.intercept_[0]); slope = float(fit.coef_[0, 0])
    except Exception:
        intercept = float("nan"); slope = float("nan")
    prevalence = float(y.mean())
    brier = float(brier_score_loss(y, p))
    uncertainty = prevalence * (1 - prevalence)
    edges = np.linspace(0, 1, n_bins + 1)
    bins = np.clip(np.digitize(p, edges[1:-1], right=False), 0, n_bins - 1)
    ece = 0.0
    reliability = 0.0
    resolution = 0.0
    for idx in range(n_bins):
        mask = bins == idx
        if mask.any():
            weight = float(mask.mean())
            mean_probability = float(p[mask].mean())
            observed_frequency = float(y[mask].mean())
            ece += weight * abs(mean_probability - observed_frequency)
            reliability += weight * (mean_probability - observed_frequency) ** 2
            resolution += weight * (observed_frequency - prevalence) ** 2
    decomposed_brier = uncertainty - resolution + reliability
    return {"n_samples": len(frame), "n_clusters": int(frame.bootstrap_cluster.nunique()), "target_prevalence": prevalence, "brier": brier, "calibration_intercept": intercept, "calibration_slope": slope, "brier_uncertainty": uncertainty, "brier_resolution_binned": resolution, "brier_reliability_binned": reliability, "brier_decomposition_binned": decomposed_brier, "brier_decomposition_residual": brier - decomposed_brier, "ece_10": ece}


def _bootstrap_calibration(frame: pd.DataFrame, n_bootstrap: int, rng: np.random.Generator) -> dict[str, float]:
    clusters = frame.bootstrap_cluster.drop_duplicates().to_numpy(); rows=[]
    for _ in range(n_bootstrap):
        sampled = rng.choice(clusters, size=len(clusters), replace=True); counts = pd.Series(sampled).value_counts(); weights = frame.bootstrap_cluster.map(counts).fillna(0).to_numpy(float)
        if frame.loc[weights > 0, "outcome_ibd"].nunique() != 2: continue
        # Use positional indexing: prediction tables retain cohort-specific labels.
        expanded = frame.iloc[np.repeat(np.arange(len(frame)), weights.astype(int))].reset_index(drop=True)
        rows.append(_calibration(expanded))
    boot = pd.DataFrame(rows); result={}
    for col in ["brier", "calibration_intercept", "calibration_slope", "ece_10", "brier_reliability_binned", "brier_resolution_binned"]:
        result[f"{col}_ci_low"] = float(boot[col].quantile(0.025)); result[f"{col}_ci_high"] = float(boot[col].quantile(0.975))
    return result


def _feature_stability(features: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows=[]; frequency_rows=[]
    keys=["direction", "fraction"]
    for values, frame in features.groupby(keys, sort=True):
        sets={int(s):set(g.feature_id) for s,g in frame.groupby("split")}
        scores=[]; kuncheva=[]; universe=126
        for a,b in itertools.combinations(sets,2):
            left,right=sets[a],sets[b]; union=left|right
            if union: scores.append(len(left&right)/len(union))
            expected=len(left)*len(right)/universe
            denominator=min(len(left),len(right))-expected
            if denominator > 0: kuncheva.append((len(left&right)-expected)/denominator)
        freq=frame.groupby("feature_id").split.nunique()/len(sets)
        rows.append(dict(zip(keys,values)) | {"n_splits":len(sets), "jaccard_mean":float(np.mean(scores)), "jaccard_sd":float(np.std(scores,ddof=1)), "kuncheva_mean":float(np.mean(kuncheva)), "kuncheva_sd":float(np.std(kuncheva,ddof=1)), "n_features_any_split":int(freq.size), "n_features_frequency_ge_0p8":int((freq>=0.8).sum())})
        frequency_rows.append(freq.rename("selection_frequency").reset_index().assign(direction=values[0], fraction=values[1], n_splits=len(sets)))
    return pd.DataFrame(rows), pd.concat(frequency_rows, ignore_index=True)


def _weight_stability(weights: pd.DataFrame) -> pd.DataFrame:
    rows=[]; keys=["direction", "fraction", "model"]
    for values, frame in weights.groupby(keys, sort=True):
        mat=frame.pivot(index="split",columns="feature_id",values="feature_weight").fillna(0.0); corr=[]; signs=[]
        for a,b in itertools.combinations(mat.index,2):
            corr.append(float(spearmanr(mat.loc[a],mat.loc[b]).statistic)); active=(mat.loc[a]!=0)&(mat.loc[b]!=0)
            if active.any(): signs.append(float((np.sign(mat.loc[a,active])==np.sign(mat.loc[b,active])).mean()))
        rows.append(dict(zip(keys,values)) | {"weight_spearman_mean":float(np.nanmean(corr)), "weight_spearman_sd":float(np.nanstd(corr,ddof=1)), "shared_sign_agreement_mean":float(np.nanmean(signs)) if signs else float("nan")})
    return pd.DataFrame(rows)


def _rank_stability(metrics: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for (direction, fraction), frame in metrics.groupby(["direction", "fraction"], sort=True):
        ranks=frame.copy()
        ranks["rank"] = ranks.groupby("split").roc_auc.rank(ascending=False, method="average")
        winner_frequency = ranks.assign(is_first=ranks.rank == 1).groupby("model").is_first.mean()
        rows.append({"direction":direction, "fraction":fraction, "top_model_by_mean_auroc":frame.groupby("model").roc_auc.mean().idxmax(), "top_model_first_rank_frequency":float(winner_frequency.max()), "n_models_with_first_rank":int((winner_frequency > 0).sum()), "rank_sd_mean":float(ranks.groupby("model")["rank"].std(ddof=1).mean())})
    return pd.DataFrame(rows)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--predictions",type=Path,required=True); ap.add_argument("--features",type=Path,required=True); ap.add_argument("--weights",type=Path,required=True); ap.add_argument("--metrics",type=Path,required=True); ap.add_argument("--output-dir",type=Path,required=True); ap.add_argument("--n-bootstrap",type=int,default=500); args=ap.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
    pred=pd.read_csv(args.predictions); rng=np.random.default_rng(20260818); rows=[]; reliability=[]
    for keys, frame in pred.groupby(["direction","fraction","model"], sort=True):
        stats=_calibration(frame); stats.update(dict(zip(["direction","fraction","model"],keys))); stats.update(_bootstrap_calibration(frame,args.n_bootstrap,rng)); rows.append(stats)
        p=frame.probability.to_numpy(); y=frame.outcome_ibd.to_numpy(); bins=np.clip(np.digitize(p,np.linspace(0,1,11)[1:-1]),0,9)
        for b in range(10):
            mask=bins==b
            if mask.any(): reliability.append({"direction":keys[0],"fraction":keys[1],"model":keys[2],"bin":b,"n_samples":int(mask.sum()),"mean_probability":float(p[mask].mean()),"observed_frequency":float(y[mask].mean())})
    pd.DataFrame(rows).to_csv(args.output_dir/"calibration_summary.csv",index=False); pd.DataFrame(reliability).to_csv(args.output_dir/"calibration_reliability_bins.csv",index=False)
    feature_stability, selection_frequency = _feature_stability(pd.read_csv(args.features))
    feature_stability.to_csv(args.output_dir/"learning_curve_feature_stability.csv",index=False)
    selection_frequency.to_csv(args.output_dir/"learning_curve_feature_selection_frequency.csv",index=False)
    _weight_stability(pd.read_csv(args.weights)).to_csv(args.output_dir/"learning_curve_weight_stability.csv",index=False)
    _rank_stability(pd.read_csv(args.metrics)).to_csv(args.output_dir/"learning_curve_model_rank_stability.csv",index=False)
    summary={"protocol":"locked_target_descriptive_calibration","n_bootstrap":args.n_bootstrap,"target_calibration_parameters_not_reused_for_training":True,"metrics":["Brier","calibration_intercept","calibration_slope","ECE10","binned Brier reliability/resolution/decomposition","feature Jaccard and chance-corrected Kuncheva","target-evaluation model rank stability"]}
    (args.output_dir/"calibration_stability_summary.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=="__main__": main()
