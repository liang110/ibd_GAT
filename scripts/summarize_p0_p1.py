#!/usr/bin/env python
from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[1] / "artifacts" / "p0_p1_20260818"
files = sorted(root.glob("contract_*/contract_sensitivity_aggregate.csv"))
frames = []
for path in files:
    frame = pd.read_csv(path)
    frame["variant_dir"] = path.parent.name
    frames.append(frame)
aggregate = pd.concat(frames, ignore_index=True)
print("P0 rows", aggregate.shape)
print("P0 variant mean deltas")
print(aggregate.groupby("variant_dir")[["delta_roc_auc_vs_registered", "delta_average_precision_vs_registered", "delta_brier_vs_registered"]].mean().round(4))
print("P0 max absolute AUROC delta")
print(aggregate.assign(abs_delta=aggregate.delta_roc_auc_vs_registered.abs()).groupby("variant_dir").abs_delta.max().round(4))
print("P0 best model by direction and variant")
print(aggregate.loc[aggregate.groupby(["variant_dir", "direction"]).roc_auc_mean.idxmax(), ["variant_dir", "direction", "model", "roc_auc_mean"]].to_string(index=False))

calibration = pd.read_csv(root / "calibration" / "calibration_summary.csv")
print("Calibration rows", calibration.shape)
print("Calibration mean Brier by fraction")
print(calibration.groupby("fraction").brier.mean().round(4))
print("Full-fraction calibration by model")
print(calibration[calibration.fraction == 1].groupby("model")[["brier", "ece_10", "calibration_slope", "calibration_intercept"]].mean().round(4))

feature_stability = pd.read_csv(root / "calibration" / "learning_curve_feature_stability.csv")
print("Full-fraction feature Jaccard by direction")
print(feature_stability[feature_stability.fraction == 1].groupby("direction").jaccard_mean.mean().round(4))
weight_stability = pd.read_csv(root / "calibration" / "learning_curve_weight_stability.csv")
print("Full-fraction weight Spearman by model")
print(weight_stability[weight_stability.fraction == 1].groupby("model").weight_spearman_mean.mean().round(4))

learning = pd.read_csv(root / "learning_curve_metrics.csv")
print("Learning curve AUROC means by fraction/model")
print(learning.groupby(["fraction", "model"]).roc_auc.mean().unstack(0).round(3).to_string())
