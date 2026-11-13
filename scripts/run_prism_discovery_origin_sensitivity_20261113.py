#!/usr/bin/env python
"""PRISM-discovery-only to HMP2 source-composition sensitivity.

This is a narrow sensitivity analysis for the pooled-source claim.  It uses
the reconciled 126-feature contract, source-training-only CLR/MI fitting and
the same participant-grouped split and model settings as the primary transfer
analysis.  It does not infer site or sequencing-run effects because those
fields are not part of the validated analysis contract.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation
from ibd_microbiome_framework.models.gat import NoPriorGAT
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


def build_graph(x: np.ndarray, mode: str, seed: int) -> torch.Tensor:
    n = x.shape[1]
    if mode == "identity":
        return torch.stack([torch.arange(n), torch.arange(n)])
    corr = pd.DataFrame(x).corr(method="spearman").to_numpy()
    np.fill_diagonal(corr, 0.0)
    edges: list[tuple[int, int]] = []
    k = min(8, n - 1)
    for i in range(n):
        for j in np.argsort(np.abs(corr[i]))[-k:]:
            if abs(corr[i, j]) > 0:
                edges.append((i, int(j)))
    if mode == "rewired":
        rng = np.random.default_rng(seed)
        edges = [(i, int((rng.integers(0, n - 1) + (rng.integers(0, n - 1) >= i)) % n)) for i, _ in edges]
    if not edges:
        return torch.empty((2, 0), dtype=torch.long)
    return torch.tensor(edges, dtype=torch.long).T.contiguous()


def graph_data(x: np.ndarray, y: np.ndarray, edge_index: torch.Tensor) -> list[Data]:
    return [Data(x=torch.tensor(row[:, None], dtype=torch.float32), y=torch.tensor([float(label)]), edge_index=edge_index)
            for row, label in zip(x, y)]


def fit_gat(x: np.ndarray, y: np.ndarray, edge_index: torch.Tensor, seed: int, epochs: int) -> NoPriorGAT:
    torch.manual_seed(seed)
    model = NoPriorGAT(add_self_loops=False)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-3)
    loader = DataLoader(graph_data(x, y, edge_index), batch_size=32, shuffle=True)
    model.train()
    for _ in range(epochs):
        for batch in loader:
            optimizer.zero_grad()
            logits = model(batch.x, batch.edge_index, batch.batch)
            loss = torch.nn.functional.binary_cross_entropy_with_logits(logits, batch.y.view(-1))
            loss.backward()
            optimizer.step()
    return model


def predict_gat(model: NoPriorGAT, x: np.ndarray, edge_index: torch.Tensor) -> np.ndarray:
    model.eval()
    loader = DataLoader(graph_data(x, np.zeros(len(x)), edge_index), batch_size=64)
    with torch.no_grad():
        return torch.cat([torch.sigmoid(model(batch.x, batch.edge_index, batch.batch)).cpu() for batch in loader]).numpy()


def metric_row(direction: str, split: int, model: str, probability: np.ndarray, y: np.ndarray, n_train: int, n_features: int, edge_index: torch.Tensor) -> dict:
    return {
        "direction": direction,
        "source_origin": "discovery",
        "target_origin": "HMP2",
        "split": int(split),
        "model": model,
        "roc_auc": float(roc_auc_score(y, probability)),
        "average_precision": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
        "n_train": int(n_train),
        "n_target": int(len(y)),
        "n_target_participants": int(len(y)),
        "n_features": int(n_features),
        "edge_count": int(edge_index.shape[1]),
        "edge_hash": hashlib.sha256(edge_index.numpy().tobytes()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-dir", type=Path, default=Path("artifacts/hpc_results/contract_v2_reconcile_20261112/data_contract"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/reconciled_submission_20261113/prism_discovery_origin_sensitivity"))
    parser.add_argument("--n-splits", type=int, default=10)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--seed", type=int, default=20261110)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(1)

    c = args.contract_dir
    features = pd.read_csv(c / "shared_feature_candidates.csv")
    features = features.loc[features.modality.eq("microbiome"), "canonical_feature_id"].astype(str).tolist()
    if len(features) != 126 or len(set(features)) != 126:
        raise ValueError("The sensitivity analysis requires exactly 126 unique shared microbial features")

    hmp = pd.read_csv(c / "hmp2_shared_microbiome_proportions.csv").rename(columns={"Unnamed: 0": "sample_id"}).set_index("sample_id")
    prism = pd.read_csv(c / "prism_shared_microbiome_proportions.csv").rename(columns={"Unnamed: 0": "sample_id"}).set_index("sample_id")
    hmp_manifest = pd.read_csv(c / "cohort_manifest_hmp2.csv").set_index("sample_id", drop=False)
    prism_manifest = pd.read_csv(c / "cohort_manifest_prism.csv").set_index("sample_id", drop=False)
    hmp_manifest = hmp_manifest.loc[hmp.index].copy()
    prism_manifest = prism_manifest.loc[prism.index].copy()
    prism_manifest = prism_manifest.loc[prism_manifest.origin.eq("discovery")].copy()
    hmp = hmp.loc[hmp_manifest.index, features]
    prism = prism.loc[prism_manifest.index, features]
    if not hmp.index.is_unique or not prism.index.is_unique:
        raise ValueError("Sample identifiers must be unique after contract alignment")
    if not set(prism_manifest.participant_id).issubset(set(prism_manifest.index)):
        # PRISM uses one operational participant per record; this check documents
        # the intended weighting unit without relying on a site field.
        raise ValueError("Unexpected PRISM participant mapping")

    rows: list[dict] = []
    predictions: list[dict] = []
    for split, train_idx, _ in repeated_group_train_validation(prism_manifest.reset_index(drop=True), args.n_splits, 0.2, args.seed):
        train_manifest = prism_manifest.reset_index(drop=True).iloc[train_idx]
        train_ids = train_manifest.sample_id.tolist()
        pipe = CLRFeaturePipeline(top_k=min(100, len(features)), random_state=args.seed + split).fit(prism.loc[train_ids], train_manifest.outcome_ibd.to_numpy())
        x_train = pipe.transform(prism.loc[train_ids]).to_numpy()
        x_target = pipe.transform(hmp).to_numpy()
        edge_index = build_graph(x_train, "real", args.seed + split)
        y_train = train_manifest.outcome_ibd.to_numpy()
        y_target = hmp_manifest.outcome_ibd.to_numpy()
        models = {
            "graph_attention": predict_gat(fit_gat(x_train, y_train, edge_index, args.seed + split, args.epochs), x_target, edge_index),
            "elastic_net": LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=.5, C=1.0, class_weight="balanced", max_iter=3000, random_state=args.seed + split).fit(x_train, y_train).predict_proba(x_target)[:, 1],
            "random_forest": RandomForestClassifier(n_estimators=300, min_samples_leaf=2, class_weight="balanced", random_state=args.seed + split, n_jobs=1).fit(x_train, y_train).predict_proba(x_target)[:, 1],
        }
        for name, prob in models.items():
            rows.append(metric_row("PRISM_discovery_to_HMP2", split, name, prob, y_target, len(train_idx), x_train.shape[1], edge_index))
            predictions.extend({"direction": "PRISM_discovery_to_HMP2", "split": split, "model": name, "sample_id": sid, "participant_id": pid, "outcome": int(y), "probability": float(p)} for sid, pid, y, p in zip(hmp_manifest.sample_id, hmp_manifest.participant_id, y_target, prob))
    metrics = pd.DataFrame(rows)
    metrics.to_csv(args.output_dir / "metrics.csv", index=False)
    pd.DataFrame(predictions).to_csv(args.output_dir / "predictions.csv", index=False)
    summary = metrics.groupby(["direction", "source_origin", "target_origin", "model"], as_index=False).agg(roc_auc_mean=("roc_auc", "mean"), roc_auc_sd=("roc_auc", "std"), roc_auc_min=("roc_auc", "min"), roc_auc_max=("roc_auc", "max"), average_precision_mean=("average_precision", "mean"), brier_mean=("brier", "mean"))
    summary.to_csv(args.output_dir / "summary.csv", index=False)
    manifest = {"analysis": "PRISM discovery-only source composition sensitivity", "source_origin": "discovery", "target_origin": "HMP2", "n_source": int(len(prism_manifest)), "n_target": int(len(hmp_manifest)), "n_source_participants": int(prism_manifest.participant_id.nunique()), "n_target_participants": int(hmp_manifest.participant_id.nunique()), "n_features": 126, "n_splits": args.n_splits, "epochs": args.epochs, "models": ["graph_attention", "elastic_net", "random_forest"], "source_only_preprocessing": True}
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
