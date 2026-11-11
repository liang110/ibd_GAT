#!/usr/bin/env python
"""Evaluate a synthetic-screened domain-aware model on clean cohort transfer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import (
    repeated_group_train_validation,
    repeated_stratified_train_validation,
)
from ibd_microbiome_framework.models.domain_interaction import DomainInteractionLogisticClassifier
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


def _metrics(y: pd.Series, probability: np.ndarray) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(y, probability)),
        "average_precision": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }


def _domain_support(manifest: pd.DataFrame, domain_column: str) -> pd.DataFrame:
    support = manifest.groupby([domain_column, "outcome_ibd"]).size().rename("n").reset_index()
    eligible = support.groupby(domain_column).outcome_ibd.nunique().eq(2).rename("eligible")
    return support.merge(eligible, on=domain_column, validate="many_to_one")


def _apply_domain_policy(
    manifest: pd.DataFrame, domain_column: str, policy: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    support = _domain_support(manifest, domain_column)
    invalid = support.loc[~support.eligible, domain_column].unique().tolist()
    if invalid and policy == "strict":
        raise ValueError(f"Source domains without both outcomes: {invalid}")
    filtered = manifest.loc[manifest[domain_column].isin(support.loc[support.eligible, domain_column])].copy()
    if filtered.empty:
        raise ValueError("No source samples remain after domain-support filtering")
    return filtered.reset_index(drop=True), support


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-dir", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--domain-policy", choices=("strict", "exclude_unsupported"), default="strict")
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    hmp_manifest = pd.read_csv(args.contract_dir / "cohort_manifest_hmp2.csv")
    prism_manifest = pd.read_csv(args.contract_dir / "cohort_manifest_prism.csv")
    shared = load_shared_feature_list(args.contract_dir)
    hmp = load_canonical_microbiome(args.legacy_root / "manual_abundance.csv", hmp_manifest, shared)
    prism = load_canonical_microbiome(args.legacy_root / "prism_abundance_clean.csv", prism_manifest, shared)
    hmp_manifest = hmp_manifest.set_index("sample_id", drop=False).loc[hmp.index].reset_index(drop=True)
    prism_manifest = prism_manifest.set_index("sample_id", drop=False).loc[prism.index].reset_index(drop=True)
    cohorts = [
        ("HMP2_to_PRISM", hmp, hmp_manifest, prism, prism_manifest, "site", True),
        ("PRISM_to_HMP2", prism, prism_manifest, hmp, hmp_manifest, "origin", False),
    ]
    metrics, predictions, support_rows = [], [], []
    for direction, source, source_manifest, target, target_manifest, domain_column, grouped in cohorts:
        source_manifest, support = _apply_domain_policy(
            source_manifest, domain_column, args.domain_policy
        )
        support_rows.extend(
            {
                "direction": direction,
                "scope": "source_cohort",
                "split": None,
                "domain_column": domain_column,
                **row,
            }
            for row in support.rename(columns={domain_column: "domain"}).to_dict(orient="records")
        )
        source = source.loc[source_manifest.sample_id]
        if grouped:
            split_iterator = repeated_group_train_validation(source_manifest, args.n_splits, 0.2, 20260816)
        else:
            split_iterator = repeated_stratified_train_validation(source_manifest, args.n_splits, 0.2, 20260816)
        for split_id, train_idx, validation_idx in split_iterator:
            train_manifest = source_manifest.iloc[train_idx]
            validation_manifest = source_manifest.iloc[validation_idx]
            # A participant-group split can leave a rare site with only one
            # outcome in its training portion. Exclude that site from fitting
            # rather than allowing its site-specific deviation to absorb class.
            train_manifest, train_support = _apply_domain_policy(
                train_manifest, domain_column, "exclude_unsupported"
            )
            support_rows.extend(
                {
                    "direction": direction,
                    "scope": "source_train",
                    "split": split_id,
                    "domain_column": domain_column,
                    **row,
                }
                for row in train_support.rename(columns={domain_column: "domain"}).to_dict(orient="records")
            )
            train_ids = train_manifest.sample_id.tolist()
            validation_ids = validation_manifest.sample_id.tolist()
            pipeline = CLRFeaturePipeline(
                min_prevalence=0.10, top_k=100, random_state=20260816 + split_id
            ).fit(source.loc[train_ids], train_manifest.outcome_ibd)
            x_train = pipeline.transform(source.loc[train_ids])
            x_validation = pipeline.transform(source.loc[validation_ids])
            x_target = pipeline.transform(target)
            candidates = []
            for c in (0.01, 0.1, 1.0):
                model = DomainInteractionLogisticClassifier(c=c, random_state=20260816 + split_id).fit(
                    x_train, train_manifest.outcome_ibd, train_manifest[domain_column]
                )
                # Validation emulates an unseen domain: no site-specific
                # deviations are available at deployment.
                validation_probability = model.predict_proba(x_validation)[:, 1]
                candidates.append((roc_auc_score(validation_manifest.outcome_ibd, validation_probability), model, c))
            _, model, selected_c = max(candidates, key=lambda item: item[0])
            validation_probability = model.predict_proba(x_validation)[:, 1]
            target_probability = model.predict_proba(x_target)[:, 1]
            for evaluation, y_eval, probability in (
                ("source_validation", validation_manifest.outcome_ibd, validation_probability),
                ("external_target", target_manifest.outcome_ibd, target_probability),
            ):
                metrics.append({
                    "direction": direction,
                    "split": split_id,
                    "model": "domain_interaction_logistic",
                    "evaluation": evaluation,
                    "selected_c": selected_c,
                    "n_train": len(train_manifest),
                    "n_validation": len(validation_manifest),
                    "n_target": len(target_manifest),
                    "selected_features": pipeline.audit_.n_selected_features,
                    **_metrics(y_eval, probability),
                })
            for sample_id, outcome, probability, cluster, origin in zip(
                target.index,
                target_manifest.outcome_ibd,
                target_probability,
                target_manifest.participant_id,
                target_manifest.origin,
            ):
                predictions.append({
                    "direction": direction,
                    "split": split_id,
                    "model": "domain_interaction_logistic",
                    "sample_id": sample_id,
                    "outcome_ibd": int(outcome),
                    "probability": float(probability),
                    "bootstrap_cluster": str(cluster),
                    "origin": origin,
                })
    metric_frame = pd.DataFrame(metrics)
    metric_frame.to_csv(out / "domain_interaction_transfer_metrics.csv", index=False)
    pd.DataFrame(predictions).to_csv(out / "domain_interaction_transfer_predictions.csv", index=False)
    pd.DataFrame(support_rows).to_csv(out / "domain_interaction_source_domain_support.csv", index=False)
    summary = metric_frame.loc[metric_frame.evaluation == "external_target"].groupby("direction")[["roc_auc", "average_precision", "brier"]].agg(["mean", "std"]).round(6)
    summary.columns = [f"{metric}_{stat}" for metric, stat in summary.columns]
    report = {"n_splits": args.n_splits, "shared_features": len(shared), "external_target_summary": summary.reset_index().to_dict(orient="records")}
    (out / "domain_interaction_transfer_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
