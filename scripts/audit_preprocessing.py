#!/usr/bin/env python
"""Build shared microbiome matrices and verify fitted-transform boundaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-dir", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-splits", type=int, default=5)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    contract = args.contract_dir.resolve()
    hmp2_manifest = pd.read_csv(contract / "cohort_manifest_hmp2.csv")
    prism_manifest = pd.read_csv(contract / "cohort_manifest_prism.csv")
    shared = load_shared_feature_list(contract)
    hmp2 = load_canonical_microbiome(args.legacy_root / "manual_abundance.csv", hmp2_manifest, shared)
    prism = load_canonical_microbiome(args.legacy_root / "prism_abundance_clean.csv", prism_manifest, shared)
    hmp2_excluded = hmp2_manifest.loc[~hmp2_manifest["sample_id"].isin(hmp2.index)].copy()
    prism_excluded = prism_manifest.loc[~prism_manifest["sample_id"].isin(prism.index)].copy()
    hmp2_manifest = (
        hmp2_manifest[hmp2_manifest["sample_id"].isin(hmp2.index)]
        .set_index("sample_id", drop=False)
        .loc[hmp2.index]
        .reset_index(drop=True)
    )
    prism_manifest = (
        prism_manifest[prism_manifest["sample_id"].isin(prism.index)]
        .set_index("sample_id", drop=False)
        .loc[prism.index]
        .reset_index(drop=True)
    )
    hmp2_excluded.to_csv(output / "hmp2_excluded_zero_shared_mass.csv", index=False)
    prism_excluded.to_csv(output / "prism_excluded_zero_shared_mass.csv", index=False)
    hmp2.to_csv(output / "hmp2_shared_microbiome_proportions.csv")
    prism.to_csv(output / "prism_shared_microbiome_proportions.csv")

    rows = []
    for split_id, train_idx, validation_idx in repeated_group_train_validation(hmp2_manifest, args.n_splits, 0.2, 42):
        train_ids = hmp2_manifest.iloc[train_idx]["sample_id"].tolist()
        validation_ids = hmp2_manifest.iloc[validation_idx]["sample_id"].tolist()
        pipeline = CLRFeaturePipeline(min_prevalence=0.10, top_k=100, random_state=42 + split_id)
        pipeline.fit(hmp2.loc[train_ids], hmp2_manifest.iloc[train_idx]["outcome_ibd"])
        train_x = pipeline.transform(hmp2.loc[train_ids])
        validation_x = pipeline.transform(hmp2.loc[validation_ids])
        target_x = pipeline.transform(prism)
        if list(train_x.columns) != list(validation_x.columns) or list(train_x.columns) != list(target_x.columns):
            raise AssertionError("Feature order diverged after transform")
        rows.append(
            {
                "split": split_id,
                **pipeline.audit_.__dict__,
                "n_validation_samples": len(validation_ids),
                "n_target_samples": len(prism),
                "target_origin_counts": json.dumps(prism_manifest["origin"].value_counts().to_dict(), sort_keys=True),
            }
        )
    pd.DataFrame(rows).to_csv(output / "preprocessing_split_audit.csv", index=False)
    summary = {
        "shared_microbiome_features": len(shared),
        "hmp2_samples": len(hmp2),
        "prism_samples": len(prism),
        "hmp2_excluded_zero_shared_mass": len(hmp2_excluded),
        "prism_excluded_zero_shared_mass": len(prism_excluded),
        "n_splits": args.n_splits,
        "metabolomics_in_primary_transfer": False,
    }
    (output / "preprocessing_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
