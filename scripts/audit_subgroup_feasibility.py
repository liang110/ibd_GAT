#!/usr/bin/env python
"""Audit participant-level feasibility for reviewer-requested IBD subgroups."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from ibd_microbiome_framework.data.matrices import (
    load_canonical_microbiome,
    load_shared_feature_list,
)


METADATA_ALIASES = {
    "diagnosis_subtype": ("diagnosis_subtype", "diagnosis", "group"),
    "activity": ("activity", "disease_activity", "hbi", "pms", "sccai"),
    "therapy": ("therapy", "medication", "treatment", "biologic"),
}


def _availability(manifest: pd.DataFrame, cohort: str) -> list[dict[str, object]]:
    rows = []
    for concept, aliases in METADATA_ALIASES.items():
        matches = [column for column in aliases if column in manifest.columns]
        rows.append(
            {
                "cohort": cohort,
                "concept": concept,
                "matching_columns": ";".join(matches),
                "available": bool(matches),
                "non_missing_samples": int(
                    manifest[matches].notna().any(axis=1).sum() if matches else 0
                ),
                "non_missing_participants": int(
                    manifest.loc[
                        manifest[matches].notna().any(axis=1), "participant_id"
                    ].nunique()
                    if matches
                    else 0
                ),
            }
        )
    return rows


def _subtype_rows(
    manifest: pd.DataFrame, cohort: str, min_per_class: int
) -> list[dict[str, object]]:
    rows = []
    strata_columns = ["cohort"]
    if cohort == "PRISM" and "origin" in manifest.columns:
        strata_columns.append("origin")
    for keys, stratum in manifest.groupby(strata_columns, dropna=False):
        key_values = keys if isinstance(keys, tuple) else (keys,)
        labels = dict(zip(strata_columns, key_values))
        controls = stratum.loc[stratum.outcome_ibd.eq(0)]
        for subtype in sorted(
            stratum.loc[stratum.outcome_ibd.eq(1), "diagnosis_subtype"].dropna().unique()
        ):
            cases = stratum.loc[
                stratum.outcome_ibd.eq(1) & stratum.diagnosis_subtype.eq(subtype)
            ]
            n_case = int(cases.participant_id.nunique())
            n_control = int(controls.participant_id.nunique())
            rows.append(
                {
                    **labels,
                    "analysis": "subtype_vs_control",
                    "subgroup": subtype,
                    "n_case_samples": int(len(cases)),
                    "n_case_participants": n_case,
                    "n_control_samples": int(len(controls)),
                    "n_control_participants": n_control,
                    "minimum_participants_per_class": min_per_class,
                    "eligible_for_descriptive_model": bool(
                        n_case >= min_per_class and n_control >= min_per_class
                    ),
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-dir", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-participants-per-class", type=int, default=20)
    args = parser.parse_args()
    if args.min_participants_per_class < 2:
        raise ValueError("Minimum participants per class must be at least two")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    manifests = {
        "HMP2": pd.read_csv(args.contract_dir / "cohort_manifest_hmp2.csv"),
        "PRISM": pd.read_csv(args.contract_dir / "cohort_manifest_prism.csv"),
    }
    shared_features = load_shared_feature_list(args.contract_dir)
    analyzable_ids = {
        "HMP2": load_canonical_microbiome(
            args.legacy_root / "manual_abundance.csv",
            manifests["HMP2"],
            shared_features,
        ).index,
        "PRISM": load_canonical_microbiome(
            args.legacy_root / "prism_abundance_clean.csv",
            manifests["PRISM"],
            shared_features,
        ).index,
    }
    availability, subgroup_rows = [], []
    for cohort, manifest in manifests.items():
        required = {"participant_id", "outcome_ibd", "diagnosis_subtype", "cohort"}
        missing = sorted(required.difference(manifest.columns))
        if missing:
            raise ValueError(f"{cohort} manifest lacks required fields: {missing}")
        manifest = (
            manifest.set_index("sample_id", drop=False)
            .loc[analyzable_ids[cohort]]
            .reset_index(drop=True)
        )
        availability.extend(_availability(manifest, cohort))
        subgroup_rows.extend(
            _subtype_rows(manifest, cohort, args.min_participants_per_class)
        )
    availability_frame = pd.DataFrame(availability)
    subgroup_frame = pd.DataFrame(subgroup_rows)
    availability_frame.to_csv(out / "subgroup_metadata_availability.csv", index=False)
    subgroup_frame.to_csv(out / "subtype_vs_control_feasibility.csv", index=False)
    eligible = subgroup_frame.loc[subgroup_frame.eligible_for_descriptive_model]
    report = {
        "minimum_participants_per_class": args.min_participants_per_class,
        "analyzable_samples": {
            cohort: int(len(sample_ids)) for cohort, sample_ids in analyzable_ids.items()
        },
        "eligible_subtype_comparisons": eligible.loc[
            :, ["cohort", "origin", "subgroup"]
        ].fillna("not_applicable").to_dict(orient="records"),
        "unavailable_metadata": availability_frame.loc[
            ~availability_frame.available, ["cohort", "concept"]
        ].to_dict(orient="records"),
        "interpretation": "Operational feasibility only; this threshold is not a power analysis.",
    }
    (out / "subgroup_feasibility_summary.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
