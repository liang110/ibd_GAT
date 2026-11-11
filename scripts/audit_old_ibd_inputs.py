#!/usr/bin/env python
"""Fail-closed Gate 0 audit for the legacy HMP2/PRISM inputs.

The audit is intentionally independent of model fitting.  It records hashes,
sample/participant/label integrity, clinical denylist hits, abundance scale,
and the exact shared feature contract before any experiment is allowed to run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from ibd_microbiome_framework.data.contracts import MatrixSpec, build_manifest, matrix_axes, sha256_file
from ibd_microbiome_framework.harmonization.features import canonical_taxon, is_forbidden_feature


def _specs(hmp2_root: Path, prism_root: Path) -> list[MatrixSpec]:
    hmp_meta = hmp2_root / "hmp2_metadata_clean.csv"
    if not hmp_meta.exists():
        hmp_meta = hmp2_root / "clean_metadata.csv"
    if not hmp_meta.exists():
        hmp_meta = hmp2_root / "hmp2_metadata_2018-08-20.csv"
    return [
        MatrixSpec("HMP2", "microbiome", hmp2_root / "manual_abundance.csv", "Species", hmp_meta),
        MatrixSpec("PRISM", "microbiome", prism_root / "prism_abundance_clean.csv", "Species", prism_root / "prism_metadata_clean.csv"),
    ]


def _matrix_audit(spec: MatrixSpec) -> tuple[pd.DataFrame, list[dict]]:
    samples, features = matrix_axes(spec.path)
    manifest = build_manifest(spec, samples)
    raw = pd.read_csv(spec.path, index_col=0, low_memory=False)
    raw = raw.loc[[not is_forbidden_feature(name) for name in raw.index]]
    numeric = raw.apply(pd.to_numeric, errors="coerce")
    invalid = int((raw.notna() & numeric.isna()).to_numpy().sum())
    sums = numeric.fillna(0).sum(axis=0)
    rows = [{"cohort": spec.cohort, "path": str(spec.path), "sha256": sha256_file(spec.path),
             "n_samples": len(samples), "n_features": len(features), "invalid_cells": invalid,
             "sum_min": float(sums.min()), "sum_median": float(sums.median()), "sum_max": float(sums.max()),
             "percent_like_fraction": bool(float(sums.median()) > 2.0)}]
    return manifest, rows


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--hmp2-root", type=Path, required=True)
    p.add_argument("--prism-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    specs = _specs(args.hmp2_root, args.prism_root)
    manifests, matrix_rows, crosswalk_rows = {}, [], []
    decision = "PASS"
    reasons: list[str] = []
    for spec in specs:
        if not spec.path.exists() or not spec.metadata_path.exists():
            reasons.append(f"missing input: {spec.path} or {spec.metadata_path}")
            continue
        try:
            manifest, rows = _matrix_audit(spec)
            manifests[spec.cohort] = manifest
            matrix_rows.extend(rows)
            for feature in matrix_axes(spec.path)[1]:
                crosswalk_rows.append({"cohort": spec.cohort, "source_feature_id": feature,
                    "canonical_feature_id": canonical_taxon(feature), "forbidden": is_forbidden_feature(feature)})
            if rows[0]["invalid_cells"]:
                reasons.append(f"non-numeric cells in {spec.path}")
        except Exception as exc:
            reasons.append(f"{spec.cohort}: {type(exc).__name__}: {exc}")
    if len(manifests) == 2:
        h = set(pd.DataFrame(crosswalk_rows).query("cohort == 'HMP2' and not forbidden").canonical_feature_id.dropna())
        q = set(pd.DataFrame(crosswalk_rows).query("cohort == 'PRISM' and not forbidden").canonical_feature_id.dropna())
        shared = sorted(h & q)
        pd.DataFrame({"canonical_feature_id": shared}).to_csv(out / "shared_feature_contract.csv", index=False)
        if not shared: reasons.append("no shared canonical microbiome features")
        for cohort, manifest in manifests.items():
            manifest.to_csv(out / f"sample_participant_label_{cohort.lower()}.csv", index=False)
            if manifest.sample_id.duplicated().any() or manifest.participant_id.isna().any():
                reasons.append(f"sample/participant integrity failed for {cohort}")
    else:
        reasons.append("both HMP2 and PRISM inputs are required")
    if reasons: decision = "FAIL"
    pd.DataFrame(matrix_rows).to_csv(out / "input_scale_audit.csv", index=False)
    pd.DataFrame(crosswalk_rows).to_csv(out / "feature_crosswalk.csv", index=False)
    (out / "input_sha256.tsv").write_text("\n".join(f"{r['cohort']}\t{r['sha256']}\t{r['path']}" for r in matrix_rows), encoding="utf-8")
    (out / "gate0_decision.json").write_text(json.dumps({"decision": decision, "reasons": reasons,
        "n_shared_microbiome_features": len(shared) if len(manifests) == 2 else 0,
        "clinical_denylist_enforced": True, "source_only_required": True}, indent=2), encoding="utf-8")
    if decision != "PASS": raise SystemExit("Gate 0 failed: " + "; ".join(reasons))


if __name__ == "__main__": main()
