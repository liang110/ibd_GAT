#!/usr/bin/env python
"""Export the fixed shared-species contract for a native SIAMCAT reference run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _aligned_manifest(path: Path, matrix: pd.DataFrame) -> pd.DataFrame:
    manifest = pd.read_csv(path).set_index("sample_id", drop=False)
    aligned = manifest.loc[matrix.index].copy()
    required = {"outcome_ibd", "participant_id", "origin"}
    missing = sorted(required.difference(aligned.columns))
    if missing:
        raise ValueError(f"{path.name} lacks required fields: {missing}")
    return aligned.reset_index(drop=True)


def _write_cohort(name: str, matrix: pd.DataFrame, manifest: pd.DataFrame, output: Path) -> dict:
    # SIAMCAT expects features x samples. These values are closure-normalized in
    # the frozen 126-species space; SIAMCAT transforms are fitted only on source.
    matrix_path = output / f"{name}_shared_species.csv"
    source_metadata_path = output / f"{name}_source_metadata.csv"
    target_metadata_path = output / f"{name}_target_metadata.csv"
    locked_outcomes_path = output / f"{name}_locked_outcomes.csv"
    matrix.T.to_csv(matrix_path, index_label="feature_id")
    # The R training process receives outcome labels only when this cohort is a
    # source. Target labels live in a separate post-prediction-only file.
    manifest.loc[:, ["sample_id", "outcome_ibd", "participant_id", "origin"]].to_csv(
        source_metadata_path, index=False
    )
    manifest.loc[:, ["sample_id", "participant_id", "origin"]].to_csv(
        target_metadata_path, index=False
    )
    manifest.loc[:, ["sample_id", "outcome_ibd", "participant_id", "origin"]].to_csv(
        locked_outcomes_path, index=False
    )
    return {
        "cohort": name,
        "samples": int(len(matrix)),
        "ibd_samples": int(manifest.outcome_ibd.sum()),
        "control_samples": int((1 - manifest.outcome_ibd).sum()),
        "features": int(matrix.shape[1]),
        "matrix_file": matrix_path.name,
        "matrix_sha256": _sha256(matrix_path),
        "source_metadata_file": source_metadata_path.name,
        "source_metadata_sha256": _sha256(source_metadata_path),
        "target_metadata_file": target_metadata_path.name,
        "target_metadata_sha256": _sha256(target_metadata_path),
        "locked_outcomes_file": locked_outcomes_path.name,
        "locked_outcomes_sha256": _sha256(locked_outcomes_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-dir", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--nielsen-projection", type=Path, required=True)
    parser.add_argument("--nielsen-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    shared = load_shared_feature_list(args.contract_dir)
    hmp_manifest_raw = pd.read_csv(args.contract_dir / "cohort_manifest_hmp2.csv")
    prism_manifest_raw = pd.read_csv(args.contract_dir / "cohort_manifest_prism.csv")
    nielsen_manifest_raw = pd.read_csv(args.nielsen_manifest)
    hmp = load_canonical_microbiome(args.legacy_root / "manual_abundance.csv", hmp_manifest_raw, shared)
    prism = load_canonical_microbiome(
        args.legacy_root / "prism_abundance_clean.csv", prism_manifest_raw, shared
    )
    nielsen = load_canonical_microbiome(args.nielsen_projection, nielsen_manifest_raw, shared)
    summaries = [
        _write_cohort("hmp2", hmp, _aligned_manifest(args.contract_dir / "cohort_manifest_hmp2.csv", hmp), output),
        _write_cohort("prism", prism, _aligned_manifest(args.contract_dir / "cohort_manifest_prism.csv", prism), output),
        _write_cohort("nielsenhb_2014", nielsen, _aligned_manifest(args.nielsen_manifest, nielsen), output),
    ]
    report = {
        "purpose": "native_siamcat_reference_input_only",
        "feature_contract": "fixed_hmp2_prism_126_species",
        "shared_features": len(shared),
        "target_labels_available_to_r_training_process": False,
        "cohorts": summaries,
    }
    (output / "siamcat_input_manifest.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
