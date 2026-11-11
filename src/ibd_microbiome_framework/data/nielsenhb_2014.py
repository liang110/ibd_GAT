"""Deterministically export the public NielsenHB_2014 RData resources.

The upstream ExperimentHub resources are immutable RData objects.  This
module deliberately performs no label recoding or feature harmonization: it
only turns the published objects into delimited inputs and records the sample
ID join needed for the later data contract.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Optional

import pandas as pd
import numpy as np

from ibd_microbiome_framework.harmonization.features import canonical_taxon


STUDY_NAME = "NielsenHB_2014"


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of a source resource."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(block_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_rdata(path: Path) -> Any:
    try:
        import rdata
    except ImportError as exc:  # pragma: no cover - depends on optional extra.
        raise RuntimeError(
            "The third-cohort importer requires the optional dependency 'rdata'. "
            "Install it with: pip install 'ibd-microbiome-framework[third-cohort]'"
        ) from exc
    loaded = rdata.read_rda(path)
    if len(loaded) != 1:
        raise ValueError(f"Expected one object in {path.name}, found {len(loaded)}")
    return next(iter(loaded.values()))


def export_public_resources(
    relative_abundance_rda: Path,
    sample_metadata_rda: Path,
    output_dir: Path,
) -> dict[str, object]:
    """Export the published taxa-by-sample matrix and curated study metadata.

    ``relative_abundance.csv`` preserves the source orientation (taxa rows,
    samples columns).  ``sample_metadata.csv`` contains every curated metadata
    row for this study; no rows are silently removed because downstream label
    eligibility must remain auditable.
    """
    abundance = _load_rdata(relative_abundance_rda)
    if not hasattr(abundance, "dims") or len(abundance.dims) != 2:
        raise ValueError("Relative-abundance object is not a two-dimensional array")
    feature_axis, sample_axis = abundance.dims
    matrix = pd.DataFrame(
        abundance.values,
        index=pd.Index(abundance.coords[feature_axis].values, name="feature_id"),
        columns=pd.Index(abundance.coords[sample_axis].values, name="sample_id"),
    )
    if matrix.index.duplicated().any() or matrix.columns.duplicated().any():
        raise ValueError("Published abundance object contains duplicate feature or sample IDs")

    metadata_object = _load_rdata(sample_metadata_rda)
    metadata = pd.DataFrame(metadata_object).copy()
    if "study_name" not in metadata.columns or "sample_id" not in metadata.columns:
        raise ValueError("Curated metadata lacks required study_name or sample_id columns")
    metadata = metadata.loc[metadata["study_name"].eq(STUDY_NAME)].copy()
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    if metadata["sample_id"].duplicated().any():
        raise ValueError("Curated metadata has duplicate sample IDs for NielsenHB_2014")

    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = output_dir / "relative_abundance.csv"
    metadata_path = output_dir / "sample_metadata.csv"
    matrix.to_csv(matrix_path, index=True)
    metadata.to_csv(metadata_path, index=False)

    matrix_ids = set(matrix.columns)
    metadata_ids = set(metadata["sample_id"])
    summary: dict[str, object] = {
        "study_name": STUDY_NAME,
        "relative_abundance_rda": str(relative_abundance_rda),
        "relative_abundance_rda_sha256": sha256_file(relative_abundance_rda),
        "sample_metadata_rda": str(sample_metadata_rda),
        "sample_metadata_rda_sha256": sha256_file(sample_metadata_rda),
        "matrix_shape_taxa_by_samples": [int(matrix.shape[0]), int(matrix.shape[1])],
        "n_metadata_rows_study": int(metadata.shape[0]),
        "n_matrix_only_sample_ids": int(len(matrix_ids - metadata_ids)),
        "n_metadata_only_sample_ids": int(len(metadata_ids - matrix_ids)),
        "n_shared_sample_ids": int(len(matrix_ids & metadata_ids)),
        "metadata_columns": metadata.columns.tolist(),
    }
    (output_dir / "export_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary


def terminal_species_label(feature_id: object) -> Optional[str]:
    """Map a MetaPhlAn lineage to a species only when it explicitly has ``s__``.

    This is intentionally separate from the framework's general canonicalizer:
    the source cohorts already supply species labels while NielsenHB_2014
    supplies complete lineages.  Collapsing a lineage to its terminal species
    is a documented representation conversion, not a feature-selection step.
    """
    species = [part for part in str(feature_id).split("|") if part.startswith("s__")]
    if len(species) != 1:
        return None
    return canonical_taxon(species[0])


def _source_shared_microbiome_features(path: Path) -> list[str]:
    registry = pd.read_csv(path)
    required = registry.loc[
        registry["modality"].eq("microbiome"), "canonical_feature_id"
    ].dropna().astype(str).tolist()
    if not required or len(required) != len(set(required)):
        raise ValueError("Source shared-feature registry must contain unique microbiome features")
    return required


def prepare_locked_external_target(
    relative_abundance_csv: Path,
    sample_metadata_csv: Path,
    source_shared_features_csv: Path,
    output_dir: Path,
) -> dict[str, object]:
    """Create a label-audited Nielsen target using the already frozen source space.

    The 126 species are loaded from a contract created before this target cohort
    is observed.  There is no target-side prevalence filter, feature ranking,
    tuning, or removal of unusable covariates in this function.
    """
    raw = pd.read_csv(relative_abundance_csv, index_col=0, low_memory=False)
    numeric = raw.apply(pd.to_numeric, errors="coerce")
    invalid = raw.notna() & numeric.isna()
    if invalid.any().any():
        bad_columns = invalid.columns[invalid.any(axis=0)].astype(str).tolist()
        raise ValueError(
            "Nielsen abundance matrix contains non-numeric values in columns: "
            + ", ".join(bad_columns[:10])
        )
    if np.isinf(numeric.to_numpy(dtype=float)).any():
        raise ValueError("Nielsen abundance matrix contains infinite values")
    if (numeric < 0).any().any():
        raise ValueError("Nielsen abundance matrix contains negative abundances")
    raw = numeric.fillna(0.0)
    raw.index = raw.index.astype(str)
    raw.columns = raw.columns.astype(str)
    if raw.index.duplicated().any() or raw.columns.duplicated().any():
        raise ValueError("Exported abundance matrix has duplicate feature or sample IDs")
    required_features = _source_shared_microbiome_features(source_shared_features_csv)
    canonical = raw.index.to_series().map(terminal_species_label)
    species_matrix = raw.loc[canonical.notna()].copy()
    species_matrix.index = canonical.loc[canonical.notna()].tolist()
    species_matrix = species_matrix.groupby(level=0, sort=True).sum()
    missing_features = sorted(set(required_features).difference(species_matrix.index))
    if missing_features:
        raise ValueError(f"NielsenHB_2014 lacks {len(missing_features)} fixed source features")
    projected = species_matrix.loc[required_features].copy()

    metadata = pd.read_csv(sample_metadata_csv, low_memory=False)
    required_metadata = {"sample_id", "subject_id", "study_condition", "disease_subtype"}
    absent_metadata = sorted(required_metadata.difference(metadata.columns))
    if absent_metadata:
        raise ValueError(f"Nielsen metadata lacks {absent_metadata}")
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    metadata = metadata.set_index("sample_id", drop=False)
    missing_metadata = sorted(set(projected.columns).difference(metadata.index))
    if missing_metadata:
        raise ValueError(f"Nielsen metadata lacks {len(missing_metadata)} matrix samples")
    manifest = metadata.loc[projected.columns].copy()
    conditions = set(manifest["study_condition"].dropna().astype(str))
    if conditions != {"control", "IBD"}:
        raise ValueError(f"Unexpected study_condition values: {sorted(conditions)}")
    manifest["cohort"] = STUDY_NAME
    manifest["participant_id"] = manifest["subject_id"].astype(str)
    manifest["outcome_ibd"] = manifest["study_condition"].eq("IBD").astype(int)
    manifest["diagnosis_subtype"] = manifest["disease_subtype"].fillna("Control")
    manifest["country"] = manifest.get("country", pd.Series("unknown", index=manifest.index))
    manifest["origin"] = manifest["country"].fillna("unknown").astype(str)
    manifest["site"] = "unknown"
    manifest["timepoint"] = manifest.get(
        "days_from_first_collection", pd.Series(pd.NA, index=manifest.index)
    )
    manifest["age_years"] = manifest.get("age", pd.Series(pd.NA, index=manifest.index))
    manifest["sex"] = manifest.get("gender", pd.Series(pd.NA, index=manifest.index))
    manifest["treatment_observed"] = manifest.get(
        "treatment", pd.Series(pd.NA, index=manifest.index)
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    projected.to_csv(output_dir / "nielsenhb_2014_source_feature_projection.csv", index=True)
    manifest.reset_index(drop=True).to_csv(
        output_dir / "cohort_manifest_nielsenhb_2014.csv", index=False
    )
    covariates = [
        "age_years",
        "sex",
        "treatment_observed",
        "calprotectin",
        "HBI",
        "SCCAI",
        "country",
        "timepoint",
    ]
    missingness = {
        name: int(manifest[name].isna().sum()) if name in manifest else int(len(manifest))
        for name in covariates
    }
    summary = {
        "cohort": STUDY_NAME,
        "role": "locked_external_target_candidate",
        "n_samples": int(manifest.shape[0]),
        "n_participants": int(manifest["participant_id"].nunique()),
        "n_repeated_samples": int(manifest["participant_id"].duplicated().sum()),
        "outcome_counts": manifest["outcome_ibd"].value_counts().sort_index().to_dict(),
        "subtype_counts": manifest["diagnosis_subtype"].value_counts().to_dict(),
        "n_fixed_source_features": len(required_features),
        "n_terminal_species_source_rows": int(canonical.notna().sum()),
        "n_missing_fixed_source_features": len(missing_features),
        "covariate_missingness": missingness,
        "feature_selection": "source_contract_only",
        "target_labels_used_for_feature_selection": False,
    }
    (output_dir / "target_preparation_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary
