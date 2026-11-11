"""Build auditable cohort and feature contracts from frozen delimited matrices."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Sequence

import numpy as np
import pandas as pd

from ibd_microbiome_framework.harmonization.features import (
    canonical_metabolite,
    canonical_taxon,
    is_forbidden_feature,
)


# Label mapping is deliberately explicit.  A missing or previously unseen label
# must stop contract construction rather than being silently treated as IBD.
ALLOWED_CONTROL_LABELS = {
    "Control", "control", "nonIBD", "non-IBD", "non-IBD/control",
}
ALLOWED_IBD_LABELS = {
    "IBD", "ibd", "CD", "cd", "UC", "uc", "Crohn's disease",
    "Crohn disease", "Crohn's disease (CD)", "Ulcerative colitis",
    "ulcerative colitis", "Ulcerative colitis (UC)",
}
ALLOWED_NIELSEN_CONTROL_LABELS = {"healthy", "Healthy", "control", "Control"}
ALLOWED_NIELSEN_IBD_LABELS = {"IBD", "ibd"}


@dataclass(frozen=True)
class MatrixSpec:
    cohort: str
    modality: str
    path: Path
    feature_column: str
    metadata_path: Path


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(block_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def matrix_axes(path: Path) -> tuple[List[str], List[str]]:
    """Read header and first-column IDs without loading the numeric matrix."""
    samples: List[str] = []
    features: List[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        samples = [cell.strip() for cell in header[1:] if cell.strip()]
        for row in reader:
            if row and row[0].strip():
                features.append(row[0].strip())
    return samples, features


def abundance_library_sizes(path: Path) -> pd.Series:
    """Compute sample sums for a modest abundance matrix using numeric coercion."""
    frame = pd.read_csv(path, index_col=0, low_memory=False)
    numeric = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return numeric.sum(axis=0)


def infer_prism_origin(sample_id: str) -> str:
    if sample_id.startswith("Validation|"):
        return "validation"
    if sample_id.startswith("PRISM|"):
        return "discovery"
    return "unknown"


def build_manifest(spec: MatrixSpec, matrix_samples: Sequence[str]) -> pd.DataFrame:
    metadata = pd.read_csv(spec.metadata_path)
    if "sample_id" not in metadata.columns:
        raise ValueError(f"{spec.metadata_path} lacks sample_id")
    manifest = metadata.copy()
    manifest = manifest[manifest["sample_id"].isin(matrix_samples)].copy()
    manifest["cohort"] = spec.cohort
    if spec.cohort == "PRISM":
        manifest["origin"] = manifest["sample_id"].map(infer_prism_origin)
        manifest["participant_id"] = manifest["sample_id"]
        manifest["site"] = "unknown"
        manifest["timepoint"] = "cross_sectional"
    elif spec.cohort == "NielsenHB_2014":
        # Nielsen origin is carried by the cohort-specific loader when present.
        # Keep the generic contract path explicit rather than falling back to HMP2.
        manifest["origin"] = manifest.get("country", spec.cohort)
        manifest["participant_id"] = manifest.get("Participant ID", manifest["sample_id"])
        manifest["site"] = manifest.get("site_name", "unknown")
        manifest["timepoint"] = manifest.get("timepoint", np.nan)
    else:
        manifest["origin"] = "HMP2"
        manifest["participant_id"] = manifest.get("Participant ID", manifest["sample_id"])
        manifest["site"] = manifest.get("site_name", "unknown")
        manifest["timepoint"] = manifest.get("week_num", np.nan)
    if "group" not in manifest.columns:
        raise ValueError(f"{spec.metadata_path} lacks the required raw 'group' label column")
    raw = manifest["group"].astype("string").str.strip()
    observed = set(raw.dropna().tolist())
    if spec.cohort == "NielsenHB_2014":
        control_labels = ALLOWED_NIELSEN_CONTROL_LABELS
        ibd_labels = ALLOWED_NIELSEN_IBD_LABELS
    else:
        control_labels = ALLOWED_CONTROL_LABELS
        ibd_labels = ALLOWED_IBD_LABELS
    allowed = control_labels | ibd_labels
    unknown = sorted(observed - allowed)
    if raw.isna().any() or unknown:
        missing_n = int(raw.isna().sum())
        raise ValueError(
            f"{spec.cohort}: label mapping is not closed; missing={missing_n}, "
            f"unknown={unknown}. Add an explicit mapping before building the contract."
        )
    manifest["outcome_ibd"] = raw.map(
        {label: 0 for label in control_labels}
        | {label: 1 for label in ibd_labels}
    ).astype(int)
    manifest["diagnosis_subtype"] = manifest.get("diagnosis", manifest["group"])
    if manifest["sample_id"].duplicated().any():
        raise ValueError(f"Duplicate sample IDs in {spec.cohort} manifest")
    return manifest.sort_values("sample_id").reset_index(drop=True)


def build_feature_registry(
    specs: Sequence[MatrixSpec],
    feature_lists: Mapping[str, Sequence[str]],
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for spec in specs:
        canonicalizer: Callable[[object], object]
        canonicalizer = canonical_taxon if spec.modality == "microbiome" else canonical_metabolite
        for raw_feature in feature_lists[spec.cohort + ":" + spec.modality]:
            canonical = canonicalizer(raw_feature)
            rows.append(
                {
                    "cohort": spec.cohort,
                    "modality": spec.modality,
                    "source_feature_id": raw_feature,
                    "canonical_feature_id": canonical,
                    "is_forbidden_clinical": is_forbidden_feature(raw_feature),
                    "mapping_status": "exact_normalized" if canonical else "excluded_or_unresolved",
                }
            )
    registry = pd.DataFrame(rows)
    valid = registry.dropna(subset=["canonical_feature_id"])
    presence = (
        valid.groupby(["modality", "canonical_feature_id"])["cohort"]
        .agg(lambda values: ";".join(sorted(set(values))))
        .rename("present_in")
        .reset_index()
    )
    presence["shared_across_hmp2_prism"] = presence["present_in"].eq("HMP2;PRISM")
    return registry.merge(presence, on=["modality", "canonical_feature_id"], how="left")


def write_contracts(
    specs: Sequence[MatrixSpec], output_dir: Path, source_root: Path
) -> Dict[str, object]:
    """Write manifests and a feature registry; never alter source matrices."""
    output_dir.mkdir(parents=True, exist_ok=True)
    axes: Dict[str, tuple[List[str], List[str]]] = {}
    feature_lists: Dict[str, Sequence[str]] = {}
    source_inventory: List[Dict[str, object]] = []
    for spec in specs:
        samples, features = matrix_axes(spec.path)
        key = spec.cohort + ":" + spec.modality
        axes[key] = (samples, features)
        feature_lists[key] = features
        source_inventory.append(
            {
                "cohort": spec.cohort,
                "modality": spec.modality,
                "path": str(spec.path),
                "metadata_path": str(spec.metadata_path),
                "sha256": sha256_file(spec.path),
                "n_samples_header": len(samples),
                "n_features_raw": len(features),
            }
        )
        manifest = build_manifest(spec, samples)
        manifest.to_csv(output_dir / f"cohort_manifest_{spec.cohort.lower()}.csv", index=False)

    registry = build_feature_registry(specs, feature_lists)
    registry.to_csv(output_dir / "feature_registry.csv", index=False)
    pd.DataFrame(source_inventory).to_csv(output_dir / "source_inventory.csv", index=False)

    library_rows: List[Dict[str, object]] = []
    for spec in specs:
        if spec.modality != "microbiome":
            continue
        sums = abundance_library_sizes(spec.path)
        for sample_id, value in sums.items():
            library_rows.append(
                {"cohort": spec.cohort, "sample_id": sample_id, "library_sum": float(value)}
            )
    library = pd.DataFrame(library_rows)
    library.to_csv(output_dir / "microbiome_library_sums.csv", index=False)

    shared = registry.loc[
        registry["shared_across_hmp2_prism"].fillna(False)
        & ~registry["is_forbidden_clinical"],
        ["modality", "canonical_feature_id"],
    ].drop_duplicates()
    shared_microbiome = shared.loc[shared["modality"] == "microbiome", "canonical_feature_id"].tolist()
    shared.to_csv(output_dir / "shared_feature_candidates.csv", index=False)
    contract_hash = hashlib.sha256("\n".join(shared_microbiome).encode("utf-8")).hexdigest()
    canonicalizer_hash = sha256_file(
        Path(__file__).resolve().parents[1] / "harmonization" / "features.py"
    )

    prism_microbiome = next(
        spec for spec in specs if spec.cohort == "PRISM" and spec.modality == "microbiome"
    )
    summary = {
        "source_root": str(source_root),
        "n_source_matrices": len(specs),
        "shared_microbiome_candidates": int((shared["modality"] == "microbiome").sum()),
        "shared_microbiome_contract_sha256": contract_hash,
        "canonicalizer_sha256": canonicalizer_hash,
        "shared_metabolite_exact_name_candidates": int((shared["modality"] == "metabolome").sum()),
        "prism_manifest_origins": build_manifest(prism_microbiome, axes["PRISM:microbiome"][0])["origin"].value_counts().to_dict(),
        "forbidden_clinical_rows": int(registry["is_forbidden_clinical"].sum()),
    }
    (output_dir / "contract_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary
