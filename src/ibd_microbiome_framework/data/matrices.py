"""Load and align microbial abundance matrices in the shared taxon space."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from ibd_microbiome_framework.harmonization.features import canonical_taxon, is_forbidden_feature


def load_canonical_microbiome(path: Path, manifest: pd.DataFrame, shared_features: Sequence[str]) -> pd.DataFrame:
    """Return samples x shared taxa after deterministic name mapping and closure."""
    raw = pd.read_csv(path, index_col=0, low_memory=False)
    # PRISM profiles include clinical metadata among the taxon rows.
    keep = [not is_forbidden_feature(name) for name in raw.index]
    raw = raw.loc[keep]
    numeric = raw.apply(pd.to_numeric, errors="coerce")
    invalid_cells = raw.notna() & numeric.isna()
    if invalid_cells.to_numpy().any():
        raise ValueError(
            f"{path.name} contains {int(invalid_cells.to_numpy().sum())} non-numeric abundance cells; "
            "review the input rather than coercing them to zero"
        )
    # Empty cells denote zero abundance.
    numeric = numeric.fillna(0.0)
    keep = [not is_forbidden_feature(name) for name in numeric.index]
    numeric = numeric.loc[keep]
    numeric.index = [canonical_taxon(name) for name in numeric.index]
    numeric = numeric.loc[numeric.index.notna()]
    numeric = numeric.groupby(level=0, sort=True).sum()
    requested = list(shared_features)
    absent = sorted(set(requested).difference(numeric.index))
    if absent:
        raise ValueError(f"{path.name} lacks {len(absent)} required shared taxa")
    sample_ids = manifest["sample_id"].tolist()
    absent_samples = sorted(set(sample_ids).difference(numeric.columns))
    if absent_samples:
        raise ValueError(f"{path.name} lacks {len(absent_samples)} manifest samples")
    matrix = numeric.loc[requested, sample_ids].T.astype(float)
    totals = matrix.sum(axis=1)
    # Closure is undefined for samples with zero mass on the shared taxa.
    matrix = matrix.loc[totals > 0].copy()
    if matrix.empty:
        raise ValueError(f"{path.name} has no samples with positive shared-space mass")
    return matrix.div(matrix.sum(axis=1), axis=0)


def load_shared_feature_list(contract_dir: Path) -> list[str]:
    shared = pd.read_csv(contract_dir / "shared_feature_candidates.csv")
    features = shared.loc[shared["modality"] == "microbiome", "canonical_feature_id"].tolist()
    if not features:
        raise ValueError("No shared microbiome features in contract")
    return features
