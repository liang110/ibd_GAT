"""Regression tests for clinical measurement names in microbial inputs."""
from pathlib import Path

import pandas as pd
import pytest

from ibd_microbiome_framework.data.contracts import MatrixSpec, build_feature_registry
from ibd_microbiome_framework.data.matrices import load_canonical_microbiome
from ibd_microbiome_framework.harmonization.features import (
    canonical_metabolite, canonical_taxon, is_forbidden_feature,
)


@pytest.mark.parametrize("name", [
    "Fecal.Calprotectin", "Faecal.Calprotectin", "faecal_calprotectin",
    "FECAL-CALPROTECTIN", "Faecal Calprotectin", "Stool Calprotectin",
    "fecal calprotectin concentration", "faecal calprotectin (ug/g)",
    "fecal-calprotectin-level [mg/kg]", "calprotectin",
])
def test_explicit_calprotectin_aliases_are_excluded(name):
    assert is_forbidden_feature(name)
    assert canonical_taxon(name) is None
    assert canonical_metabolite(name) is None


@pytest.mark.parametrize("name", [
    "Bacteroides faecalis", "Streptococcus fecalis",
    "Calprotectinibacter intestinalis", "Clostridium calprotectinum",
])
def test_taxa_are_not_removed_by_substring(name):
    assert not is_forbidden_feature(name)
    assert canonical_taxon(name) is not None


def test_registry_and_matrix_loader_exclude_aliases(tmp_path: Path):
    aliases = ["Fecal.Calprotectin", "Faecal.Calprotectin", "FECAL-CALPROTECTIN"]
    taxa = ["Bacteroides fragilis", "Bacteroides vulgatus"]
    specs = [MatrixSpec(c, "microbiome", tmp_path / (c + ".csv"), "feature", tmp_path / "metadata.csv")
             for c in ("HMP2", "PRISM")]
    registry = build_feature_registry(specs, {s.cohort + ":microbiome": taxa + aliases for s in specs})
    clinical = registry[registry.source_feature_id.isin(aliases)]
    assert clinical.is_forbidden_clinical.all()
    assert clinical.canonical_feature_id.isna().all()
    shared = registry[registry.shared_across_hmp2_prism.eq(True)].canonical_feature_id.drop_duplicates().tolist()
    assert set(shared) == {x.lower() for x in taxa}
    pd.DataFrame({"feature": taxa + aliases, "sample1": [1, 3, "x", "x", "x"]}).to_csv(specs[0].path, index=False)
    matrix = load_canonical_microbiome(specs[0].path, pd.DataFrame({"sample_id": ["sample1"]}), shared)
    assert set(matrix.columns) == set(shared)
    assert matrix.shape == (1, 2)
    assert matrix.sum(axis=1).iloc[0] == pytest.approx(1)
