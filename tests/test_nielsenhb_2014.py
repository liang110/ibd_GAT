from __future__ import annotations

import json

import pandas as pd

from ibd_microbiome_framework.data.nielsenhb_2014 import prepare_locked_external_target


def test_prepare_locked_external_target_uses_source_features_only(tmp_path):
    matrix = pd.DataFrame(
        {
            "S1": [40.0, 60.0, 10.0],
            "S2": [20.0, 80.0, 30.0],
        },
        index=[
            "k__Bacteria|g__Bacteroides|s__Bacteroides_vulgatus",
            "k__Bacteria|g__Faecalibacterium|s__Faecalibacterium_prausnitzii",
            "k__Bacteria|g__Bacteroides",
        ],
    )
    matrix.index.name = "feature_id"
    matrix_path = tmp_path / "matrix.csv"
    matrix.to_csv(matrix_path)
    metadata = pd.DataFrame(
        {
            "sample_id": ["S1", "S2"],
            "subject_id": ["P1", "P1"],
            "study_condition": ["control", "IBD"],
            "disease_subtype": [pd.NA, "UC"],
            "country": ["DNK", "DNK"],
            "age": [31, 31],
            "gender": ["female", "female"],
            "days_from_first_collection": [0, 90],
        }
    )
    metadata_path = tmp_path / "metadata.csv"
    metadata.to_csv(metadata_path, index=False)
    registry = pd.DataFrame(
        {
            "modality": ["microbiome", "microbiome"],
            "canonical_feature_id": ["bacteroides vulgatus", "faecalibacterium prausnitzii"],
        }
    )
    registry_path = tmp_path / "shared.csv"
    registry.to_csv(registry_path, index=False)

    summary = prepare_locked_external_target(matrix_path, metadata_path, registry_path, tmp_path / "out")

    projected = pd.read_csv(tmp_path / "out" / "nielsenhb_2014_source_feature_projection.csv", index_col=0)
    manifest = pd.read_csv(tmp_path / "out" / "cohort_manifest_nielsenhb_2014.csv")
    audit = json.loads((tmp_path / "out" / "target_preparation_audit.json").read_text())
    assert projected.index.tolist() == registry["canonical_feature_id"].tolist()
    assert projected.shape == (2, 2)
    assert manifest["outcome_ibd"].tolist() == [0, 1]
    assert manifest["diagnosis_subtype"].tolist() == ["Control", "UC"]
    assert summary["n_repeated_samples"] == 1
    assert audit["target_labels_used_for_feature_selection"] is False
