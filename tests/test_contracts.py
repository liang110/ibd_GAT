import pandas as pd
import pytest

from ibd_microbiome_framework.data.contracts import MatrixSpec, build_manifest
from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation


def _spec(tmp_path, cohort, group="Control"):
    metadata = tmp_path / f"{cohort.lower()}_metadata.csv"
    pd.DataFrame(
        {"sample_id": ["s1"], "group": [group], "Participant ID": ["p1"]}
    ).to_csv(metadata, index=False)
    return MatrixSpec(
        cohort=cohort,
        modality="microbiome",
        path=tmp_path / "matrix.csv",
        feature_column="feature",
        metadata_path=metadata,
    )


@pytest.mark.parametrize(
    "cohort,group,expected",
    [("HMP2", "Control", 0), ("HMP2", "CD", 1), ("PRISM", "UC", 1),
     ("NielsenHB_2014", "healthy", 0), ("NielsenHB_2014", "IBD", 1)],
)
def test_explicit_label_allow_list(tmp_path, cohort, group, expected):
    manifest = build_manifest(_spec(tmp_path, cohort, group), ["s1"])
    assert manifest.loc[0, "outcome_ibd"] == expected


@pytest.mark.parametrize("bad_label", ["", "unknown", "case"])
def test_unknown_or_missing_labels_fail_closed(tmp_path, bad_label):
    with pytest.raises(ValueError, match="label mapping is not closed"):
        build_manifest(_spec(tmp_path, "HMP2", bad_label), ["s1"])


def test_missing_labels_fail_closed(tmp_path):
    spec = _spec(tmp_path, "HMP2")
    frame = pd.read_csv(spec.metadata_path)
    frame.loc[0, "group"] = pd.NA
    frame.to_csv(spec.metadata_path, index=False)
    with pytest.raises(ValueError, match="label mapping is not closed"):
        build_manifest(spec, ["s1"])


def _longitudinal_spec(tmp_path, **columns):
    rows = 60
    frame = pd.DataFrame({"sample_id": [f"s{i}" for i in range(rows)],
                          "group": ["CD" if i // 3 % 2 else "Control" for i in range(rows)],
                          **columns})
    path = tmp_path / "hmp2.csv"
    frame.to_csv(path, index=False)
    return MatrixSpec("HMP2", "microbiome", tmp_path / "matrix.csv", "feature", path), frame


@pytest.mark.parametrize("alias", ["Participant ID", "participant_id"])
def test_hmp2_alias_groups_all_records_by_person(tmp_path, alias):
    ids = [f"p{i // 3}" for i in range(60)]
    spec, frame = _longitudinal_spec(tmp_path, **{alias: ids})
    manifest = build_manifest(spec, frame.sample_id)
    assert manifest.participant_id.nunique() == 20
    for _, train, validation in repeated_group_train_validation(manifest, 10, .2, 20261110):
        assert set(manifest.iloc[train].participant_id).isdisjoint(manifest.iloc[validation].participant_id)


def test_hmp2_both_aliases_are_normalized_and_coalesced(tmp_path):
    ids = [f"p{i // 3}" for i in range(60)]
    a = [f" {v} " for v in ids]
    a[0] = None
    spec, frame = _longitudinal_spec(tmp_path, **{"Participant ID": a, "participant_id": ids})
    assert set(build_manifest(spec, frame.sample_id).participant_id) == set(ids)


def test_hmp2_conflicting_aliases_fail(tmp_path):
    ids = [f"p{i // 3}" for i in range(60)]
    wrong = ids.copy()
    wrong[0] = "p999"
    spec, frame = _longitudinal_spec(tmp_path, **{"Participant ID": ids, "participant_id": wrong})
    with pytest.raises(ValueError, match="identifier columns conflict"):
        build_manifest(spec, frame.sample_id)


@pytest.mark.parametrize("values", [None, ["" if i == 0 else f"p{i // 3}" for i in range(60)]])
def test_hmp2_missing_identifier_fails(tmp_path, values):
    spec, frame = _longitudinal_spec(tmp_path, **({} if values is None else {"participant_id": values}))
    with pytest.raises(ValueError, match="participant identifier"):
        build_manifest(spec, frame.sample_id)


def test_prism_retains_cross_sectional_sample_identifier(tmp_path):
    spec = _spec(tmp_path, "PRISM")
    assert build_manifest(spec, ["s1"]).loc[0, "participant_id"] == "s1"
