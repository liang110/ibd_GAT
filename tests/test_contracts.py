import pandas as pd
import pytest

from ibd_microbiome_framework.data.contracts import MatrixSpec, build_manifest


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
