from pathlib import Path
import pandas as pd
import pytest
from ibd_microbiome_framework.data.contracts import MatrixSpec, write_contracts


def make_specs(tmp_path):
    specs = []
    for cohort in ("HMP2", "PRISM"):
        ids = [f"{cohort}|s{i}" for i in range(3)]
        metadata = tmp_path / f"{cohort}_metadata.csv"
        pd.DataFrame({"sample_id": ids, "group": ["CD", "Control", "CD"]}).to_csv(metadata, index=False)
        for modality, samples in (("microbiome", ids), ("metabolome", ids[:1])):
            matrix = tmp_path / f"{cohort}_{modality}.csv"
            pd.DataFrame([[1.0] * len(samples)], index=["Bacteroides vulgatus" if modality == "microbiome" else "citrate"], columns=samples).to_csv(matrix, index_label="feature")
            specs.append(MatrixSpec(cohort, modality, matrix, "feature", metadata))
    return specs


@pytest.mark.parametrize("reverse", [False, True])
def test_mixed_modality_manifests_are_separate_and_alias_is_microbiome(tmp_path, reverse):
    specs = make_specs(tmp_path)
    if reverse:
        specs = list(reversed(specs))
    output = tmp_path / "out"
    write_contracts(specs, output, tmp_path)
    for cohort in ("hmp2", "prism"):
        micro = output / f"cohort_manifest_{cohort}_microbiome.csv"
        metabolome = output / f"cohort_manifest_{cohort}_metabolome.csv"
        alias = output / f"cohort_manifest_{cohort}.csv"
        assert len(pd.read_csv(micro)) == 3
        assert len(pd.read_csv(metabolome)) == 1
        assert alias.read_bytes() == micro.read_bytes()
        assert set(pd.read_csv(micro).sample_id) != set(pd.read_csv(metabolome).sample_id)


def test_order_swap_produces_same_modality_manifests(tmp_path):
    specs = make_specs(tmp_path)
    a, b = tmp_path / "forward", tmp_path / "backward"
    write_contracts(specs, a, tmp_path)
    write_contracts(list(reversed(specs)), b, tmp_path)
    for cohort in ("hmp2", "prism"):
        for suffix in ("microbiome", "metabolome"):
            name = f"cohort_manifest_{cohort}_{suffix}.csv"
            assert (a / name).read_bytes() == (b / name).read_bytes()
        assert (a / f"cohort_manifest_{cohort}.csv").read_bytes() == (b / f"cohort_manifest_{cohort}.csv").read_bytes()


def test_duplicate_cohort_modality_fails_before_output(tmp_path):
    specs = make_specs(tmp_path)
    output = tmp_path / "out"
    with pytest.raises(ValueError, match="Duplicate cohort/modality"):
        write_contracts([*specs, specs[0]], output, tmp_path)
    assert not output.exists()


def test_duplicate_pair_ignores_case(tmp_path):
    specs = make_specs(tmp_path)
    first = specs[0]
    duplicate = MatrixSpec(first.cohort.lower(), first.modality.upper(), first.path, first.feature_column, first.metadata_path)
    with pytest.raises(ValueError, match="Duplicate cohort/modality"):
        write_contracts([*specs, duplicate], tmp_path / "out", tmp_path)
