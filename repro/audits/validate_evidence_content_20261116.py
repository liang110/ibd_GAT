#!/usr/bin/env python3
"""Content-level checks for the compact manuscript evidence package."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def fail(message: str) -> None:
    raise SystemExit(f"FAIL: {message}")


def source(root: Path, relative: str) -> Path:
    path = root / "results" / "reconciled_evidence" / relative
    if not path.exists():
        fail(f"missing evidence source: {relative}")
    return path


def check_registry(root: Path) -> pd.DataFrame:
    path = root / "results/reconciled_evidence/evidence_registry.csv"
    if not path.exists():
        fail("missing evidence_registry.csv")
    registry = pd.read_csv(path)
    required = {"object", "primary_sources", "status", "limitation"}
    missing = required.difference(registry.columns)
    if missing:
        fail(f"registry missing columns: {sorted(missing)}")
    row = registry.loc[registry["object"].eq("Supplementary Figure S4")]
    if len(row) != 1:
        fail("registry must contain exactly one Supplementary Figure S4 row")
    s4_sources = str(row.iloc[0]["primary_sources"])
    if "prism_discovery_validation_transfer/" not in s4_sources:
        fail("Supplementary Figure S4 is not registered to direct stratum transfer")
    for item in s4_sources.split(";"):
        source(root, item.strip())
    return registry


def check_direct_transfer(root: Path) -> None:
    directory = root / "results/reconciled_evidence/prism_discovery_validation_transfer"
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    expected_directions = {
        "PRISM_discovery_to_validation",
        "PRISM_validation_to_discovery",
    }
    if set(manifest.get("directions", [])) != expected_directions:
        fail("direct-transfer manifest directions do not match the declared experiment")
    if manifest.get("n_splits") != 10 or manifest.get("n_bootstrap") != 5000:
        fail("direct-transfer manifest must declare 10 source splits and 5,000 draws")
    if set(manifest.get("models", [])) != {
        "graph_attention",
        "elastic_net",
        "random_forest",
        "no_message_mlp",
    }:
        fail("direct-transfer manifest model set is incomplete or mismatched")

    metrics = pd.read_csv(directory / "metrics.csv")
    if len(metrics) != 80:
        fail(f"direct-transfer metrics must contain 80 rows, found {len(metrics)}")
    if set(metrics["direction"]) != expected_directions:
        fail("direct-transfer metrics contain an unexpected direction")
    if set(metrics["model"]) != {
        "graph_attention",
        "elastic_net",
        "random_forest",
        "no_message_mlp",
    }:
        fail("direct-transfer metrics model set is incomplete")
    if metrics.groupby("direction")["split"].nunique().to_dict() != {
        "PRISM_discovery_to_validation": 10,
        "PRISM_validation_to_discovery": 10,
    }:
        fail("direct-transfer split counts are not ten per direction")

    summary = pd.read_csv(directory / "summary.csv")
    if len(summary) != 8 or set(summary["direction"]) != expected_directions:
        fail("direct-transfer summary must contain eight direction/model rows")
    contrasts = pd.read_csv(directory / "paired_bootstrap_summary_b5000.csv")
    if len(contrasts) != 6:
        fail("direct-transfer paired summary must contain six contrasts")
    if set(contrasts["direction"]) != expected_directions:
        fail("direct-transfer paired summary direction mismatch")
    if set(contrasts["n_bootstrap"]) != {5000}:
        fail("direct-transfer paired summary is not a 5,000-draw result")


def check_units(root: Path) -> None:
    path = root / "results/reconciled_evidence/paired_bootstrap_reconciled_by_split.csv"
    df = pd.read_csv(path)
    required = {"scoring_unit", "bootstrap_unit", "aggregation_unit", "draw_scope"}
    missing = required.difference(df.columns)
    if missing:
        fail(f"paired bootstrap table missing explicit unit columns: {sorted(missing)}")
    if set(df["scoring_unit"]) != {"visit"} or set(df["bootstrap_unit"]) != {"participant"}:
        fail("paired bootstrap unit convention is not visit-scored/participant-clustered")
    if set(df["aggregation_unit"]) != {"source_split"}:
        fail("paired bootstrap aggregation unit is not source_split")
    if set(df["n_bootstrap"]) != {5000}:
        fail("paired bootstrap table must declare 5,000 draws")

    common = pd.read_csv(root / "results/reconciled_evidence/common_draw/gate_a_common_draw_summary.csv")
    if set(common["n_bootstrap"]) != {5000}:
        fail("common-draw summary must declare 5,000 draws")
    if set(common["analysis_level"]) != {"visit", "participant"}:
        fail("common-draw summary must expose visit and participant scoring levels")


def check_figure3(root: Path) -> None:
    path = root / "manuscript/source_data/Figure3_repeated_leakage.csv"
    df = pd.read_csv(path)
    if len(df) != 10:
        fail("Figure 3 source table must contain ten seeds")
    if abs(df["random_split_subject_overlap"].mean() - 43.0) > 1e-9:
        fail("Figure 3 mean sample-split participant overlap is not 43")
    if df["group_split_subject_overlap"].max() != 0:
        fail("Figure 3 grouped split contains participant overlap")
    if df["acceptance_pass"].astype(bool).sum() != 10:
        fail("Figure 3 acceptance record is not 10/10")


def check_current_contract_registry(root: Path) -> None:
    table = pd.read_csv(root / "results/submission_tables/Table_1_cohort_contract.csv")
    if set(table["Cohort"]) != {"HMP2/IBDMDB", "PRISM"}:
        fail("Table 1 contains a non-current cohort or is missing HMP2/PRISM")
    hmp2 = table.loc[table["Cohort"].eq("HMP2/IBDMDB")].iloc[0]
    prism = table.loc[table["Cohort"].eq("PRISM")].iloc[0]
    included_column = "Included observations" if "Included observations" in table.columns else "Valid prediction samples"
    if (int(hmp2[included_column]), int(hmp2["Participants"])) != (1626, 130):
        fail("Table 1 HMP2 current-contract counts are incorrect")
    if (int(prism[included_column]), int(prism["Participants"])) != (220, 220):
        fail("Table 1 PRISM current-contract counts are incorrect")

    registry_path = root / "manuscript/registry/figure_table_source_registry.csv"
    registry = pd.read_csv(registry_path)
    required_ids = {"Figure1", "Figure2", "Figure3", "Figure4", "Figure5", "Figure6", "Figure7", "Table1", "Table2", "Table3", "Table4", "Supplementary Figure S4"}
    if not required_ids.issubset(set(registry["figure_id"])):
        fail("public source registry is missing one or more current manuscript objects")
    for row in registry.itertuples():
        for raw in str(row.source_artifact).split(";"):
            rel = raw.strip()
            if rel and not (root / rel).exists():
                fail(f"public source registry points to a missing artifact: {rel}")
    by_id = {row.figure_id: str(row.source_artifact) for row in registry.itertuples()}
    if "common_draw/gate_a_common_draw_summary.csv" not in by_id.get("Table3", ""):
        fail("Table 3 registry does not point to the current common-draw summary")
    if "prism_discovery_validation_transfer" not in by_id.get("Supplementary Figure S4", ""):
        fail("Supplementary Figure S4 registry does not point to direct stratum transfer")
    if "formal_topology_runtime_manifest.json" not in by_id.get("Figure6", ""):
        fail("Figure 6 registry does not include formal topology provenance")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()
    check_registry(root)
    check_direct_transfer(root)
    check_units(root)
    check_figure3(root)
    check_current_contract_registry(root)
    print("PASS: registry, direct transfer, unit semantics, common draws and Figure 3 content checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
