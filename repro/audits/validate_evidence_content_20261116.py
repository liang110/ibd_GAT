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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()
    check_registry(root)
    check_direct_transfer(root)
    check_units(root)
    check_figure3(root)
    print("PASS: registry, direct transfer, unit semantics, common draws and Figure 3 content checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
