#!/usr/bin/env python3
"""Verify score summaries, model comparisons and target denominators."""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path


SCHEMAS = {
    "known_truth_acceptance_summary.csv": (10, {"experiment", "scenario", "disease_truth_defined", "n_seeds", "pass_count", "mean_auroc", "mean_truth_recovery"}),
    "primary_common_draw_summary.csv": (72, {"direction", "target_stratum", "analysis_level", "candidate", "reference", "metric", "split_mean_point_delta", "common_draw_bootstrap_mean", "ci_low", "ci_high", "n_source_splits", "n_bootstrap", "n_target_samples", "n_target_participants"}),
    "prism_direct_transfer_metrics.csv": (8, {"direction", "source_origin", "target_origin", "model", "n_splits", "auroc_mean", "auroc_sd", "ap_mean", "ap_sd", "brier_mean", "brier_sd", "target_n"}),
    "prism_direct_transfer_paired.csv": (6, {"direction", "reference", "point_mean_delta_auroc", "bootstrap_mean_delta_auroc", "ci_low", "ci_high", "n_bootstrap", "n_source_splits", "bootstrap_unit", "scoring_unit"}),
    "topology_intervention_paired.csv": (48, {"direction", "target_stratum", "scoring_unit", "reference", "metric", "n_draws", "mean", "ci_low", "ci_high"}),
    "topology_taxon_aligned_summary.csv": (2, {"direction", "topology", "n_pairs", "taxon_mean", "local_mean"}),
    "xgboost_metrics_by_scoring_unit.csv": (12, {"direction", "target_stratum", "scoring_unit", "model", "bootstrap_unit", "n_splits", "n_scored", "roc_auc_mean", "roc_auc_sd", "average_precision_mean", "average_precision_sd", "brier_mean", "brier_sd"}),
    "xgboost_paired_common_draw.csv": (18, {"direction", "target_stratum", "analysis_level", "candidate", "reference", "metric", "split_mean_point_delta", "common_draw_bootstrap_mean", "ci_low", "ci_high", "n_source_splits", "n_bootstrap", "n_target_samples", "n_target_participants"}),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def near(actual, expected, label, tol=0.00005):
    actual = float(actual)
    require(math.isfinite(actual) and abs(actual - expected) < tol,
            f"{label}: got {actual}, expected approximately {expected}")


def one(rows, **key):
    selected = [r for r in rows if all(r.get(k) == v for k, v in key.items())]
    require(len(selected) == 1, f"Expected one row for {key}; got {len(selected)}")
    return selected[0]


def validate_comparisons(directory):
    data = {}
    for name, (count, required) in SCHEMAS.items():
        with (directory / name).open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            require(reader.fieldnames is not None, f"{name}: missing header")
            require(len(reader.fieldnames) == len(set(reader.fieldnames)), f"{name}: duplicate columns")
            missing = required - set(reader.fieldnames)
            require(not missing, f"{name}: missing columns {sorted(missing)}")
            rows = list(reader)
        require(len(rows) == count, f"{name}: expected {count} rows; got {len(rows)}")
        data[name] = rows
    print("PASS: eight CSV schemas and row counts")

    primary = data["primary_common_draw_summary.csv"]
    for direction, expected in (("HMP2_to_PRISM", (0.099078, 0.047623, 0.149063)), ("PRISM_to_HMP2", (-0.083925, -0.150056, -0.018179))):
        r = one(primary, direction=direction, target_stratum="pooled", analysis_level="visit", reference="random_forest", candidate="no_prior_gat", metric="roc_auc")
        for field, val in zip(("common_draw_bootstrap_mean", "ci_low", "ci_high"), expected): near(r[field], val, f"primary {direction} {field}")
        require(r["n_bootstrap"] == "5000" and r["n_source_splits"] == "10", "primary draw/split count")
        require(int(r["n_target_participants"]) == (220 if direction == "HMP2_to_PRISM" else 130), "primary participant denominator")
    print("PASS: bidirectional GAT–RF summary and conditional interval counts")

    xmetrics = data["xgboost_metrics_by_scoring_unit.csv"]
    for direction, unit, expected_n, vals in (("HMP2_to_PRISM", "visit", 220, (0.698955, 0.860016, 0.178552)), ("PRISM_to_HMP2", "visit", 1626, (0.702835, 0.840825, 0.241456)), ("PRISM_to_HMP2", "participant", 130, (0.721143, 0.880780, 0.187841))):
        r = one(xmetrics, direction=direction, target_stratum="pooled", scoring_unit=unit, model="xgboost")
        require(r["bootstrap_unit"] == "participant" and int(r["n_scored"]) == expected_n and r["n_splits"] == "10", "S7 units or denominator")
        for field, value in zip(("roc_auc_mean", "average_precision_mean", "brier_mean"), vals): near(r[field], value, f"S7 {direction}/{unit} {field}")
    xp = data["xgboost_paired_common_draw.csv"]
    for direction, unit, expected in (("HMP2_to_PRISM", "visit", (0.157589, 0.157263, 0.083774, 0.231026)), ("PRISM_to_HMP2", "visit", (-0.089063, -0.088636, -0.161680, -0.016870)), ("PRISM_to_HMP2", "participant", (-0.108738, -0.109168, -0.188918, -0.027451))):
        r = one(xp, direction=direction, target_stratum="pooled", analysis_level=unit, reference="xgboost", candidate="no_prior_gat", metric="roc_auc")
        require(r["n_bootstrap"] == "5000" and r["n_source_splits"] == "10", "XGBoost draws/splits")
        for field, value in zip(("split_mean_point_delta", "common_draw_bootstrap_mean", "ci_low", "ci_high"), expected): near(r[field], value, f"XGBoost {direction}/{unit} {field}")
        require(abs(float(r["split_mean_point_delta"]) - float(r["common_draw_bootstrap_mean"])) > 0.0001, "point/boot mean unexpectedly identical")
        g = one(xmetrics, direction=direction, target_stratum="pooled", scoring_unit=unit, model="no_prior_gat")
        x = one(xmetrics, direction=direction, target_stratum="pooled", scoring_unit=unit, model="xgboost")
        near(float(g["roc_auc_mean"]) - float(x["roc_auc_mean"]), float(r["split_mean_point_delta"]), "absolute split mean vs observed paired difference", 1e-10)
    for stratum, n, absolute, paired in (("discovery", 155, .696330, (.174064, .173010, .088993, .258724)), ("validation", 65, .716913, (.114905, .115424, -.027806, .265528))):
        x = one(xmetrics, direction="HMP2_to_PRISM", target_stratum=stratum, scoring_unit="visit", model="xgboost")
        g = one(xmetrics, direction="HMP2_to_PRISM", target_stratum=stratum, scoring_unit="visit", model="no_prior_gat")
        require(int(x["n_scored"]) == n and int(g["n_scored"]) == n and x["bootstrap_unit"] == "participant", "origin scoring unit/denominator")
        near(x["roc_auc_mean"], absolute, "origin XGBoost absolute AUROC")
        r = one(xp, direction="HMP2_to_PRISM", target_stratum=stratum, analysis_level="visit", candidate="no_prior_gat", reference="xgboost", metric="roc_auc")
        require(int(r["n_target_participants"]) == n and r["n_bootstrap"] == "5000", "origin target participants/draw count")
        for field, value in zip(("split_mean_point_delta", "common_draw_bootstrap_mean", "ci_low", "ci_high"), paired): near(r[field], value, "origin " + field)
        near(float(g["roc_auc_mean"]) - float(x["roc_auc_mean"]), float(r["split_mean_point_delta"]), "origin absolute means vs observed difference", 1e-10)
    print("PASS: S7 visit/participant metrics, target-origin values and observed-vs-bootstrap point distinction")

    direct = data["prism_direct_transfer_metrics.csv"]
    for direction, target_n in (("PRISM_discovery_to_validation", "65"), ("PRISM_validation_to_discovery", "155")):
        for model in ("graph_attention", "elastic_net", "random_forest", "no_message_mlp"):
            r = one(direct, direction=direction, model=model)
            require(r["n_splits"] == "10" and r["target_n"] == target_n, "S3 target denominator/split count")
    f = one(direct, direction="PRISM_discovery_to_validation", model="graph_attention")
    v = one(direct, direction="PRISM_validation_to_discovery", model="graph_attention")
    near(f["auroc_mean"], .552907, "S3 forward GAT AUROC"); near(f["ap_mean"], .714927, "S3 forward GAT AP")
    near(v["auroc_mean"], .697679, "S3 reverse GAT AUROC"); near(v["auroc_sd"], .243583, "S3 reverse GAT SD")
    for r in data["prism_direct_transfer_paired.csv"]:
        require(r["scoring_unit"] == r["bootstrap_unit"] == "participant" and r["n_bootstrap"] == "5000" and r["n_source_splits"] == "10", "direct transfer unit/draw count")
        g = one(direct, direction=r["direction"], model="graph_attention")
        ref = one(direct, direction=r["direction"], model=r["reference"])
        near(float(g["auroc_mean"]) - float(ref["auroc_mean"]), float(r["point_mean_delta_auroc"]), "direct transfer paired point difference", 1e-10)
        require(float(r["ci_low"]) < float(r["ci_high"]), "direct transfer reversed interval")
    print("PASS: Table S3 directions, model set, AUROC, SD and paired observed differences")

    topo = data["topology_taxon_aligned_summary.csv"]
    for direction, mean in (("HMP2_to_PRISM", .522595), ("PRISM_to_HMP2", .530638)):
        r = one(topo, direction=direction, topology="spearman")
        require(r["n_pairs"] == "45", "topology pair count")
        near(r["taxon_mean"], mean, "taxon-aligned Jaccard")
    for direction, expected in (("HMP2_to_PRISM", (.097872, .070440, .125322)), ("PRISM_to_HMP2", (.097146, -.003346, .191511))):
        r = one(data["topology_intervention_paired.csv"], direction=direction, target_stratum="pooled", scoring_unit="visit", reference="rewired", metric="roc_auc")
        require(r["n_draws"] == "5000", "topology draw count")
        for field, val in zip(("mean", "ci_low", "ci_high"), expected): near(r[field], val, "topology " + field)
    print("PASS: taxon-aligned edge summary and Real–Rewired conditional intervals")

    truth = data["known_truth_acceptance_summary.csv"]
    for scenario, count in (("null_control", 10), ("batch_only", 10), ("signal_only", 8), ("signal_reversing_batch", 7), ("signal_topology_shift", 10)):
        r = one(truth, experiment="end_to_end_truth_panel", scenario=scenario)
        require(int(r["n_seeds"]) == 10 and int(r["pass_count"]) == count, "truth acceptance count")
        if scenario in ("null_control", "batch_only"):
            require(r["disease_truth_defined"].lower() == "false" and not r["mean_truth_recovery"].strip(), "negative-control recall must be NA")
        else: require(r["disease_truth_defined"].lower() == "true", "positive truth missing")
    print("PASS: end-to-end truth acceptance counts and negative-control NA")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    validate_comparisons(parser.parse_args().root / "results/comparisons")
