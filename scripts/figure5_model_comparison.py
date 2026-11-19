"""Assemble Supplementary Table 12 and Supplementary Figure 3 from the released observed and bootstrap summaries."""
import argparse
import csv
from pathlib import Path


SETTINGS = [
    ("HMP2", "PRISM discovery", "record", "HMP2_to_PRISM", "discovery", "visit", 155, 155),
    ("HMP2", "Netherlands validation", "record", "HMP2_to_PRISM", "validation", "visit", 65, 65),
    ("Franzosa", "HMP2", "record", "PRISM_to_HMP2", "pooled", "visit", 1626, 130),
    ("Franzosa", "HMP2", "participant", "PRISM_to_HMP2", "pooled", "participant", 1626, 130),
]
MODELS = {
    "Graph attention": "no_prior_gat",
    "Elastic net": "elastic_net",
    "Random forest": "random_forest",
    "No-message MLP": "capacity_matched_no_message_passing",
    "XGBoost": "xgboost",
}


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def unique(rows, **keys):
    matches = [row for row in rows if all(row[k] == v for k, v in keys.items())]
    if len(matches) != 1:
        raise ValueError(f"Expected one source row for {keys}, found {len(matches)}")
    return matches[0]


def assemble(root):
    results = root / "results"
    final = results / "p0_rechecked_20261117"
    main = read(results / "reconciled_evidence/absolute_performance_summary.csv")
    boosted = read(final / "xgboost_metrics_by_scoring_unit.csv")
    contrasts = read(final / "primary_common_draw_summary.csv") + read(final / "xgboost_paired_common_draw.csv")
    output = []
    for source, target, unit, direction, stratum, level, records, people in SETTINGS:
        for model, code in MODELS.items():
            is_xgb = model == "XGBoost"
            if is_xgb:
                absolute = unique(boosted, direction=direction, target_stratum=stratum, scoring_unit=level, model=code)
                if int(absolute["n_scored"]) != (records if unit == "record" else people):
                    raise ValueError("XGBoost target denominator mismatch")
            else:
                origin = "HMP2" if target == "HMP2" else stratum
                absolute = unique(main, direction=direction, target_origin=origin, analysis_level=level, model=model)
                if float(absolute["n_participants"]) != people:
                    raise ValueError("Primary target denominator mismatch")
            if int(absolute["n_splits" if is_xgb else "n_source_splits"]) != 10:
                raise ValueError(f"Unexpected source-fit count: {source}, {target}, {model}")
            row = dict(source=source, target=target, scoring_unit=unit, model=model,
                       records=records, participants=people, n_fits=10,
                       training_scope="complete-source refit" if is_xgb else "grouped source training partition")
            for name in ("auroc", "average_precision", "brier"):
                key = "roc_auc" if is_xgb and name == "auroc" else name
                for statistic in ("mean", "sd"):
                    row[f"{name}_{statistic}"] = absolute[f"{key}_{statistic}"]
            paired = {}
            if model != "Graph attention":
                paired = unique(contrasts, direction=direction, target_stratum=stratum, analysis_level=level,
                                candidate="no_prior_gat", reference=code, metric="roc_auc")
                if int(paired["n_bootstrap"]) != 5000 or int(paired["n_source_splits"]) != 10:
                    raise ValueError("Unexpected conditional bootstrap configuration")
            for name, key in [("gat_minus_model_observed", "split_mean_point_delta"),
                              ("gat_minus_model_bootstrap_mean", "common_draw_bootstrap_mean"),
                              ("conditional_95_low", "ci_low"), ("conditional_95_high", "ci_high")]:
                row[name] = paired.get(key, "")
            output.append(row)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--check", action="store_true", help="Compare with the released Table S12.")
    action.add_argument("--output", type=Path, help="Write a new CSV file.")
    args = parser.parse_args()
    rows = assemble(args.root)
    if args.check:
        current = read(args.root / "manuscript/source_data/figure5_model_comparison.csv")
        if current != [{key: str(value) for key, value in row.items()} for row in rows]:
            raise SystemExit("Released S12 differs from the frozen aggregate sources")
        print(f"Verified {len(rows)} absolute-score rows and 16 paired AUROC contrasts.")
    else:
        with args.output.open("x", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
