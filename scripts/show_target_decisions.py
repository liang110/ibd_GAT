"""Print the ten archived GAT-minus-random-forest contrasts in manuscript Box 1."""

import argparse
import csv
import math
import sys
from pathlib import Path


SUMMARY = Path("results/p0_rechecked_20261117/primary_common_draw_summary.csv")
# Target label, stored direction and stratum, scoring level, scored units, people, control units.
TARGETS = (
    ("Combined Franzosa", "HMP2_to_PRISM", "pooled", "visit", 220, 220, 56),
    ("PRISM discovery", "HMP2_to_PRISM", "discovery", "visit", 155, 155, 34),
    ("Netherlands validation", "HMP2_to_PRISM", "validation", "visit", 65, 65, 22),
    ("HMP2", "PRISM_to_HMP2", "pooled", "visit", 1626, 130, 427),
    ("HMP2", "PRISM_to_HMP2", "pooled", "participant", 130, 130, 27),
)


def assemble(root):
    with (root / SUMMARY).open(encoding="utf-8-sig", newline="") as handle:
        contrasts = [r for r in csv.DictReader(handle)
                     if r["candidate"] == "no_prior_gat" and r["reference"] == "random_forest"]
    rows = []
    for target, direction, stratum, level, units, people, controls in TARGETS:
        for metric in ("roc_auc", "brier"):
            key = (direction, stratum, level, metric)
            matches = [r for r in contrasts if
                       (r["direction"], r["target_stratum"], r["analysis_level"], r["metric"]) == key]
            if len(matches) != 1:
                raise ValueError(f"{key}: expected one paired contrast, found {len(matches)}")
            row = matches[0]
            fields = ("n_source_splits", "n_bootstrap", "n_target_samples", "n_target_participants",
                      "class_0_count", "class_1_count")
            if tuple(int(row[f]) for f in fields) != (10, 5000, units, people, controls, units-controls):
                raise ValueError(f"{key}: fixed-fit count, bootstrap size or target denominator changed")
            estimate, lo, hi = (float(row[f]) for f in ("split_mean_point_delta", "ci_low", "ci_high"))
            excludes_zero = lo > 0 or hi < 0
            if not all(math.isfinite(x) and -1 <= x <= 1 for x in (estimate, lo, hi)) or lo > hi:
                raise ValueError(f"{key}: invalid contrast or interval")
            if row["ci_excludes_zero"] != str(excludes_zero):
                raise ValueError(f"{key}: interval and exclusion flag disagree")
            choice = "unresolved"
            if excludes_zero:
                choice = "GAT" if (lo > 0) == (metric == "roc_auc") else "RF"
            rows.append(dict(
                target=target, scoring_unit="record" if level == "visit" else "participant",
                scored_units=units, participants=people, metric=metric,
                observed_gat_minus_rf=row["split_mean_point_delta"],
                conditional_95_low=row["ci_low"], conditional_95_high=row["ci_high"],
                supported_procedure=choice,
            ))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1],
                        help="Repository or Supplementary Software 1 root.")
    args = parser.parse_args()
    rows = assemble(args.root)
    writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
