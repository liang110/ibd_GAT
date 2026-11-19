"""Rebuild Supplementary Tables S14–S15 from the archived GAT and random-forest probabilities."""

import argparse
import csv
import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
from statistics import fmean


MODELS = ("no_prior_gat", "random_forest")
FORWARD = "HMP2_to_PRISM"
REVERSE = "PRISM_to_HMP2"
ORIGINS = ("discovery", "validation")
STAGES = ("record_weighted_loss", "participant_equal_record_loss", "participant_mean_probability")
TABLES = ("pooled_auroc_components.csv", "participant_brier_components.csv")
ROOT = Path(__file__).resolve().parents[1]


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as file:
        yield from csv.DictReader(file)


def read_predictions(path):
    metadata = {FORWARD: {}, REVERSE: {}}
    fits = {FORWARD: defaultdict(dict), REVERSE: defaultdict(dict)}
    for row in read_csv(path):
        direction, model = row["direction"], row["model"]
        if row["experiment_id"] != "reconciled_transfer" or model not in MODELS or direction not in fits:
            continue
        fit, sample = row["fit_id"], row["sample_id"]
        participant, origin = row["participant_id"], row["target_origin"]
        outcome, score = int(row["outcome"]), float(row["probability"])
        identity = (participant, origin, outcome)
        if not all((fit, sample, participant, origin)) or outcome not in (0, 1) or not 0 <= score <= 1:
            raise ValueError(f"Invalid probability, outcome or identifier for {direction}, {fit}, {sample}")
        if metadata[direction].setdefault(sample, identity) != identity:
            raise ValueError(f"Target identity differs between predictions for {sample}")
        scores = fits[direction][fit].setdefault(sample, {})
        if model in scores:
            raise ValueError(f"Duplicate prediction for {fit}, {model}, {sample}")
        scores[model] = score
    for direction, by_fit in fits.items():
        if len(by_fit) != 10:
            raise ValueError(f"Expected 10 archived fits in {direction}, found {len(by_fit)}")
        for fit, panel in by_fit.items():
            if panel.keys() != metadata[direction].keys() or any(set(s) != set(MODELS) for s in panel.values()):
                raise ValueError(f"Incomplete paired target panel for {fit}")
        people = {}
        for participant, origin, outcome in metadata[direction].values():
            if people.setdefault(participant, (origin, outcome)) != (origin, outcome):
                raise ValueError(f"Inconsistent recruitment origin or outcome for {participant}")
        sizes = (len(metadata[direction]), len(people), sum(v[2] for v in metadata[direction].values()),
                 sum(v[1] for v in people.values()))
        if sizes != ({FORWARD: (220, 220, 164, 164), REVERSE: (1626, 130, 1199, 103)}[direction]):
            raise ValueError(f"Target panel counts differ in {direction}: {sizes}")
    return fits, metadata


def pair_auroc(case_scores, control_scores):
    """Mean case-control concordance with half credit for equal probabilities."""
    controls = sorted(control_scores)
    doubled_credits = sum(bisect_left(controls, case) + bisect_right(controls, case) for case in case_scores)
    return doubled_credits / (2 * len(case_scores) * len(controls))


def decompose_auroc(fits, metadata):
    groups = defaultdict(list)
    for sample, (_, origin, outcome) in metadata.items():
        groups[(origin, outcome)].append(sample)
    expected = {("discovery", 1): 121, ("discovery", 0): 34, ("validation", 1): 43, ("validation", 0): 22}
    if {k: len(v) for k, v in groups.items()} != expected:
        raise ValueError("The four forward recruitment/outcome groups do not match the archived panel")
    n_cases = sum(len(groups[origin, 1]) for origin in ORIGINS)
    n_controls = sum(len(groups[origin, 0]) for origin in ORIGINS)
    n_pairs = n_cases * n_controls
    result = []
    for case_origin in ORIGINS:
        for control_origin in ORIGINS:
            cases, controls = groups[case_origin, 1], groups[control_origin, 0]
            means = [fmean(pair_auroc([panel[s][model] for s in cases], [panel[s][model] for s in controls])
                           for panel in fits.values()) for model in MODELS]
            pairs = len(cases) * len(controls)
            weight, delta = pairs / n_pairs, means[0] - means[1]
            result.append(dict(case_origin=case_origin, control_origin=control_origin, n_cases=len(cases),
                               n_controls=len(controls), case_control_pairs=pairs, pair_weight=weight,
                               n_fits=len(fits), gat_auroc=means[0], rf_auroc=means[1],
                               gat_minus_rf=delta, weighted_delta=weight * delta))
    return result


def decompose_brier(fits, metadata):
    members = defaultdict(list)
    for sample, (participant, _, outcome) in metadata.items():
        members[(participant, outcome)].append(sample)
    by_stage = defaultdict(list)
    for panel in fits.values():
        for model in MODELS:
            losses = [(scores[model] - metadata[sample][2]) ** 2 for sample, scores in panel.items()]
            within_person, person_mean = [], []
            for (_, outcome), samples in members.items():
                scores = [panel[s][model] for s in samples]
                within_person.append(fmean((p - outcome) ** 2 for p in scores))
                person_mean.append((fmean(scores) - outcome) ** 2)
            for stage, value in zip(STAGES, (fmean(losses), fmean(within_person), fmean(person_mean))):
                by_stage[stage, model].append(value)
    result, previous = [], None
    for stage in STAGES:
        gat, rf = [fmean(by_stage[stage, model]) for model in MODELS]
        delta = gat - rf
        result.append(dict(stage=stage, n_records=len(metadata), n_participants=len(members), n_fits=len(fits),
                           gat_brier=gat, rf_brier=rf, gat_minus_rf=delta,
                           change_from_previous="" if previous is None else delta - previous))
        previous = delta
    return result


def verify(tables, summary):
    forward, reverse = tables
    archived = {(r["direction"], r["target_stratum"], r["analysis_level"], r["metric"]): r
                for r in read_csv(summary) if (r["candidate"], r["reference"]) == MODELS}

    def equal(a, b):
        if not math.isclose(float(a), float(b), rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"Numerical mismatch: {a} versus {b}")

    def observed(key):
        return archived[key]["split_mean_point_delta"]

    if [(r["case_origin"], r["control_origin"]) for r in forward] != [(a, b) for a in ORIGINS for b in ORIGINS]:
        raise ValueError("Missing or repeated source-pair group")
    for expected_pairs, r in zip((4114, 2662, 1462, 946), forward):
        equal(r["n_fits"], 10)
        equal(r["case_control_pairs"], expected_pairs)
        equal(int(r["n_cases"]) * int(r["n_controls"]), expected_pairs)
        equal(r["pair_weight"], expected_pairs / 9184)
        equal(float(r["gat_auroc"]) - float(r["rf_auroc"]), r["gat_minus_rf"])
        equal(float(r["pair_weight"]) * float(r["gat_minus_rf"]), r["weighted_delta"])
        if r["case_origin"] == r["control_origin"]:
            equal(r["gat_minus_rf"], observed((FORWARD, r["case_origin"], "visit", "roc_auc")))
    equal(sum(float(r["weighted_delta"]) for r in forward), observed((FORWARD, "pooled", "visit", "roc_auc")))
    if [r["stage"] for r in reverse] != list(STAGES) or reverse[0]["change_from_previous"] != "":
        raise ValueError("Unexpected reverse Brier stages")
    for i, r in enumerate(reverse):
        if (int(r["n_records"]), int(r["n_participants"]), int(r["n_fits"])) != (1626, 130, 10):
            raise ValueError("Unexpected reverse denominators")
        equal(float(r["gat_brier"]) - float(r["rf_brier"]), r["gat_minus_rf"])
        if i > 0:
            equal(float(r["gat_minus_rf"]) - float(reverse[i-1]["gat_minus_rf"]), r["change_from_previous"])
    for r, unit in ((reverse[0], "visit"), (reverse[-1], "participant")):
        equal(r["gat_minus_rf"], observed((REVERSE, "pooled", unit, "brier")))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--predictions", type=Path, help="Study prediction_ledger_long.csv")
    modes.add_argument("--check", action="store_true", help="Verify packaged tables against the archived observed contrasts")
    parser.add_argument("--table-dir", type=Path, default=ROOT / "manuscript/source_data")
    parser.add_argument("--summary", type=Path, default=ROOT / "results/p0_rechecked_20261117/primary_common_draw_summary.csv")
    args = parser.parse_args()
    if args.check:
        tables = [list(read_csv(args.table_dir / name)) for name in TABLES]
    else:
        fits, metadata = read_predictions(args.predictions)
        tables = [decompose_auroc(fits[FORWARD], metadata[FORWARD]), decompose_brier(fits[REVERSE], metadata[REVERSE])]
    verify(tables, args.summary)
    if args.predictions:
        args.table_dir.mkdir(parents=True, exist_ok=True)
        for filename, rows in zip(TABLES, tables):
            with (args.table_dir / filename).open("w", newline="", encoding="utf-8") as file:
                writer = csv.DictWriter(file, fieldnames=rows[0].keys(), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
    print("Verified four source-pair AUROCs, three reverse Brier stages and five archived observed contrasts.")


if __name__ == "__main__":
    main()
