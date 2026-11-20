"""Recompute acceptance from each seed's target AUROC and feature recovery."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def validate_truth(ledger_path, rules_path, summary_path):
    runs = pd.read_csv(ledger_path)
    rules = json.loads(rules_path.read_text(encoding="utf8"))["end_to_end_truth_panel"]
    positive = rules["positive_rule"]
    negative = rules["negative_rule"]
    reported = pd.read_csv(summary_path).set_index(["experiment", "scenario"], verify_integrity=True)
    scenarios = set(rules["observed_acceptance"])
    keys = {(exp, s) for exp in ("end_to_end_truth_panel", "selector_audit") for s in scenarios}
    if set(reported.index) != keys or set(runs.groupby(["experiment","scenario"]).groups) != keys:
        raise ValueError("Missing or unexpected true-signal control")
    for key, group in runs.groupby(["experiment", "scenario"]):
        counts = reported.loc[key]
        if len(group) != rules["seeds_per_scenario"] or group.seed.nunique() != len(group):
            raise ValueError(f"Seed coverage: {key}")
        if key[1] in {"batch_only","null_control"}:
            passes = group.target_roc_auc.between(negative["auroc_min_inclusive"],negative["auroc_max_inclusive"])
            valid_rules = (group.acceptance_logic.eq("inclusive_interval") & group.truth_denominator.isna()
                           & group.acceptance_rule.eq("negative_auc_0.35_to_0.65"))
        else:
            passes = ((group.target_roc_auc >= positive["auroc_min_inclusive"]) &
                      (group.n_disease_features_recovered >= positive["recovered_features_min_inclusive"]))
            valid_rules = (group.acceptance_logic.eq(positive["logic"]) &
                           group.truth_denominator.eq(positive["truth_denominator"]) &
                           group.acceptance_rule.eq("positive_auc_0.70_and_recovery_3_of_6"))
        if not valid_rules.all() or not passes.equals(group.acceptance_pass):
            raise ValueError(f"Seed-level rule or decision: {key}")
        if counts.pass_count != passes.sum() or counts.n_seeds != len(group):
            raise ValueError(f"Reported acceptance count: {key}")
        if not np.isclose(counts.mean_auroc,group.target_roc_auc.mean(),atol=1e-12,rtol=0):
            raise ValueError(f"Reported target AUROC: {key}")
        if key[0] == "end_to_end_truth_panel" and rules["observed_acceptance"][key[1]] != f"{passes.sum()}/{len(group)}":
            raise ValueError(f"Acceptance rule summary: {key}")
    print(f"Verified {len(runs)} seed-level decisions and {len(reported)} scenario summaries.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    root = parser.parse_args().root
    validate_truth(root / "results/controls/per_seed_acceptance.csv",
                   root / "results/controls/truth_acceptance_rules.json",
                   root / "results/comparisons/known_truth_acceptance_summary.csv")
