from __future__ import annotations

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "artifacts" / "reconciled_submission_20261113" / "synthetic_truth_ledger.csv"
OUT = ROOT / "artifacts" / "reconciled_submission_20261113"


def main() -> None:
    df = pd.read_csv(LEDGER)
    audit = df.loc[df.experiment.eq("selector_audit")].copy()
    rows = []
    for scenario, g in audit.groupby("scenario", sort=False):
        truth_defined = bool(g.disease_truth_defined.iloc[0])
        rows.append({
            "experiment": "selector_audit",
            "scenario": scenario,
            "n_seeds": int(g.seed.nunique()),
            "selected_budget_mean": float(g.n_selected_features.mean()),
            "disease_truth_defined": truth_defined,
            "disease_tp_mean": float(g.disease_tp.mean()) if truth_defined else None,
            "disease_fn_mean": float(g.disease_fn.mean()) if truth_defined else None,
            "precision_mean": float(g.disease_precision.mean()) if truth_defined else None,
            "recall_mean": float(g.disease_recall.mean()) if truth_defined else None,
            "batch_fp_mean": float(g.batch_fp.mean()),
            "batch_fpr_mean": float(g.batch_fpr.mean()),
            "random_same_budget_expected_disease_tp": float(g.random_same_budget_expected_disease_tp.mean()),
            "random_same_budget_expected_batch_fp": float(g.random_same_budget_expected_batch_fp.mean()),
            "designated_feature_overlap_mean": float(g.disease_tp.mean()) if not truth_defined else None,
            "target_auroc_mean": float(g.target_roc_auc.mean()),
            "acceptance_rate": float(g.acceptance_pass.mean()),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "selector_audit_summary.csv", index=False)
    (OUT / "selector_audit_summary.md").write_text(
        "# Synthetic selector audit\n\n"
        "TP/FP metrics are reported only for scenarios with a declared disease-feature truth set. "
        "For null and batch-only scenarios, the designated-feature overlap is a random overlap with inert IDs and is not a causal disease TP.\n",
        encoding="utf-8",
    )
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
