# Reconciled evidence tables

These compact CSV files support the manuscript's evidence audit. They contain
summary statistics only and no raw cohort matrix, clinical identifier or machine
execution record.

- `synthetic_truth_summary.csv`: truth-aware selector and end-to-end control summaries. No-signal scenarios report disease recovery as missing and list false-positive candidate selection separately.
- `absolute_performance_summary.csv`: absolute AUROC, average precision and Brier summaries across source splits.
- `paired_bootstrap_reconciled_by_split.csv`: paired participant-draw contrasts with 5,000 draws per source split.
- `topology_split_stability.csv`, `topology_split_stability.json`, `topology_split_stability.md`: split-level real-versus-rewired topology contrasts and their summary.
- `selector_audit_summary.csv`, `selector_audit_summary.md`: truth-aware selector stress-test summaries; disease TP/FP fields are used only when a declared truth set exists.
- `prism_discovery_origin_sensitivity/`: discovery-only PRISM-to-HMP2 source-composition sensitivity (`summary.csv`, `metrics.csv`, `manifest.json`).

The tables are derived from the permitted cohort analysis outputs and should be
used together with the manuscript figure registry. Brier is lower-is-better.
