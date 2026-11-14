# Reconciled evidence tables

Compact machine-readable sources for the manuscript claims. These files contain summary statistics and audit metadata, not raw cohort matrices or participant identifiers.

- `absolute_performance_summary.csv`: AUROC, average precision and Brier summaries across source splits.
- `paired_bootstrap_reconciled_by_split.csv`: paired participant-draw contrasts with 5,000 draws per source split.
- `topology_split_stability.csv`: split-level real-versus-rewired topology contrasts.
- `identity_aware_topology/`: compact metrics and paired split contrasts for the identity-aware topology sensitivity analysis.
- `prism_discovery_origin_sensitivity/`: compact metrics and summaries for the PRISM provenance-stratum sensitivity analysis.
- `synthetic_truth_summary.csv` and `truth_acceptance/per_seed_acceptance.csv`: truth-panel summaries and per-seed acceptance records.
- `topology_constraints/`: canonical formal audit tables, hashes, manifests and mismatch resolution.
- `../../releases/formal_topology_edge_export_20261114.zip`: the 140 canonical edge lists and the complete edge-level audit payload.

Brier is lower-is-better. The figure/table registry identifies which source table supports each reported object.
