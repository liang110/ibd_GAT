# Reconciled XGBoost ledger

This directory contains the current-contract XGBoost predictions and summaries used for manuscript-facing comparisons. The ledger has ten participant-grouped source panels per direction. Hyperparameter selection is source-only within each panel; the selected configuration is refit on all rows of that source panel before locked external scoring. This is a ten-panel prediction ledger with a full-source refit inside each panel, not a single pooled model and not a prospective preregistration.

Files:

- `metrics_by_split.csv`: split-level metrics recomputed from the canonical 18,460-row ledger.
- `common_draw_summary.csv`: paired common-draw summaries used by the manuscript.
- `origin_participant_summary.csv`: target-origin and participant-scoring summaries.
- `gat_xgb_split_deltas.csv`: split-level GAT−XGBoost contrasts.
- `xgboost_metrics_by_split_recomputed.csv`: metrics recomputed directly from the ledger.
- `xgboost_common_draw_summary.csv`: canonical paired common-draw summaries.
- `xgboost_common_draw_deltas.csv`: draw-level deltas used for the pooled summaries.
- `xgboost_origin_participant_summary.csv`: target-origin and participant-scoring summary source.
- `xgboost_vs_gat_split_metric_deltas.csv`: per-split GAT-minus-XGBoost metric contrasts reconstructed against the final No-Prior GAT ledger.
- `xgboost_run_manifest.json`: contract, source hashes, training scope and validation checks.

Target labels are present only for locked scoring and summary validation; they were not used for feature selection, hyperparameter selection or refitting.

The full long-format prediction ledger is retained in the accompanying evidence archive rather than duplicated in the public repository.
