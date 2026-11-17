# Compact evidence summaries

These tables accompany the manuscript's fixed two-cohort analyses. They contain aggregate statistics, not patient-level predictions or identifiers. The `validate_p0_rechecked.py` script in this directory checks the schemas and selected reported values with the Python standard library:

```bash
python results/p0_rechecked_20261117/validate_p0_rechecked.py
# If the tables are installed elsewhere:
python validate_p0_rechecked.py --evidence-dir path/to/public_evidence
```

| File | Contents |
|---|---|
| `primary_common_draw_summary.csv` | GAT–elastic net, random forest and no-message MLP, including pooled and PRISM target-stratum results. |
| `xgboost_metrics_by_scoring_unit.csv` | Ten-fit means and source-fit standard deviations for GAT and XGBoost by target and scoring unit. |
| `xgboost_paired_common_draw.csv` | Paired GAT–XGBoost differences, including target-stratum summaries. |
| `prism_direct_transfer_metrics.csv` | Direct discovery–validation transfer for four learners, mean and source-fit SD. |
| `prism_direct_transfer_paired.csv` | Direct-transfer paired point differences and conditional target-participant intervals. |
| `topology_intervention_paired.csv` | Real-versus-rewired or identity intervention summaries, by direction and scoring unit. |
| `topology_taxon_aligned_summary.csv` | Pairwise directed Spearman edge overlap after taxon-ID alignment. |
| `known_truth_acceptance_summary.csv` | End-to-end and separate selector audit scenario summaries; negative-scenario disease recall is undefined. |

`analysis_level` or `scoring_unit` denotes whether metrics are computed over visits or participant-averaged predictions. Where specified, the resampling unit is the participant. `split_mean_point_delta` (or `point_mean_delta_auroc`) is the mean of the observed split-wise differences; `common_draw_bootstrap_mean` (or `bootstrap_mean_delta_auroc`) is the mean across paired target-bootstrap draws. They are different quantities. The percentile intervals describe target resampling conditional on the fixed source fits; they do not quantify retraining or transport to unseen cohorts. XGBoost was finally refitted on all source-cohort rows after source-only selection, whereas the primary GAT used each source training partition.

The validator checks declared values, row/column structure, selected denominator and unit labels and internal point-difference arithmetic. It does **not** retrain models, rerun the bootstrap or independently reconstruct all statistics from individual predictions. Input acquisition, privacy and reuse remain subject to the original data-provider terms.
