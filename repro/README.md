# Minimal reproducibility entry point

This directory contains the reviewer-facing commands for the compact evidence package. It deliberately omits raw cohort data and machine-specific execution details.

## 1. Validate the package

From the repository root:

```bash
PYTHONPATH=src python -m pytest -q
python repro/validate_package.py
```

## 2. Inspect the evidence layer

The primary review files are:

- `results/submission_tables/Table_1_cohort_contract.csv`
- `results/submission_tables/Table_2_primary_benchmark_metrics.csv`
- `results/submission_tables/Supplementary_Table_1_known_truth_controls.csv`
- `manuscript/source_data/`
- `manuscript/registry/figure_table_source_registry.csv`
- `results/reconciled_evidence/` — compact truth-aware controls, absolute performance and paired bootstrap summaries.

These files contain compact summaries only. They do not contain protected sample identifiers or raw abundance matrices.

## 3. Full rerun boundary

A full rerun requires permitted HMP2/IBDMDB and PRISM inputs, the declared feature mapping and the optional graph dependencies. The manuscript reports the preprocessing boundary, participant grouping, directional transfer and uncertainty estimands; the repository keeps those components in `src/` and `scripts/` so they can be rerun without copying restricted data into the repository.

## 4. Reproducibility principle

Every reported figure or table should be regenerated from a machine-readable source table. Do not edit values manually in figure files or manuscript text.
