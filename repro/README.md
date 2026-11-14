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

## Reviewer-facing audit layer

The compact audit layer exposes only the evidence needed to check the declared protocol:

- `results/reconciled_evidence/topology_constraints/edge_constraint_audit.csv` records per-direction and per-split edge counts, self-loop counts and identity Jaccard values.
- `results/reconciled_evidence/topology_constraints/edge_hash_audit.csv` records edge-set reconstruction hashes and rewiring seeds.
- `results/reconciled_evidence/truth_acceptance/per_seed_acceptance.csv` records the per-seed truth-panel acceptance inputs and outcomes.
- `repro/audits/training_configuration.json` declares the loss, optimizer, batching, seeds and uncertainty units used by the reported runs.

Run the validators from the repository root:

```bash
python repro/audits/validate_truth_acceptance.py --ledger results/reconciled_evidence/truth_acceptance/per_seed_acceptance.csv
python repro/audits/validate_topology_constraints.py --edge-audit results/reconciled_evidence/topology_constraints/edge_constraint_audit.csv --hash-audit results/reconciled_evidence/topology_constraints/edge_hash_audit.csv --status results/reconciled_evidence/topology_constraints/topology_audit_status.json
```

The topology validator is intentionally fail-closed: a hash mismatch is reported rather than silently treated as a successful reconstruction. The compact files contain no raw abundance matrix or participant identifier.


## Reviewer-facing evidence status

`results/reconciled_evidence/evidence_registry.csv` maps Table 3, Figure 6, Supplementary Figure S3 and Supplementary Figure S4 to compact machine-readable sources. `results/reconciled_evidence/topology_constraints/topology_audit_status.json` records the current formal edge-audit limitation; summary metrics are not treated as proof of a complete edge-level reconstruction. The package excludes raw cohort matrices, participant identifiers, full prediction ledgers, checkpoints and scheduler logs.

```bash
python repro/audits/validate_evidence_registry.py --root .
```
