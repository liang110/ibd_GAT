# Reproducibility entry point

This directory documents the lightweight checks for the review release. Raw cohort data and machine-specific run products are not distributed.

## Validate

From the repository root:

```bash
PYTHONPATH=src python -m pytest -q
python repro/validate_package.py
python repro/audits/validate_truth_acceptance.py --ledger results/reconciled_evidence/truth_acceptance/per_seed_acceptance.csv
python repro/audits/validate_evidence_registry.py --root .
python repro/audits/validate_topology_constraints.py \
  --edge-audit results/reconciled_evidence/topology_constraints/formal_edge_constraint_audit.csv \
  --hash-audit results/reconciled_evidence/topology_constraints/formal_gateb_hash_audit.csv \
  --status results/reconciled_evidence/topology_constraints/topology_audit_status.json
python repro/audits/validate_reconstructed_edge_lists.py \
  --archive releases/formal_topology_edge_export_20261114.zip
```

## Evidence boundary

`manuscript/source_data/` and `results/submission_tables/` contain compact figure/table sources. `results/reconciled_evidence/` contains the paired bootstrap summaries, topology summaries, truth-acceptance ledger and formal edge audit. Full reruns require permitted HMP2/IBDMDB and PRISM inputs and the optional graph dependencies.

The public release excludes raw abundance matrices, participant identifiers, full prediction ledgers, checkpoints and scheduler logs. Formal edge lists remain in the single release archive because they directly audit the declared graph interventions.
