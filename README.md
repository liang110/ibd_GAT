# IBD microbiome transfer framework

Compact code and evidence sources for **A bidirectional cross-cohort benchmark of graph-attention and tabular models for microbiome-based classification of inflammatory bowel disease**.

The release implements a frozen feature contract, source-only preprocessing, participant-aware uncertainty, bidirectional external transfer and same-architecture topology controls. It evaluates a protocol under declared data and graph conditions; it does not claim universal graph-model superiority or clinical utility.

## Repository layout

- `src/ibd_microbiome_framework/` — reusable contracts, preprocessing, splits, models and metrics.
- `scripts/` and `workflows/` — analysis entry points and scheduler templates.
- `manuscript/source_data/` and `manuscript/registry/` — compact figure/table source tables.
- `results/submission_tables/` — manuscript table sources.
- `results/reconciled_evidence/` — compact audit tables for the main claims.
- `repro/audits/` — fail-closed validators for evidence, truth acceptance and formal topology constraints.
- `tests/` — lightweight contract and synthetic-control tests.

Raw cohort matrices, clinical identifiers, full prediction ledgers, checkpoints and scheduler logs are intentionally excluded. Public inputs should be obtained from their original repositories under their access terms.

## Minimal setup

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
PYTHONPATH=src python -m pytest -q
```

The graph branch requires a compatible PyTorch/PyTorch-Geometric installation; these wheels are environment-specific and are not pinned in the lightweight package.

## Evidence validation

From the repository root:

```bash
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

The formal topology export contains 140 edge lists (20 real, 20 identity and 100 degree-preserving rewired graphs) in `releases/formal_topology_edge_export_20261114.zip`. All official edge hashes and selected-feature hashes match; non-identity graphs have no self-loops or duplicate edges, and the former 24 local reconstruction mismatches are documented in `results/reconciled_evidence/topology_constraints/mismatch_resolution.csv`.

## Release scope

The repository keeps only machine-readable files that directly support a manuscript table, figure, truth-panel acceptance rule or topology audit. Raw matrices, participant-level predictions and transient run products remain outside the public review release. The formal edge lists are retained because they are direct evidence for the topology claim rather than intermediate training output.

## Citation

Please cite the accompanying manuscript and the references in `manuscript/clean/References.md`.

