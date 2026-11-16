# IBD microbiome transfer protocol

Compact code and evidence sources for **A bidirectional cross-cohort benchmark of graph-attention and tabular models for microbiome-based classification of inflammatory bowel disease**.

The release implements a frozen microbiome feature contract, source-only preprocessing, participant-aware uncertainty, bidirectional transfer and controlled topology comparisons. Its contribution is an auditable evaluation protocol under declared data and graph conditions.

## Repository layout

- `src/ibd_microbiome_framework/` — reusable contracts, preprocessing, splits, models and metrics.
- `scripts/` and `workflows/` — analysis entry points and scheduler templates.
- `manuscript/source_data/` and `manuscript/registry/` — compact figure/table source tables.
- `results/submission_tables/` — manuscript table sources.
- `repro/audits/` — fail-closed validators for evidence, truth acceptance and topology constraints.
- `tests/` — lightweight contract and synthetic-control tests.

Raw cohort matrices, clinical identifiers, full prediction ledgers, checkpoints and scheduler logs are excluded. Public inputs are acquired from their original repositories under their access terms.

## Minimal setup

```bash
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
PYTHONPATH=src python -m pytest -q
```

The graph branch requires a compatible PyTorch/PyTorch-Geometric installation; cluster-specific wheels are not pinned in the lightweight package.

## Focused evidence checks

```bash
python repro/validate_package.py
python repro/audits/validate_evidence_registry.py --root .
python repro/audits/validate_evidence_content_20261116.py --root .
python scripts/audit_cohort_independence_20261116.py \
  --source-root path/to/canonical_contract/data_contract \
  --output-dir artifacts/accession_audit
```

The current HMP2/PRISM contract records 1,638 HMP2 manifest observations, 1,626 included observations from 130 participants, and 220 PRISM observations from 220 operational participants. The current Figure 3 source table records 80 synthetic participants, 240 observations, 43 mean participant overlaps under sample-level splitting, zero overlap under participant-grouped splitting, and current-contract mean AUROC values of 0.960 and 0.552182. The HMP2/PRISM audit records zero exact, normalized and substring overlap for sample and participant identifiers in both manifest and included panels. The direct PRISM discovery-to-validation transfer evidence is stored separately under `results/reconciled_evidence/prism_discovery_validation_transfer/`; its recruitment-source confounding remains explicit. The local HMP2 manifest contains no raw accession column, so accession-to-sample mapping remains a provider-level metadata boundary.

## Reproducibility tiers

1. **Local validation:** package tests and synthetic known-truth controls.
2. **Evidence audit:** compact source tables, manifests and figure/table registries.
3. **Full rerun:** optional after acquiring the public cohort inputs and matching the declared environment.

The formal Figure/Table source registry is the sole authoritative mapping for this manuscript; historical contracts remain outside the current registry. The repository retains only files that directly support manuscript tables, figures and audit claims. It does not redistribute raw cohort data or transient scheduler outputs.

## Formal topology evidence
The formal topology results are bound to the compact current-contract exports under `results/reconciled_evidence/topology_constraints/` and the runtime-boundary manifest `formal_topology_runtime_manifest.json`. The exported Real, Rewired and Identity objects are pre-convolution input configurations. The generic transfer launcher is retained as a reusable example and is not presented as the exact generator of every formal artifact.


## Canonical evidence closure (2026-11-17)

The manuscript-facing reconstruction uses the reconciled contract recorded in manuscript/manifests/canonical_contract_manifest_20261117.json (126 shared taxa; HMP2 1,638 manifest and 1,626 included observations; PRISM 220 observations). Compact evidence sources are under esults/reconciled_evidence/xgboost_current_contract/ and esults/reconciled_evidence/topology_contract_freeze_20261117/. Historical source tables are marked in manuscript/source_data/source_data_registry_20261117.csv and are not used for current figures. The canonical clean manuscript source is under manuscript/clean/canonical_20261117/.

