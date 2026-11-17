# IBD microbiome transfer protocol

Compact code and evidence sources for **An auditable framework for directional transfer of graph and tabular microbiome classifiers in inflammatory bowel disease**.

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

Figure and table sources are identified by analysis and scoring unit; current P0-rechecked summaries are listed below. The repository retains only files that directly support manuscript tables, figures and audit claims. It does not redistribute raw cohort data or transient scheduler outputs.

## Formal topology evidence
The formal topology results are bound to the compact current-contract exports under `results/reconciled_evidence/topology_constraints/` and the runtime-boundary manifest `formal_topology_runtime_manifest.json`. The exported Real, Rewired and Identity objects are pre-convolution input configurations. The generic transfer launcher is retained as a reusable example and is not presented as the exact generator of every formal artifact.


## Data-contract outputs

The contract builder writes `cohort_manifest_<cohort>_<modality>.csv` for each input modality. For compatibility with existing microbiome analysis commands, `cohort_manifest_<cohort>.csv` is an exact alias of the microbiome manifest. Use a fresh output directory when rebuilding contracts. Raw input matrices are never modified.

HMP2/IBDMDB resources are available through https://ibdmdb.org/; the original PRISM sequence project is NCBI BioProject PRJNA400072. The processed HMP2 manifest used here does not contain a raw-accession-to-local-sample crosswalk.

## P0-rechecked analysis summaries

Compact aggregate tables are available in `results/p0_rechecked_20261117/`. They cover the principal common-target paired summaries, XGBoost scoring-unit metrics, direct PRISM provenance transfer, taxon-aligned topology overlap, graph interventions and known-truth acceptance. Values and intervals use the units documented in the tables and manuscript. The data do not include individual participant IDs or raw prediction ledgers.

Older source-table and registry entries predate this recheck and may retain superseded scoring-unit, topology-index or figure-coordinate values. Use the P0-rechecked aggregate tables for the revised manuscript; the existing older validators do not by themselves validate these new tables.

Run `python results/p0_rechecked_20261117/validate_p0_rechecked.py` for read-only checks of these compact tables. This check does not regenerate predictions or bootstrap samples.

In the synthetic acceptance table, `mean_hits_among_six_designated_candidates` is a restricted candidate-hit count and is not the overall false-positive count. Topology-overlap summaries include only the Spearman edge sets whose node order and formal real-edge exports were independently reconciled.

The XGBoost procedure uses source-only inner selection followed by a final refit on the complete source cohort. The primary GAT is trained on each grouped source-training split. Their external targets and candidate registry are shared, but their final training samples and fitted preprocessing matrices differ. Topology intervention full-batch training uses each source-training partition, not a complete-cohort refit. Participant bootstrap intervals are conditional on the fitted source models.
