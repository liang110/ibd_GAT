# Reproducibility runbook

This document is the short operational companion to the root README. It separates the checks that can run without raw cohort data from the full HPC rerun.

## 1. Local code and synthetic checks

```bash
python -m pip install -e ".[dev]"
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=src python -m pytest -q
PYTHONPATH=src python scripts/run_synthetic_validation.py --help
PYTHONPATH=src python scripts/run_synthetic_control_suite.py --help
```

These checks validate explicit label contracts, participant-aware splits, source-only compositional preprocessing and the known-truth control machinery. They do not establish biological validity.

## 2. Gate0 and contract

Set `IBD_LEGACY_ROOT` to a local copy of the public input files. The expected filenames are documented in `configs/data_catalog.json` and the Gate0 launcher. Run the input audit before any model fitting. The contract output must include the feature registry, cohort manifests, shared feature candidates, library sums, input inventory and SHA-256 values.

```bash
export IBD_GAT_ROOT="$PWD"
export IBD_LEGACY_ROOT=/path/to/acquired/public-input-root
export PYTHONPATH="$IBD_GAT_ROOT/src"
python scripts/build_data_contract.py \
  --legacy-root "$IBD_LEGACY_ROOT" \
  --output-dir artifacts/runs/gate0_contract
```

## 3. Primary transfer and audits

The primary graph run is `scripts/run_no_prior_gat_transfer.py`. It uses a frozen shared feature list, source-only closure/CLR/scaling/MI selection, participant-grouped source splits, source-derived topology and a No-Prior GATv2 branch paired with elastic-net and random-forest references. The target matrix is transformed but never used to fit labels, thresholds, feature rankings or topology.

The companion jobs are:

- `scripts/run_prism_origin_transfer.py` and `scripts/run_origin_paired.py`: discovery/validation origin audit;
- `scripts/run_repeated_leakage_controls.py`: repeated-participant and planted leakage controls;
- `scripts/run_synthetic_control_suite.py` and `scripts/run_domain_stable_synthetic_sweep.py`: known-truth stress tests;
- `scripts/audit_topology_runs.py` and `scripts/export_topology_edges.py`: topology, edge Jaccard and fixed-k diagnostics;
- `scripts/run_capacity_control.py`: approximate no-message comparator;
- `scripts/paired_bootstrap_confirm.py` and `scripts/run_origin_paired.py`: participant-cluster paired uncertainty.

## 4. Outputs and audit trail

Every run should write a manifest containing the Git revision, input/output hashes, feature-contract hash, split identifiers, seeds, topology hash, model parameters and prediction paths. Figure and table source files must be generated from the frozen output tables. Do not edit figure numbers manually after a run.

The public package includes compact source tables and figure registries. Large raw matrices, checkpoints, logs and transient Slurm outputs remain external.
