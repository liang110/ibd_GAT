# An auditable framework for directional microbiome transfer in IBD

This repository contains the code and frozen evidence layer for the manuscript **“An auditable framework for directional transferability of relational and tabular microbiome representations in inflammatory bowel disease.”** It is a framework for testing when a relational microbiome representation adds information beyond like-for-like tabular references under an explicitly declared source-to-target direction.

The contribution is the **auditable transfer contract**, not a universally superior GAT implementation. The primary relational branch is a source-only, No-Prior GATv2 representation. The repository is organized so that input harmonization, leakage controls, directional transfer, topology audits, capacity controls and uncertainty are inspectable independently.

> **Release boundary.** Raw HMP2/IBDMDB and PRISM files are not redistributed here. The repository contains acquisition expectations, hashes, manifests, source tables and figure registries. The public repository URL and archival DOI will be added when the authors create the final GitHub release.

## What the framework tests

The analysis has three linked phases:

1. **Contract and validity.** Freeze the shared feature space and denylist, harmonize compositional scale, fit CLR/scaling/MI selection on source-training participants only, and test participant leakage with known-truth controls.
2. **Directional transfer.** Fit on one cohort and evaluate on a locked target cohort in both directions. PRISM discovery and validation are also kept as separate origin strata.
3. **Interpretation boundary.** Use topology nulls, edge Jaccard, the fixed- *k* degree diagnostic, participant-cluster paired bootstrap and a capacity-matched no-message comparator to determine what the result does and does not establish.

The frozen primary contract contains 126 shared microbial features. The exact-name metabolite intersection is zero, so the cross-cohort claim is microbiome-only. HMP2/IBDMDB contributes 1,638 samples from 130 participants; PRISM contributes 220 samples from 220 participants, split into 155 discovery and 65 validation samples. The endpoint is IBD versus non-IBD; subtype, activity and therapy effects are outside the primary estimand.

The principal directional result is conditional: the No-Prior Spearman graph reaches AUROC 0.846 for HMP2 → PRISM, compared with 0.715 for elastic net and 0.768 for random forest. In the reverse PRISM → HMP2 direction, the graph reaches 0.612 while elastic net and random forest reach 0.711 and 0.680. These values are evidence about this declared contract and direction; they are not a universal GAT ranking or a validated biomarker result.

## Repository map

```text
src/ibd_microbiome_framework/   Production package: contracts, harmonization,
                                compositional preprocessing, splits, models and metrics
scripts/                         Reproducible audit, benchmark, transfer, bootstrap and export commands
configs/                         Immutable experiment YAML and data-catalog templates
workflows/slurm/                 Portable Slurm launchers for the final audit sequence
slurm/                           Earlier/compatibility launchers retained for provenance
tests/                            Contract, leakage, synthetic and model tests
docs/                            Architecture, implementation and reviewer-closure documentation
manuscript/clean/                Clean manuscript sections used for the current draft
manuscript/figures/              Figure 1–7 source plots and editable PNG/PDF/SVG outputs
manuscript/source_data/          Figure source CSV/JSON files
manuscript/registry/             Figure/table source registry
manuscript/manifests/            Frozen evidence register
results/submission_tables/       Compact supplementary/result tables behind the manuscript claims
sources/                          Literature and benchmark notes
```

The Word manuscript is intentionally not part of this code release. It remains a submission artifact outside the reproducibility core.

## Installation

The supported runtime is Python 3.10 or newer. For the contract, synthetic and CPU checks:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

For the graph and boosting branches:

```bash
python -m pip install -e ".[model]"
```

PyTorch and PyTorch Geometric must match the CUDA/CPU runtime on the HPC system. The local validation environment used for this release recorded the core versions in `requirements-reproducibility.txt`; the cluster-specific torch wheel should be installed according to the cluster module or conda environment. SIAMCAT is an optional R reference benchmark and is installed by `scripts/install_siamcat.R` in an isolated R library.

## Quick checks

Run the package tests without loading unrelated globally installed pytest plugins:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=src python -m pytest -q
```

Run the synthetic validation entry points before touching public cohort data:

```bash
PYTHONPATH=src python scripts/run_synthetic_validation.py --help
PYTHONPATH=src python scripts/run_synthetic_control_suite.py --help
PYTHONPATH=src python scripts/run_domain_stable_synthetic_sweep.py --help
```

## Reproducing the evidence layers

Set the two paths explicitly on a workstation or HPC login node. The raw inputs must be acquired from their public source repositories under the applicable terms:

```bash
export IBD_GAT_ROOT="$PWD"
export IBD_LEGACY_ROOT=/path/to/acquired/public-input-root
export CONDA_SH="$HOME/miniforge3/etc/profile.d/conda.sh"
export PYTHONPATH="$IBD_GAT_ROOT/src:$PYTHONPATH"
```

Build the immutable contract first:

```bash
python scripts/build_data_contract.py \
  --legacy-root "$IBD_LEGACY_ROOT" \
  --output-dir artifacts/runs/gate0_contract
```

The contract records canonical feature IDs, participant manifests, input SHA-256 hashes, library sums, the shared feature list and zero-mass decisions. All later runs must point to this contract directory. No target labels or target-derived feature statistics may enter a source object.

For a complete Slurm sequence, submit the validity gates before the transfer jobs:

```bash
sbatch workflows/slurm/00_gate0_input_audit_cpu.slurm
sbatch workflows/slurm/01_build_data_contract_cpu.slurm
sbatch workflows/slurm/02_audit_preprocessing_cpu.slurm
sbatch workflows/slurm/03_synthetic_validation_cpu.slurm
sbatch workflows/slurm/10_synthetic_controls_cpu.slurm
sbatch workflows/slurm/12_repeated_leakage_controls_cpu.slurm
sbatch workflows/slurm/13_prism_origin_transfer_cpu.slurm
sbatch workflows/slurm/14_fair_cpu_baselines.slurm
export TOPOLOGY=spearman
sbatch workflows/slurm/22_no_prior_gat_transfer_normal2.slurm
sbatch workflows/slurm/23_paired_bootstrap_confirm.slurm
sbatch workflows/slurm/25_topology_audit.slurm
sbatch workflows/slurm/26_capacity_control.slurm
sbatch workflows/slurm/27_origin_paired_bootstrap.slurm
sbatch workflows/slurm/28_export_topology_edges.slurm
```

The launchers use `IBD_GAT_ROOT`, `IBD_LEGACY_ROOT`, `CONDA_SH` and, for the R reference, `SIAMCAT_RSCRIPT`. They declare `normal,normal1,normal2,normal3` where the job is intended to run across the normal queues. Review every `#SBATCH` resource request against the target cluster before submission.

## Frozen evidence and figure regeneration

The source of each manuscript object is recorded in:

- `manuscript/manifests/frozen_evidence_register_v1.md`
- `manuscript/registry/figure_table_source_registry.csv`
- `manuscript/source_data/`
- `results/submission_tables/`

The final figures are generated from source tables rather than manually edited values. The plotting scripts are in `manuscript/figures/`; they write PNG, PDF and SVG outputs. TIFF files, Word files, transient renderings and HPC logs are excluded from the public code copy.

## Reviewer-facing evidence map

| Reviewer concern | Repository evidence |
|---|---|
| **Reviewer 1: graph branch underperforms externally** | `manuscript/clean/Discussion.md`, directional results, topology controls and the capacity comparator. The graph is framed as conditional and diagnostic, not universally superior. |
| **Reviewer 1: MI lookahead boundary** | `src/ibd_microbiome_framework/preprocessing/compositional.py`, `scripts/run_no_prior_gat_transfer.py`, source-only manifests and `manuscript/clean/Methods.md`. CLR, prevalence, scaling and MI selection are fitted inside source-training participants. |
| **Reviewer 1: baseline equivalence and edge rationale** | Shared denylist and input contract in `src/.../data/contracts.py`, topology scripts, Figure 6 source data and the fixed top- *k*=8 rule. |
| **Reviewer 2: known-truth validation** | `scripts/run_synthetic_control_suite.py`, `scripts/run_domain_stable_synthetic_sweep.py`, `results/submission_tables/Supplementary_Table_1_known_truth_controls.csv`. Synthetic controls validate the analysis procedure, not IBD biology. |
| **Reviewer 2: small heterogeneous cohorts and unclear narrative** | Participant counts, endpoint boundary and three-phase framework in `manuscript/clean/Introduction.md`, `Methods.md`, `Results.md` and `Discussion.md`. |
| **Reviewer 2: graph should be exploratory/contextual** | No-Prior graph is the declared relational branch; tabular baselines remain co-primary references and reverse transfer is reported explicitly. |
| **Reviewer 3: abundance-scale mismatch and shared feature space** | Contract construction, CLR/closure code, Figure 2 source data, 126-feature manifest and explicit zero metabolite intersection. |
| **Reviewer 3: pooled PRISM origins** | `scripts/run_prism_origin_transfer.py`, `scripts/run_origin_paired.py`, Figure 4 source data and discovery↔validation audit outputs. |
| **Reviewer 3: clinical leakage and whole-cohort priors** | Shared denylist, planted leakage controls and the primary No-Prior branch. The legacy whole-cohort biological prior is excluded from the primary claim. |
| **Reviewer 3: uncertainty, PR-AUC and effective sample size** | Participant-grouped splits, participant-cluster bootstrap, AUROC/average precision/Brier outputs, split manifests and supplementary tables. |
| **Reviewer 3: reproducibility** | `pyproject.toml`, `requirements-reproducibility.txt`, portable Slurm launchers, manifests, hashes, source registries and tests. |

The detailed closure notes are in `docs/review/reviewer_audit.md`, `docs/review/final_consistency_audit.md` and `docs/review/reference_integrity_audit.md`.

## Reproducibility tiers

- **Tier 1: code validity.** Run unit tests and synthetic known-truth controls without downloading cohort data.
- **Tier 2: evidence audit.** Inspect the frozen source tables, figure source data, manifests and registry included in this repository.
- **Tier 3: full rerun.** Acquire the public raw inputs, set the environment variables, build Gate0 and submit the Slurm sequence. Full retraining depends on the cluster's PyTorch Geometric/CUDA environment and the original public data versions.

## Data, privacy and provenance

No raw HMP2/IBDMDB or PRISM matrices are stored here. The repository records expected filenames, acquisition boundaries, feature counts, zero-mass decisions and SHA-256 hashes where the frozen evidence package contains them. Raw data must be obtained from the original public repositories and must not be silently replaced by a differently processed table. Identifiers and clinical fields are removed by the shared denylist before any model branch is built.

## Scientific boundaries

This release does **not** claim:

- a universal GAT or GATv2 advantage;
- a validated microbiome biomarker panel or causal ecological interaction;
- prospective clinical utility;
- a third independent biological cohort validation;
- a cross-cohort metabolome transfer result (the exact-name metabolite intersection is zero).

The primary result is a conditional, directional relational gain under a frozen microbiome contract. A tabular model is the appropriate deployment choice when it transfers better.

## Citation and license

Please cite the accompanying manuscript and the methods papers listed in `manuscript/clean/References.md`. A `CITATION.cff` template is included for the final author/DOI update.

The authors have not yet selected a public software license for this local release candidate. Until `LICENSE` is replaced with the authors' approved license, treat the code as all rights reserved and use it for internal research review only.
