# IBD microbiome transfer framework

Code and compact evidence tables for the manuscript **A bidirectional cross-cohort benchmark of graph-attention and tabular models for microbiome-based classification of inflammatory bowel disease**.

The contribution is an auditable analysis protocol: a frozen microbiome feature contract, source-only preprocessing, participant-aware uncertainty, bidirectional external transfer and same-architecture topology controls. The repository does not claim a universally superior graph model or clinical utility.

## Repository contents

- `src/ibd_microbiome_framework/` — reusable contracts, preprocessing, splits, models and metrics.
- `scripts/` — lightweight validation, analysis and export entry points.
- `configs/` — cohort and experiment configuration templates.
- `results/submission_tables/` — compact tables supporting the manuscript claims.
- `manuscript/clean/` — section-level manuscript text, captions and references.
- `manuscript/source_data/` — small figure source tables and registries.
- `tests/` — contract, synthetic-control and model tests.

Raw cohort matrices, clinical identifiers and transient run products are not distributed. Obtain public inputs from their original repositories under the applicable access terms.

## Minimal setup

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Run the lightweight validation suite:

```bash
PYTHONPATH=src python -m pytest -q
```

The graph branch requires a compatible PyTorch and PyTorch-Geometric installation. Those packages are intentionally not pinned here because the correct wheel depends on the local CPU/CUDA runtime.

## Evidence regeneration

The compact tables in `results/submission_tables/` and the figure source files in `manuscript/source_data/` are the preferred review entry point. Use the scripts documented in `repro/README.md` to validate the package and regenerate the summary tables from permitted inputs. Full cohort reruns require the original public data and are intentionally separated from the lightweight checks.

## Scientific scope

The primary endpoint is binary IBD versus non-IBD classification. Subtype, disease activity, treatment response, causal microbial interactions and prospective clinical utility are outside the primary estimand. The exact-name cross-cohort metabolite intersection is zero; the reported transfer analysis is microbiome-only.

## Citation

Please cite the accompanying manuscript and the references listed in `manuscript/clean/References.md`. The software is distributed as a review-release package; publication metadata will be added with the final archival release.

