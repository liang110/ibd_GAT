#!/bin/bash
#SBATCH --job-name=ibd_domain_paired
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/18_domain_paired_%j.out
#SBATCH --error=logs/18_domain_paired_%j.err

set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
domain_run=artifacts/runs/domain_transfer_20260816_050740
baseline_run=artifacts/runs/fair_baselines_20260816_044331
out="${domain_run}/paired_elastic_net"
mkdir -p "$out" logs
PYTHONPATH=src python scripts/prepare_paired_predictions.py \
  --reference-predictions "$baseline_run/fair_baseline_external_predictions.csv" \
  --candidate-predictions "$domain_run/domain_interaction_transfer_predictions.csv" \
  --reference-model elastic_net \
  --output-dir "$out"
PYTHONPATH=src python scripts/paired_model_bootstrap.py \
  --predictions "$out/paired_predictions.csv" \
  --output-dir "$out" \
  --reference elastic_net \
  --n-bootstrap 500

