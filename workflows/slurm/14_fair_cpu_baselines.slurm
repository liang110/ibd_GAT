#!/bin/bash
#SBATCH --job-name=ibd_fair_baselines
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/14_fair_baselines_%j.out
#SBATCH --error=logs/14_fair_baselines_%j.err

set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
contract=$(find artifacts/runs -maxdepth 1 -type d -name 'contract_*' | sort | tail -n 1)/data_contract
run_id="fair_baselines_$(date +%Y%m%d_%H%M%S)"
out="artifacts/runs/${run_id}"
mkdir -p "$out/bootstrap" "$out/paired" logs
PYTHONPATH=src python scripts/run_fair_cpu_baselines.py \
  --contract-dir "$contract" \
  --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
  --output-dir "$out" \
  --n-splits 5
PYTHONPATH=src python scripts/bootstrap_transfer_uncertainty.py \
  --predictions "$out/fair_baseline_external_predictions.csv" \
  --output-dir "$out/bootstrap" \
  --n-bootstrap 500
PYTHONPATH=src python scripts/paired_model_bootstrap.py \
  --predictions "$out/fair_baseline_external_predictions.csv" \
  --output-dir "$out/paired" \
  --reference elastic_net \
  --n-bootstrap 500

