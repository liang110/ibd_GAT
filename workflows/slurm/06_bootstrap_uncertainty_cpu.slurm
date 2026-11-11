#!/bin/bash
#SBATCH --job-name=ibd_transfer_ci
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/06_bootstrap_uncertainty_%j.out
#SBATCH --error=logs/06_bootstrap_uncertainty_%j.err

set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
baseline=$(find artifacts/runs -maxdepth 1 -type d -name 'baseline_*' | sort | tail -n 1)
run_id="bootstrap_$(date +%Y%m%d_%H%M%S)"
mkdir -p "artifacts/runs/${run_id}" logs
PYTHONPATH=src python scripts/bootstrap_transfer_uncertainty.py --predictions "$baseline/cpu_baseline_external_predictions.csv" --output-dir "artifacts/runs/${run_id}" --n-bootstrap 500

