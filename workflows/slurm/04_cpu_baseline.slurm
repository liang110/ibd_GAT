#!/bin/bash
#SBATCH --job-name=ibd_cpu_baseline
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/04_cpu_baseline_%j.out
#SBATCH --error=logs/04_cpu_baseline_%j.err

set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
contract=$(find artifacts/runs -maxdepth 1 -type d -name 'contract_*' | sort | tail -n 1)/data_contract
preprocess=$(find artifacts/runs -maxdepth 1 -type d -name 'preprocess_*' | sort | tail -n 1)
run_id="baseline_$(date +%Y%m%d_%H%M%S)"
mkdir -p "artifacts/runs/${run_id}" logs
PYTHONPATH=src python scripts/run_cpu_baseline.py --contract-dir "$contract" --preprocess-dir "$preprocess" --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} --output-dir "artifacts/runs/${run_id}" --n-splits 5

