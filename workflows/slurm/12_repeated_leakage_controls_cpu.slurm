#!/bin/bash
#SBATCH --job-name=ibd_repeat_leakage
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/12_repeated_leakage_%j.out
#SBATCH --error=logs/12_repeated_leakage_%j.err

set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
background=$(find artifacts/runs -path '*/hmp2_shared_microbiome_proportions.csv' | sort | tail -n 1)
run_id="repeated_leakage_$(date +%Y%m%d_%H%M%S)"
mkdir -p "artifacts/runs/${run_id}" logs
PYTHONPATH=src python scripts/run_repeated_leakage_controls.py \
  --background "$background" \
  --output-dir "artifacts/runs/${run_id}"

