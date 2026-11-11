#!/bin/bash
#SBATCH --job-name=ibd_synthetic_sweep
#SBATCH --partition=normal,normal1,normal2,normal3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/07_synthetic_sweep_%j.out
#SBATCH --error=logs/07_synthetic_sweep_%j.err

set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
background=$(find artifacts/runs -path '*/hmp2_shared_microbiome_proportions.csv' | sort | tail -n 1)
run_id="synthetic_sweep_$(date +%Y%m%d_%H%M%S)"
mkdir -p "artifacts/runs/${run_id}" logs
PYTHONPATH=src python scripts/run_synthetic_sweep.py --background "$background" --output-dir "artifacts/runs/${run_id}"

