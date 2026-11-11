#!/bin/bash
#SBATCH --job-name=ibd_synthetic_truth
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/03_synthetic_validation_%j.out
#SBATCH --error=logs/03_synthetic_validation_%j.err

set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
run_id="synthetic_$(date +%Y%m%d_%H%M%S)"
background=$(ls -dt artifacts/runs/preprocess_*/hmp2_shared_microbiome_proportions.csv | head -n 1)
mkdir -p "artifacts/runs/${run_id}" logs
PYTHONPATH=src python scripts/run_synthetic_validation.py \
  --background "$background" \
  --output-dir "artifacts/runs/${run_id}"

