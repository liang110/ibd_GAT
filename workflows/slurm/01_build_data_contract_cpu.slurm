#!/bin/bash
#SBATCH --job-name=ibd_contract
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=artifacts/runs/contract_%j.out
#SBATCH --error=artifacts/runs/contract_%j.err

# The HPC conda activation hook references optional shell variables.
set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

REPO=${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
cd "$REPO"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"
RUN_ID="contract_$(date +%Y%m%d_%H%M%S)"
OUT="artifacts/runs/$RUN_ID/data_contract"
mkdir -p "$OUT"
python scripts/build_data_contract.py \
  --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
  --output-dir "$OUT" | tee "$OUT/console_summary.json"

