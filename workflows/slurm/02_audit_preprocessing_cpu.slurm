#!/bin/bash
#SBATCH --job-name=ibd_preprocess_audit
#SBATCH --partition=normal,normal1,normal2,normal3
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=artifacts/runs/preprocess_%j.out
#SBATCH --error=artifacts/runs/preprocess_%j.err

set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
REPO=${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
cd "$REPO"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"
CONTRACT=$(find artifacts/runs -maxdepth 1 -type d -name 'contract_*' | sort | tail -n 1)/data_contract
RUN_ID="preprocess_$(date +%Y%m%d_%H%M%S)"
OUT="artifacts/runs/$RUN_ID"
mkdir -p "$OUT"
python scripts/audit_preprocessing.py \
  --contract-dir "$CONTRACT" \
  --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
  --output-dir "$OUT" \
  --n-splits 5 | tee "$OUT/console_summary.json"

