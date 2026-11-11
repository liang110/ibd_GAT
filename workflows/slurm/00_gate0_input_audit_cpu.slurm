#!/bin/bash
#SBATCH --job-name=ibd_gate0
#SBATCH --partition=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=artifacts/runs/gate0_%j.out
#SBATCH --error=artifacts/runs/gate0_%j.err
set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
REPO=${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
cd "$REPO"; export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"
RUN_ID="gate0_$(date +%Y%m%d_%H%M%S)"
OUT="artifacts/runs/$RUN_ID"; mkdir -p "$OUT"
python scripts/audit_old_ibd_inputs.py \
  --hmp2-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
  --prism-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
  --output-dir "$OUT"

