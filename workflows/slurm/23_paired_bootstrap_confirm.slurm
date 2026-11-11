#!/bin/bash
#SBATCH --job-name=ibd_paired_ci
#SBATCH --partition=normal,normal1,normal2,normal3
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=artifacts/runs/paired_ci_%j.out
#SBATCH --error=artifacts/runs/paired_ci_%j.err
set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"; conda activate ckd_env
REPO=${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}; cd "$REPO"; export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"
PRED=$(find artifacts/runs -type f -name paired_predictions.csv | sort | tail -1); OUT="artifacts/runs/paired_bootstrap_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$OUT"
python scripts/paired_bootstrap_confirm.py --predictions "$PRED" --output-dir "$OUT" --n-bootstrap 1000

