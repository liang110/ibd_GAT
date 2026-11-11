#!/bin/bash
#SBATCH --job-name=ibd_origin_gat
#SBATCH --partition=normal,normal1,normal2,normal3
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=artifacts/runs/origin_gat_%j.out
#SBATCH --error=artifacts/runs/origin_gat_%j.err
set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"; conda activate ckd_env
REPO=${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}; cd "$REPO"; export PYTHONPATH="$REPO/src:$REPO:${PYTHONPATH:-}"
CONTRACT=$(find artifacts/runs -type f -name contract_summary.json | sort | tail -1 | xargs dirname); OUT="artifacts/runs/origin_gat_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$OUT"
python scripts/run_origin_gat_transfer.py --contract-dir "$CONTRACT" --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} --output-dir "$OUT" --n-splits "${N_SPLITS:-5}"

