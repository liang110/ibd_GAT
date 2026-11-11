#!/bin/bash
#SBATCH --job-name=ibd_capacity_control
#SBATCH --partition=normal,normal1,normal2,normal3
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=artifacts/runs/capacity_%j.out
#SBATCH --error=artifacts/runs/capacity_%j.err
set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}";conda activate ckd_env
REPO=${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)};cd "$REPO";export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"
C=$(find artifacts/runs -type f -name contract_summary.json|sort|tail -1|xargs dirname);O="artifacts/runs/capacity_$(date +%Y%m%d_%H%M%S)";mkdir -p "$O";python scripts/run_capacity_control.py --contract-dir "$C" --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} --output-dir "$O" --n-splits 10

