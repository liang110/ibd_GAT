#!/bin/bash
#SBATCH --job-name=ibd_origin_paired
#SBATCH --partition=normal,normal1,normal2,normal3
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=artifacts/runs/origin_paired_%j.out
#SBATCH --error=artifacts/runs/origin_paired_%j.err
set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}";conda activate ckd_env
R=${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)};cd "$R";export PYTHONPATH="$R/src:$R:${PYTHONPATH:-}";C=$(find artifacts/runs -type f -name contract_summary.json|sort|tail -1|xargs dirname);O="artifacts/runs/origin_paired_$(date +%Y%m%d_%H%M%S)";mkdir -p "$O";python scripts/run_origin_paired.py --contract-dir "$C" --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} --output-dir "$O" --n-splits 10 --n-bootstrap 1000

