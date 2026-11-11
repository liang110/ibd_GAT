#!/bin/bash
#SBATCH --job-name=ibd_topology_audit
#SBATCH --partition=normal,normal1,normal2,normal3
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=artifacts/runs/topology_audit_%j.out
#SBATCH --error=artifacts/runs/topology_audit_%j.err
set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"; conda activate ckd_env
REPO=${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}; cd "$REPO"; OUT="artifacts/runs/topology_audit_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$OUT"
python scripts/audit_topology_runs.py --runs-root artifacts/runs --output-dir "$OUT"

