#!/bin/bash
#SBATCH --job-name=ibd_subgroup_audit
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/19_subgroup_audit_%j.out
#SBATCH --error=logs/19_subgroup_audit_%j.err

set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
contract=artifacts/runs/contract_20260815_225635/data_contract
run_id="subgroup_feasibility_$(date +%Y%m%d_%H%M%S)"
out="artifacts/runs/${run_id}"
mkdir -p "$out" logs
PYTHONPATH=src python scripts/audit_subgroup_feasibility.py \
  --contract-dir "$contract" \
  --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
  --output-dir "$out" \
  --min-participants-per-class 20

