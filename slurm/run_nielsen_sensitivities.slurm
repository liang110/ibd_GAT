#!/bin/bash
#SBATCH -J nielsen_sensitivity
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH -o logs/nielsen_sensitivity.%j.out
#SBATCH -e logs/nielsen_sensitivity.%j.err

set -eo pipefail
set +u
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
set -u
export PYTHONPATH=src
export PYTHONUNBUFFERED=1

cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
run_id="nielsen_sensitivity_${SLURM_JOB_ID}"
mkdir -p logs "artifacts/runs/${run_id}"

srun python -u scripts/audit_nielsen_sensitivities.py \
  --predictions artifacts/runs/nielsen_transfer_286472/nielsen_external_transfer_predictions.csv \
  --manifest artifacts/runs/nielsen_target_20260816_083000/cohort_manifest_nielsenhb_2014.csv \
  --output-dir "artifacts/runs/${run_id}" \
  --n-bootstrap 500 \
  --seed 20260816

