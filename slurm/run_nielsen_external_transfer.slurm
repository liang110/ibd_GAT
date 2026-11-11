#!/bin/bash
#SBATCH -J nielsen_transfer
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH -o logs/nielsen_transfer.%j.out
#SBATCH -e logs/nielsen_transfer.%j.err

set -eo pipefail

set +u
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
set -u
export PYTHONPATH=src
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PYTHONUNBUFFERED=1

cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
run_id="nielsen_transfer_${SLURM_JOB_ID}"
mkdir -p logs "artifacts/runs/${run_id}"

srun python -u scripts/run_nielsen_external_transfer.py \
  --contract-dir artifacts/runs/contract_20260815_225635/data_contract \
  --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
  --nielsen-projection artifacts/runs/nielsen_target_20260816_083000/nielsenhb_2014_source_feature_projection.csv \
  --nielsen-manifest artifacts/runs/nielsen_target_20260816_083000/cohort_manifest_nielsenhb_2014.csv \
  --output-dir "artifacts/runs/${run_id}" \
  --n-splits 10

