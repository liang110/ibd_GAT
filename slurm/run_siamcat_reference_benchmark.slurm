#!/bin/bash
#SBATCH -J ibd_siamcat_reference
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH -o logs/ibd_siamcat_reference.%j.out
#SBATCH -e logs/ibd_siamcat_reference.%j.err

set -eo pipefail
set +u
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
set -u
export PYTHONPATH=src
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PYTHONUNBUFFERED=1

cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
run_id="siamcat_reference_${SLURM_JOB_ID}"
run_dir="artifacts/runs/${run_id}"
input_dir="${run_dir}/inputs"
mkdir -p logs "${input_dir}"

python -u scripts/export_siamcat_inputs.py \
  --contract-dir artifacts/runs/contract_20260815_225635/data_contract \
  --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
  --nielsen-projection artifacts/runs/nielsen_target_20260816_083000/nielsenhb_2014_source_feature_projection.csv \
  --nielsen-manifest artifacts/runs/nielsen_target_20260816_083000/cohort_manifest_nielsenhb_2014.csv \
  --output-dir "${input_dir}"

python -u scripts/run_siamcat_reference_benchmark.py \
  --inputs-dir "${input_dir}" \
  --output-dir "${run_dir}" \
  --config configs/experiments/siamcat_reference_v1.yaml \
  --reference-script scripts/run_siamcat_reference.R \
  --rscript ${SIAMCAT_RSCRIPT:-Rscript} \
  --n-splits 10

python -u scripts/bootstrap_transfer_uncertainty.py \
  --predictions "${run_dir}/siamcat_reference_predictions.csv" \
  --output-dir "${run_dir}" \
  --n-bootstrap 500 \
  --seed 20260816

