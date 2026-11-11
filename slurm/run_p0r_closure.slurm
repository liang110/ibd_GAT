#!/bin/bash
#SBATCH -J ibd_p0r
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=18:00:00
#SBATCH -o logs/ibd_p0r.%j.out
#SBATCH -e logs/ibd_p0r.%j.err
set -eo pipefail
set +u
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
set -u
cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
export PYTHONPATH=src
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PYTHONUNBUFFERED=1

RUN=artifacts/runs/p0r_closure_20260818
mkdir -p "${RUN}"

srun python -u scripts/audit_p0r_protocol.py \
  --project-root . \
  --output-dir "${RUN}/protocol_audit"

srun python -u scripts/run_p0r_closure_analysis.py \
  --predictions artifacts/runs/reciprocal_nielsen_standard_289766/combined/six_direction_predictions.csv \
  --hmp2-manifest artifacts/runs/contract_20260815_225635/data_contract/cohort_manifest_hmp2.csv \
  --prism-manifest artifacts/runs/contract_20260815_225635/data_contract/cohort_manifest_prism.csv \
  --nielsen-manifest artifacts/runs/nielsen_target_20260816_083000/cohort_manifest_nielsenhb_2014.csv \
  --output-dir "${RUN}/analysis" \
  --n-bootstrap 500 \
  --seed 20260818 \
  --expected-features 126

srun pytest -q tests

