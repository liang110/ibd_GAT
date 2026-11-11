#!/bin/bash
#SBATCH -J ibd_p0p1
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=1-12:00:00
#SBATCH -o logs/ibd_p0p1.%A_%a.out
#SBATCH -e logs/ibd_p0p1.%A_%a.err
#SBATCH --array=0-5
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
mkdir -p logs artifacts/runs/p0_p1_20260818
case ${SLURM_ARRAY_TASK_ID} in
  0) CAP=50; PREV=0.10 ;;
  1) CAP=100; PREV=0.10 ;;
  2) CAP=126; PREV=0.10 ;;
  3) CAP=126; PREV=0.05 ;;
  4) CAP=126; PREV=0.20 ;;
  5)
    srun python -u scripts/run_learning_curves.py \
      --contract-dir artifacts/runs/contract_20260815_225635/data_contract \
      --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
      --nielsen-projection artifacts/runs/nielsen_target_20260816_083000/nielsenhb_2014_source_feature_projection.csv \
      --nielsen-manifest artifacts/runs/nielsen_target_20260816_083000/cohort_manifest_nielsenhb_2014.csv \
      --output-dir artifacts/runs/p0_p1_20260818/learning_curve \
      --n-splits 10
    srun python -u scripts/analyze_calibration_and_stability.py \
      --predictions artifacts/runs/reciprocal_nielsen_standard_289766/combined/six_direction_predictions.csv \
      --features artifacts/runs/p0_p1_20260818/learning_curve/learning_curve_features.csv \
      --weights artifacts/runs/p0_p1_20260818/learning_curve/learning_curve_weights.csv \
      --output-dir artifacts/runs/p0_p1_20260818/calibration \
      --n-bootstrap 500
    exit 0
    ;;
esac
RUN=artifacts/runs/p0_p1_20260818/contract_${CAP}_${PREV}
mkdir -p "${RUN}"
srun python -u scripts/run_contract_sensitivity.py \
  --contract-dir artifacts/runs/contract_20260815_225635/data_contract \
  --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
  --nielsen-projection artifacts/runs/nielsen_target_20260816_083000/nielsenhb_2014_source_feature_projection.csv \
  --nielsen-manifest artifacts/runs/nielsen_target_20260816_083000/cohort_manifest_nielsenhb_2014.csv \
  --output-dir "${RUN}" --feature-cap ${CAP} --min-prevalence ${PREV} --n-splits 10

