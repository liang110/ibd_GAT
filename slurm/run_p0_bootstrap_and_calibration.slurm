#!/bin/bash
#SBATCH -J ibd_p0boot
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=18:00:00
#SBATCH -o logs/ibd_p0boot.%A_%a.out
#SBATCH -e logs/ibd_p0boot.%A_%a.err
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
    srun python -u scripts/analyze_calibration_and_stability.py \
      --predictions artifacts/runs/p0_p1_20260818/learning_curve/learning_curve_predictions.csv \
      --features artifacts/runs/p0_p1_20260818/learning_curve/learning_curve_features.csv \
      --weights artifacts/runs/p0_p1_20260818/learning_curve/learning_curve_weights.csv \
      --metrics artifacts/runs/p0_p1_20260818/learning_curve/learning_curve_metrics.csv \
      --output-dir artifacts/runs/p0_p1_20260818/calibration \
      --n-bootstrap 500
    exit 0
    ;;
esac
RUN=artifacts/runs/p0_p1_20260818/contract_${CAP}_${PREV}
srun python -u scripts/bootstrap_grouped_predictions.py \
  --predictions "${RUN}/contract_sensitivity_predictions.csv" \
  --output-dir "${RUN}/bootstrap" \
  --group-columns variant_id,direction,split,model \
  --aggregate-columns variant_id,direction,model \
  --n-bootstrap 500 --seed 20260818

