#!/bin/bash
#SBATCH -J ibd_calib
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=06:00:00
#SBATCH -o logs/ibd_calib.%A.out
#SBATCH -e logs/ibd_calib.%A.err
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
mkdir -p logs artifacts/runs/p0_p1_20260818/calibration
srun python -u scripts/analyze_calibration_and_stability.py \
  --predictions artifacts/runs/p0_p1_20260818/learning_curve/learning_curve_predictions.csv \
  --features artifacts/runs/p0_p1_20260818/learning_curve/learning_curve_features.csv \
  --weights artifacts/runs/p0_p1_20260818/learning_curve/learning_curve_weights.csv \
  --metrics artifacts/runs/p0_p1_20260818/learning_curve/learning_curve_metrics.csv \
  --output-dir artifacts/runs/p0_p1_20260818/calibration \
  --n-bootstrap 500

