#!/bin/bash
# Post-freeze reporting sensitivity: consumes only the immutable six-edge predictions.
# It neither refits models nor accesses abundance matrices or target labels for fitting.
#SBATCH -J nielsen_target_participant
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -o logs/nielsen_target_participant.%j.out
#SBATCH -e logs/nielsen_target_participant.%j.err

set -euo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export PYTHONUNBUFFERED=1

cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
run_id="nielsen_target_participant_unit_20260821_${SLURM_JOB_ID}"
outdir="artifacts/runs/${run_id}"
mkdir -p logs "${outdir}"

sha256sum artifacts/runs/reciprocal_nielsen_standard_289766/combined/six_direction_predictions.csv \
  > "${outdir}/input_predictions.sha256"

python -u scripts/analyze_participant_unit_sensitivity_nielsen_target_20260821.py \
  --predictions artifacts/runs/reciprocal_nielsen_standard_289766/combined/six_direction_predictions.csv \
  --outdir "${outdir}" \
  --bootstrap 500 \
  --seed 20260819 \
  --direction HMP2_to_NielsenHB_2014 \
  --direction PRISM_to_NielsenHB_2014

