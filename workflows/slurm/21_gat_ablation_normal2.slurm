#!/bin/bash
#SBATCH --job-name=ibd_gat_ablation
#SBATCH --partition=normal2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=artifacts/runs/gat_%j.out
#SBATCH --error=artifacts/runs/gat_%j.err
set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
REPO=${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
LEGACY=${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root}
cd "$REPO"; mkdir -p artifacts/runs
MODE="${ABLATION_MODE:?ABLATION_MODE must be exported}"
TAG="gat_${MODE}_$(date +%Y%m%d_%H%M%S)"
export FILE_ABD="$LEGACY/manual_abundance.csv"
export FILE_METAB="$LEGACY/metabolites.csv"
export FILE_META="$LEGACY/clean_metadata.csv"
export ABLATION_MODE="$MODE"
export EVAL_MODE=intra
export N_SPLITS=5
export SPLIT_METHOD=stratified_group
export XGB_USE_META=0
export OUTPUT_ROOT="$REPO/artifacts/runs/$TAG"
python "$LEGACY/step4_results_bundle_20260424_144620/code/train_gat_v20_ablations_journal_incremental_light.py"

