#!/bin/bash
#SBATCH --job-name=ibd_prism_origin
#SBATCH --partition=normal,normal1,normal2,normal3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=logs/13_prism_origin_%j.out
#SBATCH --error=logs/13_prism_origin_%j.err

set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
contract=$(find artifacts/runs -maxdepth 1 -type d -name 'contract_*' | sort | tail -n 1)/data_contract
run_id="prism_origin_$(date +%Y%m%d_%H%M%S)"
out="artifacts/runs/${run_id}"
mkdir -p "$out/bootstrap" logs
PYTHONPATH=src python scripts/run_prism_origin_transfer.py \
  --contract-dir "$contract" \
  --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
  --output-dir "$out" \
  --n-splits 5
PYTHONPATH=src python scripts/bootstrap_transfer_uncertainty.py \
  --predictions "$out/prism_origin_external_predictions.csv" \
  --output-dir "$out/bootstrap" \
  --n-bootstrap 500

