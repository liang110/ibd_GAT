#!/bin/bash
#SBATCH --job-name=ibd_subtype_transfer
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=logs/20_subtype_transfer_%j.out
#SBATCH --error=logs/20_subtype_transfer_%j.err

set -eo pipefail
source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
conda activate ckd_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd ${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}
contract=artifacts/runs/contract_20260815_225635/data_contract
run_id="subtype_transfer_$(date +%Y%m%d_%H%M%S)"
out="artifacts/runs/${run_id}"
mkdir -p "$out" logs
for subtype in CD UC; do
  subtype_out="$out/$subtype"
  mkdir -p "$subtype_out/bootstrap"
  PYTHONPATH=src python scripts/run_subtype_transfer.py \
    --contract-dir "$contract" \
    --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
    --output-dir "$subtype_out" \
    --case-subtype "$subtype" \
    --n-splits 5
  PYTHONPATH=src python scripts/bootstrap_transfer_uncertainty.py \
    --predictions "$subtype_out/subtype_transfer_predictions.csv" \
    --output-dir "$subtype_out/bootstrap" \
    --n-bootstrap 500
done

