#!/bin/bash
#SBATCH -J ibd_nielsen_reciprocal
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH -o logs/ibd_nielsen_reciprocal.%j.out
#SBATCH -e logs/ibd_nielsen_reciprocal.%j.err

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
run_id="reciprocal_nielsen_standard_${SLURM_JOB_ID}"
mkdir -p logs "artifacts/runs/${run_id}" "artifacts/runs/${run_id}/combined"

srun python -u scripts/run_reciprocal_nielsen_standard_benchmark.py \
  --contract-dir artifacts/runs/contract_20260815_225635/data_contract \
  --legacy-root ${IBD_LEGACY_ROOT:?Set IBD_LEGACY_ROOT to the public-input root} \
  --nielsen-projection artifacts/runs/nielsen_target_20260816_083000/nielsenhb_2014_source_feature_projection.csv \
  --nielsen-manifest artifacts/runs/nielsen_target_20260816_083000/cohort_manifest_nielsenhb_2014.csv \
  --output-dir "artifacts/runs/${run_id}" \
  --n-splits 10

# Merge the immutable four-direction result with the two newly evaluated edges.
python - "${run_id}" <<'PY'
import shutil
import sys
from pathlib import Path
import pandas as pd

run_id = sys.argv[1]
root = Path("${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}")
out = root / "artifacts/runs" / run_id
old = root / "artifacts/runs/standard_model_benchmark_286983"
combined = out / "combined"
metrics = pd.concat([
    pd.read_csv(old / "standard_model_benchmark_metrics.csv"),
    pd.read_csv(out / "reciprocal_nielsen_metrics.csv"),
], ignore_index=True)
predictions = pd.concat([
    pd.read_csv(old / "standard_model_benchmark_predictions.csv"),
    pd.read_csv(out / "reciprocal_nielsen_predictions.csv"),
], ignore_index=True)
metrics.to_csv(combined / "six_direction_metrics.csv", index=False)
predictions.to_csv(combined / "six_direction_predictions.csv", index=False)
PY

srun python -u scripts/bootstrap_transfer_uncertainty.py \
  --predictions "artifacts/runs/${run_id}/combined/six_direction_predictions.csv" \
  --output-dir "artifacts/runs/${run_id}/combined" \
  --n-bootstrap 500 \
  --seed 20260816

srun python -u scripts/analyze_six_direction_interactions.py \
  --metrics "artifacts/runs/${run_id}/combined/six_direction_metrics.csv" \
  --bootstrap-summary "artifacts/runs/${run_id}/combined/aggregate_cluster_bootstrap.csv" \
  --output-dir "artifacts/runs/${run_id}/combined"

python - "${run_id}" <<'PY'
import json
import sys
from pathlib import Path
run_id = sys.argv[1]
root = Path("${IBD_GAT_ROOT:-$(git rev-parse --show-toplevel)}")
out = root / "artifacts/runs" / run_id
summary = json.loads((out / "reciprocal_nielsen_summary.json").read_text())
summary["combined_existing_run"] = "standard_model_benchmark_286983"
summary["combined_analysis_dir"] = str(out / "combined")
(out / "run_manifest.json").write_text(json.dumps(summary, indent=2))
PY

