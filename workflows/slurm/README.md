# Slurm launcher conventions

The launchers are intentionally thin. They should resolve one contract and one output directory, record the environment, and invoke a Python or R entry point. They should not contain hidden feature selection, label mapping or result aggregation.

Before submission:

```bash
export IBD_GAT_ROOT=/path/to/ibd_GAT
export IBD_LEGACY_ROOT=/path/to/acquired/public-input-root
export CONDA_SH=$HOME/miniforge3/etc/profile.d/conda.sh
export PYTHONPATH="$IBD_GAT_ROOT/src"
```

The normal-queue launchers use `normal,normal1,normal2,normal3` so Slurm can place the job on an available queue. For a strict queue policy, edit the `#SBATCH --partition` line in a local copy. The repository does not assume one cluster's conda path, account, wall time or memory policy.

All raw-input paths in the launchers are environment variables. A missing `IBD_LEGACY_ROOT` must fail closed rather than silently using a personal path.
