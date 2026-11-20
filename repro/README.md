# Check the study results

Run these commands from the repository root. The checks use aggregate data, recorded graphs and synthetic examples.

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python repro/audits/validate_current_manuscript_index.py --root .
```

The index check locates the registered sources and verifies the target-specific GAT–random forest comparisons, pooled scores by fit, the contributions to pooled AUROC, and the three participant-weighting stages for Brier error. It also compares the 45 distinct source-fit graph pairs in each direction after matching node indices to taxa. [Manuscript source index](../manuscript/registry/figure_table_source_registry_20261117.csv)

## Individual checks

```bash
python scripts/decompose_fixed_predictions.py --check
python scripts/figure5_model_comparison.py --check
python repro/audits/validate_truth_acceptance.py --ledger results/reconciled_evidence/truth_acceptance/per_seed_acceptance.csv
python repro/audits/validate_reconstructed_edge_lists.py --archive releases/formal_topology_edge_export_20261114.zip
```

The first command checks the source-pair AUROC and participant-weighted Brier components. The second verifies the absolute scores and paired AUROC contrasts used in Supplementary Table 12 and Supplementary Figure 3. The remaining commands check the known-truth decisions for individual seeds and the archived input edges for the topology experiment.

With Matplotlib installed, the plotting script renders the target-specific comparison from the three aggregate source tables.

```bash
python -m pip install matplotlib
python scripts/plot_target_decisions.py --root . --output-dir figures_out
```

The [scoring example](../README.md#quick-start) checks the score calculation with synthetic participant records. To reproduce external predictions from the published cohort inputs, follow the [study guide](../docs/reproducibility.md).
