# Verify the study results

Run the following commands from the repository root with the Python dependencies installed.

```bash
python repro/audits/validate_current_manuscript_index.py --root .
python -m pytest -q
```

The [results index](../results/display_sources.csv) identifies the data for all 28 manuscript and supplementary displays. The first command checks source-fit scores, shared target-participant intervals, the ten GAT–random forest comparisons, the AUROC and Brier decompositions, graph inputs and known-truth controls.

For a focused check:

```bash
python repro/audits/validate_truth_acceptance.py --root .
python repro/audits/validate_topology.py --root .
```

The truth check applies the published acceptance rules to both sets of fifty seed-level results. The topology check reads each graph's directed edge list, split-specific taxon order and recorded fit metrics. It verifies the five degree-preserving rewires against each associated Spearman graph, checks self-loop controls and recalculates the 45 between-fit overlaps per direction.
