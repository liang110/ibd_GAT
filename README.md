# A framework for benchmarking microbiome classifiers across heterogeneous inflammatory bowel disease cohorts

This software compares classifiers on a specified target population, scoring unit and metric. In the IBD study, abundance scale and taxon coverage determine which inputs can be aligned between HMP2 and the Franzosa panel. Recruitment source and repeat visits determine how the same external predictions are scored. The framework fits models on the source cohort and compares their fixed predictions using paired target-participant draws.

## Quick start

Run these commands from the repository root with Python 3.10 or later. The example requires NumPy and pandas. Two toy classifiers each provide two sets of predictions for 12 records from eight fictional participants.

```bash
python -m pip install "numpy>=1.26" "pandas>=2.1"
python scripts/score_fixed_predictions.py --predictions examples/fixed_predictions/predictions.csv --target-manifest examples/fixed_predictions/targets.csv --candidate procedure_a --references procedure_b --expected-fits 2 --n-bootstrap 5000 --seed 73 --output toy_scores.csv
python examples/fixed_predictions/check_output.py toy_scores.csv
```

Use a fresh path for the output. The last command checks all 18 comparisons against the expected scores. For example, the pooled participant-level AUROC is 0.6406 for procedure A and 0.7500 for procedure B.

## Use the IBD data or your own predictions

The study guide covers downloading the public inputs, preparing the 126 shared candidate taxa, fitting the primary GAT and tabular references, and scoring their predictions in either transfer direction. To use your own classifier, supply its predicted probabilities and a target manifest to the same scorer used in the quick start. [Study guide and input schema](docs/reproducibility.md)

This command prints the study's ten graph attention–random forest comparisons for AUROC and Brier score across five target settings.

```bash
python scripts/show_target_decisions.py
```

Each row gives the target, scoring unit, metric, observed difference and conditional interval. The source index maps the manuscript results to their underlying data. [Figure and table sources](manuscript/registry/figure_table_source_registry_20261117.csv) · [Verification commands](repro/README.md)

## Citation

Please cite the accompanying manuscript when using this software or the study results. [Citation metadata](CITATION.cff) · [License](LICENSE)
