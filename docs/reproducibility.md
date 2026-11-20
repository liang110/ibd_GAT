# Study inputs, training and external scoring

Run commands from the repository root in an activated Python 3.10+ environment. Install the package with `python -m pip install -e ".[input]"` to read the provider files. The core dependency versions used for local validation are listed in [requirements-reproducibility.txt](../requirements-reproducibility.txt).

## Prepare the public inputs

Download these files to a local directory named `provider_files`. The release catalogue records the expected filenames and fingerprints. [Input catalogue](../configs/data_catalog.json)

| Filename | Provider and contents |
| --- | --- |
| `taxonomic_profiles_3.tsv.gz` | HMP2 merged metagenomic taxonomic profiles. [IBDMDB release, 4 May 2018](https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/2018-05-04/taxonomic_profiles_3.tsv.gz) |
| `hmp2_metadata_2018-08-20.csv` | Assay, sample, participant, visit, diagnosis and read depth. [IBDMDB metadata, 20 August 2018](https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/metadata/hmp2_metadata_2018-08-20.csv) |
| `41564_2018_306_MOESM6_ESM.xlsx` | Species abundances and diagnosis for PRISM discovery and Netherlands validation, sheet `dataset_s4`. [Franzosa et al., Supplementary Dataset 4](https://www.nature.com/articles/s41564-018-0306-4#Sec35) |
| `dysbiosis_scores.tsv` (optional) | Included sample IDs for the depth sensitivity analysis. [Published HMP2 scores](https://forum.biobakery.org/uploads/short-url/umwfR0kDJ6s5RXHwtIMgLaOEoOI.tsv) |

```bash
python scripts/prepare_published_profiles.py --input-dir provider_files --output-dir prepared
```

Choose a new output directory. The importer writes the sample manifests and abundance matrices and checks their sample-keyed values against the ordered 126-taxon panel. HMP2 contributes 1,626 eligible profiles from 130 participants after exclusion of 12 profiles with zero mass over the shared taxa. The Franzosa panel contributes one profile per person from 155 PRISM discovery participants and 65 Netherlands validation participants. To verify the 1,595 published depth-qualified IDs, add `--published-dysbiosis-ids provider_files/dysbiosis_scores.tsv`.

## Fit the source models

Install PyTorch and PyTorch Geometric versions compatible with your CPU or CUDA environment. The recorded local stack used torch 2.7.1+cu118 and torch-geometric 2.4.0.

```bash
python -m pip install -e ".[model]"
python scripts/run_no_prior_gat_transfer.py --contract-dir prepared/contract --legacy-root prepared --output-dir run_primary --topology spearman --n-splits 10
```

This entry point fits GAT, elastic net and random forest in both directions. Each source split groups records by participant. The primary GAT uses minibatches of 32 and a source-derived Spearman graph. The separate model comparisons in the paper include a no-message MLP and XGBoost. XGBoost uses source-grouped parameter selection and then refits on the complete source panel. The independent topology experiment uses full-batch updates and a degree-preserving rewiring procedure. Its formal record identifies the exact input graphs. [Training configurations](../repro/audits/training_configuration.json) · [Topology run manifest](../results/reconciled_evidence/topology_constraints/formal_topology_runtime_manifest.json)

## Score fixed predictions

```bash
python scripts/score_fixed_predictions.py --predictions run_primary/paired_predictions.csv --prepared-dir prepared --candidate no_prior_gat --references random_forest elastic_net --expected-fits 10 --n-bootstrap 5000 --output fixed_prediction_scores.csv
```

For an additional classifier, supply its CSV after `--predictions` alongside the other model files. Use one row per direction, fit, classifier and target record.

| Input | Required columns |
| --- | --- |
| Predictions | `direction,split,model,sample_id,participant_id,outcome,probability` |
| Target manifest for `--target-manifest`, an alternative to `--prepared-dir` | `direction,sample_id,participant_id,outcome,origin` |

Code `outcome` as 0 or 1 and use `probability` for outcome 1. All compared procedures must cover the same target records in each fit. For a target without recruitment groups leave `origin` blank throughout. Outputs include AUROC, average precision, Brier score and paired differences. A positive AUROC or average precision difference, or a negative Brier difference, favours the candidate. Each conditional interval uses shared participant draws across models and fixed source fits in the specified setting.

Record-level scoring weights each profile equally. Participant-level scoring averages the probabilities for each participant before calculating the metric. The following command assesses weighting within HMP2 collections and restriction to the published high-depth panel.

```bash
python scripts/score_hmp2_collection_sensitivity.py --predictions run_primary/paired_predictions.csv --metadata provider_files/hmp2_metadata_2018-08-20.csv --published-scores provider_files/dysbiosis_scores.tsv --output hmp2_sensitivity.csv
```

The manuscript index identifies the result sources, including the five-procedure pooled scores, the separate SIAMCAT reference and the HMP2 collection sensitivity analysis. [Study results](../manuscript/registry/figure_table_source_registry_20261117.csv) · [Checks on released results](../repro/README.md)
