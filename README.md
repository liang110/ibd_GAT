# A framework for benchmarking microbiome classifiers across heterogeneous inflammatory bowel disease cohorts

This repository compares classifiers for IBD versus non-IBD discrimination across HMP2 and the Franzosa panel. Abundance scale and reported taxa define comparable inputs. Repeated sampling determines participant grouping and scoring units. Recruitment defines the populations in which the same external predictions are evaluated.

The [framework diagram](manuscript/figures/Figure1_population_benchmark.png) is the Figure 1 embedded in the accompanying manuscript. The section drafts in `manuscript/clean/` and the older plotting files document earlier editorial versions. The [display index](manuscript/registry/figure_table_source_registry_20261117.csv) identifies the current results.

## Read the model comparison

[Figure 5 and Supplementary Table 12](manuscript/source_data/figure5_model_comparison.csv) compare five fitted procedures in four settings. Rows give the ten-fit observed means and source-fit SDs for AUROC, average precision and Brier score. The paired columns give the observed mean GAT-minus-reference AUROC difference, its bootstrap mean (the Figure 5 marker) and the percentile interval from 5,000 target-participant draws shared across the ten fixed fits. The GAT row has no self-contrast. [Ten paired GAT/RF fits per setting](manuscript/source_data/figure5_gat_rf_by_fit.csv) show variation over source fits.

For example, GAT and random forest trained on HMP2 have observed mean AUROCs of 0.8318 and 0.7616 in the 65-person Netherlands validation group. Their observed paired difference is +0.0702. The conditional 95% interval, −0.0289 to +0.1736, includes a lead for either procedure in that group. The same fits scored on the 155-person PRISM discovery group give +0.1030 [0.0487, 0.1590]. With HMP2 as target, one score per person is the average of all their available record probabilities.

```bash
python scripts/figure5_model_comparison.py --check
python scripts/figure5_model_comparison.py --output table_s12.csv
```

These commands use the published aggregate tables and Python's standard library. Use a new filename for the export. To find any main or supplementary display, see the [manuscript source index](manuscript/registry/figure_table_source_registry_20261117.csv). [Supplementary Table 11](manuscript/source_data/siamcat_observed_reference.csv) reports the observed scores of a separate native SIAMCAT LASSO reference, with 385 HMP2 records from 105 people.

## Prepare the study inputs

Download these provider files into a single directory, here called `provider_files`. The [catalogue](configs/data_catalog.json) gives the exact release filenames and file fingerprints. The new input directory and its derived tables stay local.

| Source | File | Contents |
| --- | --- | --- |
| [IBDMDB MGX release, 4 May 2018](https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/2018-05-04/taxonomic_profiles_3.tsv.gz) | `taxonomic_profiles_3.tsv.gz` | Merged abundance profiles. This is the file with the `_3` suffix. |
| [IBDMDB metadata, 20 August 2018](https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/metadata/hmp2_metadata_2018-08-20.csv) | `hmp2_metadata_2018-08-20.csv` | Assay type, external sample ID, participant, visit, diagnosis, site and read depth. |
| [Franzosa et al., Supplementary Dataset 4](https://www.nature.com/articles/s41564-018-0306-4#Sec35) | `41564_2018_306_MOESM6_ESM.xlsx` | Sheet `dataset_s4`, header in row 2. Sample names carry recruitment group. The Diagnosis row supplies CD, UC or Control labels. |
| [Published HMP2 dysbiosis scores](https://forum.biobakery.org/uploads/short-url/umwfR0kDJ6s5RXHwtIMgLaOEoOI.tsv) (optional) | `dysbiosis_scores.tsv` | Provider's headerless three-column file. The importer compares its first-column IDs with the read-depth-qualified subset. |

```bash
python -m venv .venv
# Activate the virtual environment before installing the package.
python -m pip install -e ".[input,dev]"
python scripts/prepare_published_profiles.py --input-dir provider_files --output-dir prepared
python repro/audits/validate_current_manuscript_index.py --root .
```

To verify the 1,595 published sample IDs, add `--published-dysbiosis-ids provider_files/dysbiosis_scores.tsv` to the preparation command. The command creates four local CSV inputs and a `prepared/contract/` directory, checks the ordered 126-taxon registry, and verifies the sample-keyed projected abundances against the study panels. It requires a fresh output path and reports an error if the release fails these checks. For HMP2, species rows are matched to metagenomic External IDs after removing the profile suffix. For Franzosa, the loader excludes eight clinical/identifier rows before forming the taxon intersection. All fitted transformations and feature selection subsequently use source data.

### Which observations entered the benchmark?

HMP2 supplies 1,638 profiles matched to metadata. Twelve have zero mass over the 126 shared taxa. The remaining 1,626 profiles represent 130 participants (65 CD, 38 UC, 27 controls) and 1,584 participant–visit combinations. In the [original HMP2 analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC6650278/), downstream MGX analyses used a one-million-read threshold after host filtering and required detected microbial abundance. All 1,595 scored sample IDs match included profiles here. The 31 additional included profiles from 18 participants fall below the read-depth threshold. All 12 zero-shared-mass exclusions also fall below it. The model results use the full 1,626-record panel.

Forty included participant–visit combinations have multiple profiles linked to the same collection. Of the 42 additional profiles at these keys, 20 carry technical-replicate and 22 pilot profile designations. Participant grouping keeps all such profiles together on the source side. On the target side, record-level metrics weight each profile equally, whereas participant-level metrics average a person's probabilities before scoring. The archived results describe these specified weights.

The Franzosa panel has one profile per participant. PRISM discovery has 68 CD, 53 UC and 34 controls. Netherlands validation has 20 CD, 23 UC and 22 controls. Dutch cases were recruited through UMCG and controls through LifeLines-DEEP. Local IDs show no exact or normalized matches between studies. Study sample linkage and recruitment information do not establish whether any person enrolled in both studies.

## Run a classifier or extend the comparison

With compatible PyTorch and PyTorch Geometric installations, the primary transfer entry point accepts the prepared inputs directly. The core dependency versions are listed in [requirements-reproducibility.txt](requirements-reproducibility.txt).

```bash
python -m pip install -e ".[model]"
python scripts/run_no_prior_gat_transfer.py --contract-dir prepared/contract --legacy-root prepared --output-dir run_primary --topology spearman --n-splits 10
```

This entry point fits GAT, elastic net and random forest on grouped source partitions and writes their metrics and probabilities. The no-message MLP and XGBoost are separate procedures. XGBoost selects parameters in grouped inner source folds, then refits on the complete source panel. [Supplementary model configurations](repro/audits/training_configuration.json) describe the primary and full-batch learners. The [formal graph intervention](results/reconciled_evidence/topology_constraints/formal_topology_runtime_manifest.json) uses an archived degree-preserving rewiring procedure and compares pre-convolution inputs. The generic launcher's shuffled-graph option uses a different randomization.

For an additional classifier, generate one out-of-cohort probability for each source-fit ID and target record; preserve the shared target record, participant, outcome and recruitment-group keys. Evaluate the same target records under each paired procedure. The source-fitting and target-sampling implementations are in [evaluation](src/ibd_microbiome_framework/evaluation/) and the [common-target analysis script](scripts/run_common_target_bootstrap_reanalysis.py). For a new evaluation population, record its measured taxa, outcome, participant IDs and recruitment definitions before computing these comparisons. Source-fit SD and target-participant uncertainty describe different sources of variability.

Compact aggregate results are in [the paired summaries](results/p0_rechecked_20261117/) and [the absolute performance tables](results/reconciled_evidence/absolute_performance_summary.csv). Older dated artifacts retain their original conventions. The current display index selects the registered estimates for the manuscript, including the 5,000 shared draws for direct recruitment-group transfer.
