# Architecture

## Design objective

The package separates data identity, fitted transforms, model fitting,
evaluation, and reporting. This prevents the legacy failure mode in which
cleaning, feature selection, and model-specific code diverged between the main
model and baselines.

## Dependency flow

```text
data/catalog -> data/ingest -> preprocessing/canonicalize -> harmonization/shared_space
                                                           -> preprocessing/fitted_transform
shared feature matrix + split plan -> models/* -> evaluation/* -> statistics/* -> reporting/*
```

`data` owns immutable input metadata and identifiers. `preprocessing` owns
training-fit transforms. `harmonization` owns taxonomy/metabolite crosswalks and
the persisted shared-feature contract. `models` must accept an already prepared
matrix. `evaluation` owns split plans and leakage assertions. `statistics`
computes uncertainty. `reporting` reads immutable result tables only.

## Required contracts

1. `CohortManifest`: sample ID, participant ID, outcome, diagnosis subtype,
   collection origin, batch/site, timepoint, and source checksum.
2. `FeatureRegistry`: canonical feature ID, source feature ID, modality,
   mapping confidence, and inclusion reason.
3. `SplitPlan`: training/validation/test sample IDs, grouping column, seed,
   repeated split ID, and provenance hash.
4. `TransformState`: training-only normalization, zero-replacement, CLR/scaler,
   feature filter, and any training-only prior.
5. `RunManifest`: Git revision, config hash, environment lock, input checksums,
   feature counts, class counts, timings, and output paths.

## Experiment classes

- `primary_external`: repeated source-group fits on a fixed untouched target.
- `within_cohort_audit`: participant/origin-aware confound sensitivity analyses.
- `synthetic_validation`: known disease, confound, leakage, and topology truth.
- `exploratory_graph`: graph stability and biological hypothesis generation.
- `target_adaptation`: labelled target data, never mixed with source-only claims.

## Drug-synergy extension

Drug synergy is a downstream, separate task. It may consume validated,
cohort-stable host-microbe-metabolite features but must use its own drug-pair
response matrix, treatment-context split groups, synergy metric (for example
Bliss/Loewe/HSA), and external validation. It must not infer synergy from the
current binary IBD classification labels.
