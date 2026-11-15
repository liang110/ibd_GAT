# Supplementary Information

## Supplementary index

Supplementary Table 1 records the frozen input and feature contract; Table 2 records repeated-participant and known-truth controls; Table 3 records PRISM origin transfer; Table 4 records within-cohort and configuration sensitivity; Table 5 records contract-aligned topology stability; Table 6 records the capacity comparator and paired bootstrap manifest. Supplementary Figure 1 records selector operating characteristics under the separately versioned `v2_prevalent_truth_features` benchmark. Machine-readable source CSV files, immutable manifests and the figure/table source registry are stored beside the manuscript source.

## Supplementary Methods

The supplementary package is a provenance layer for the main manuscript. Each artifact is linked to one experiment identifier, one source-target direction, one source split and one frozen feature-contract definition. The package does not introduce a second primary estimand.

The input-contract manifest records the HMP2 and PRISM feature registries, canonical crosswalk, 126-feature order, denylist, sample and participant manifests, zero-mass decisions and the exact-name metabolite intersection. Twelve HMP2 samples were excluded because their total mass in the shared feature space was zero; no PRISM sample was excluded. A separate manifest audit found zero exact, normalized and substring overlap for sample or participant identifiers between HMP2 and PRISM in both the manifest and included panels.

The repeated-leakage package contains ten seeds for sample-level and participant-grouped splitting, overlap counts and grouped-split assertions. The current contract uses 80 synthetic participants, three observations per participant, 240 observations in total, and 180/60 train/test observations per policy; grouped splitting uses 60/20 train/test participants. The planted clinical-field package contains age, faecal-calprotectin, medication and label-proxy cases; the validator rejected all 10/10 cases. The known-truth benchmark is versioned as `v2_prevalent_truth_features`; null and batch-only controls remained near chance, while signal, reversing-batch and topology-shift scenarios recovered planted features under their declared acceptance gates.

The PRISM origin package keeps discovery and validation labels separate and stores predictions for each source split. The final topology package stores Pearson, Spearman and destination-shuffled edge exports for ten repeats in each direction, together with edge counts and Jaccard summaries. Fixed k=8 produces zero out-degree variance, so degree correlation is NA under the prespecified policy. Exported Spearman edge records match the transfer metrics for all 60 rows checked.

For each source split, participant-cluster bootstrap resamples target participants and applies the same target indices to all models in a paired comparison. Baseline GAT-versus-elastic-net/random-forest comparisons use B=5,000 with zero invalid draws; the capacity-matched GAT-versus-MLP comparison also uses B=5,000 with zero invalid draws. Ten-split summaries average split-specific point estimates and interval endpoints; they are not pooled intervals.

## Supplementary table descriptions

**Supplementary Table 1. Frozen input and evidence contract.** Manifest and included sample counts, participant counts, shared feature count, zero-mass exclusions, exact-name metabolite intersection and independence-audit status.

**Supplementary Table 2. Repeated-participant and known-truth controls.** Naive and grouped AUROC, participant overlap, planted-field rejection, five synthetic scenarios, recovery counts, acceptance rates and denominators.

**Supplementary Table 3. PRISM origin-held-out transfer.** Ten-split No-Prior/elastic-net/random-forest metrics and the separate tabular origin audit.

**Supplementary Table 4. Within-cohort and configuration sensitivity.** Full-Model, No-Prior, tabular and capacity-control records, with the legacy whole-cohort prior excluded from the primary transfer estimand.

**Supplementary Table 5. Contract-aligned topology stability.** Ten-split edge counts, Jaccard summaries and the fixed-k degree diagnostic.

**Supplementary Table 6. Capacity control and paired uncertainty.** Model IDs, split IDs, source and target counts, parameter budget, point metrics and B=5,000 paired bootstrap intervals.

## Supplementary Figure 1

**Supplementary Figure 1. Selector operating characteristics.** The figure uses the separately versioned local selector audit `v2_prevalent_truth_features` with ten seeds. It reports disease-feature recall, batch false-positive rate and the expected disease true-positive count under a random selector with the same budget. These operating characteristics are diagnostic controls for the selector and are not biological validation or a second transfer estimand.

## Reproducibility checklist

The public release package contains the environment specification, production code, compact source tables, split manifests, model configurations, audit summaries, figure source registry and table source registry. Full prediction ledgers, participant identifiers, checkpoints and restricted raw data are not redistributed.

## Reference linkage

Supplementary scientific statements use the numbered reference list in the main manuscript. No supplementary result introduces a second set of transfer performance values.


