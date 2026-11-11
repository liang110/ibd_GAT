# Supplementary Information

## Supplementary Methods

The supplementary package is a provenance layer for the main manuscript. Each artifact is linked to one frozen run identifier, one source-target direction and one feature-contract hash. The package does not introduce a second estimand. It records the input contract, participant manifests, known-truth controls, origin-held-out predictions, topology exports, capacity control and source-split paired bootstrap outputs used in the main text.

### Input and contract records

Gate0 records the HMP2 and PRISM input SHA-256 hashes, canonical feature crosswalk, 126-feature order, denylist, sample and participant manifests, zero-mass decisions and the exact-name metabolite intersection. Three HMP2 samples were excluded because their total mass in the shared feature space was zero; no PRISM sample was excluded. Clinical, outcome and identifier fields were rejected before any model branch was fitted.

### Leakage and known-truth controls

The repeated-leakage package contains ten seeds for naive sample splitting and participant-grouped splitting, overlap counts and grouped-split assertions. The planted clinical-leakage package contains age, faecal-calprotectin, medication and label-proxy cases; the validator rejected all planted cases. Synthetic null, batch-only, signal-only, reversing-batch and topology-shift scenarios are known-truth stress tests. Their purpose is to test pipeline behaviour under controlled perturbation, not to validate IBD biology.

### Origin, topology and capacity packages

The PRISM origin package keeps discovery and validation labels separate and stores common target predictions for each source repeat. The topology package stores Pearson, Spearman, rewired, identity and no-edge exports for the five-repeat ablation. The capacity package stores the approximate no-message-passing comparator using the same source-only transformed input boundary. The out-degree diagnostic is reported as a construction condition: fixed k=8 gives zero variance in the out-degree vector, so conventional degree correlation is undefined; the artifact value 1.0 is not ecological evidence.

### Paired uncertainty

For each source split, bootstrap resampling was performed at the participant level and predictions were held fixed. Split-specific AUROC, average precision, Brier and paired differences are stored as source CSV files. The main manuscript reports arithmetic means of split-specific point estimates and interval endpoints where a ten-split summary is required. The supplementary tables retain the split-level records so that a future reanalysis can recompute pooled or hierarchical intervals without changing the primary estimand.

## Supplementary Tables

**Supplementary Table 1. Frozen input and evidence contract.** Raw sample and feature dimensions, shared feature count, zero-mass exclusions, exact-name metabolite intersection, hashes and denylist status.

**Supplementary Table 2. Repeated-participant and known-truth controls.** Naive and grouped AUROC, participant overlap, planted leakage rejection and synthetic scenario summaries.

**Supplementary Table 3. PRISM origin-held-out transfer.** Discovery-to-validation and validation-to-discovery model metrics, source-repeat variability and participant-cluster paired differences.

**Supplementary Table 4. Within-cohort and configuration sensitivity.** Full-Model, No-Prior and No-Context records, with the legacy whole-cohort prior explicitly excluded from the primary claim.

**Supplementary Table 5. Topology edge-set diagnostics.** Pearson–Spearman and Spearman–rewired Jaccard overlap, topology hashes and fixed-k degree diagnostic.

**Supplementary Table 6. Capacity control and source-split manifests.** Capacity-matched no-message results, model configuration, optimizer settings, split identifiers and output hashes.

## Supplementary Figure

**Supplementary Figure 1. Participant-leakage and known-truth control panels.** The figure reproduces the audited leakage, overlap, synthetic stress-test and domain-stability panels. It defines the selector’s operating boundary and is not a biological validation figure.

## Reproducibility checklist

The release package will contain the environment lockfile, Python and PyTorch versions, license, acquisition instructions for third-party public inputs, Git revision, Slurm scripts, input/output SHA-256 tables, feature-contract manifest, split manifests, model configurations, checkpoints where redistribution is permitted, figure source registry and table source registry. The GitHub repository and archival DOI will be inserted when the repository is created.

## Reference linkage

All scientific references used by the supplementary text are numbered in the main manuscript reference list. No supplementary result introduces an uncited external claim or a second set of performance values.
