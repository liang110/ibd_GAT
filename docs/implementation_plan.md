# Executable Implementation Plan

## Operating policy

Work only in this repository. `legacy/` and `artifacts/legacy_*` are evidence,
not executable inputs. Every milestone is developed on a short-lived branch,
reviewed against its acceptance gate, then merged into `main`. No new result is
used in the manuscript until its run manifest and source table are committed or
archived under a content-addressed run ID.

## Milestone 0: governance and reproducibility

**Input:** frozen legacy evidence and the experimental ledger.

**Implement:**

1. Add `CODEOWNERS`, contributor instructions, a change log, and a decision
   record template.
2. Add a lockable environment (`conda-lock` or pip lock), formatter/linter,
   pytest, and a small CI workflow if the repository gains a remote host.
3. Create a data catalog schema and a run-manifest schema.

**Accept only when:** `pytest`, linting, and manifest-schema validation run on a
clean checkout; a run records Git revision, config hash, checksums, seeds,
class counts, group counts, and feature counts.

## Milestone 1: data contracts and shared feature space

**Input:** original HMP2 and PRISM matrices plus source publications/tables.

**Implement:**

1. Retain PRISM discovery/validation origin, geography, technical batch, and
   treatment/clinical covariates in `CohortManifest`.
2. Convert microbial identifiers to one versioned NCBI/GTDB taxonomy registry;
   record one-to-many and unresolved mappings rather than guessing.
3. Resolve metabolites to HMDB/KEGG/InChIKey with mapping confidence and
   provenance. Create a separate microbiome-only contract when a defensible
   shared metabolite universe is unavailable.
4. Persist the exact shared feature set per modality and enforce it at model
   input. Absent target features cause a contract error, not silent zero fill.

**Accept only when:** mapping reports list all input, mapped, excluded, and
shared features; sample-library-sum diagnostics are stored; PRISM origin is
non-null; each external feature matrix has identical ordered feature IDs.

## Milestone 2: leak-proof preprocessing and splits

**Implement:**

1. Put filtering, zero replacement, compositional transform (CLR or justified
   alternative), scaling, MI selection, and clinical-variable removal in one
   fitted `Pipeline` interface used by all models.
2. Create `SplitPlan` generators for HMP2 participant grouping, PRISM origin
   grouping, leave-one-origin/study-out, and repeated source-to-fixed-target
   transfer.
3. Add automated tests that plant `Age`, `Fecal.Calprotectin`, outcome labels,
   and target-only features; tests must fail when any reaches an omics model.

**Accept only when:** fitted state is created only with training IDs; all
comparators receive byte-identical prepared matrices; the leakage tests pass.

## Milestone 3: known-truth method validation

**Implement:**

1. Generate semi-synthetic cohorts preserving real compositional covariance
   while independently injecting disease, site/batch, repeated-subject,
   clinical-leakage, and topology-shift signals.
2. Pre-register recovery metrics: disease-effect ranking recovery, confounder
   detection sensitivity/specificity, calibration, and false portable-signal
   rate.
3. Run enough deterministic seeds to form uncertainty intervals.

**Accept only when:** the audit correctly identifies injected effects at
pre-specified thresholds and does not label pure site/leakage signal as
portable biology.

## Milestone 4: primary IBD benchmark

**Implement:**

1. Establish the shared-space, source-only transfer as the primary endpoint.
2. Benchmark CLR elastic net, tuned XGBoost/RF, shallow MLP, ComBat/MMUPHin
   variants where appropriate, and a SIAMCAT-style study-aware reference.
3. Use nested source-side tuning/early stopping consistently. Report AUROC,
   PR-AUC, calibration, subject/origin class counts, repeated-source-split
   distributions, cluster bootstrap CIs, and paired comparisons.

**Accept only when:** every table and figure is generated from one immutable
long-form metric table; no qualitative label is assigned to an effect whose CI
crosses zero or whose paired test fails the pre-specified threshold.

## Milestone 5: graph and biological interpretation

**Implement:**

1. Make the default graph analysis prior-free. Any supervised prior is fitted
   inside each training fold and reported as a nested sensitivity analysis.
2. Compare association types and pre-specified edge-density/top-k settings;
   quantify node/edge stability, graph overlap, and transfer performance.
3. Limit graph output to exploratory ecological hypotheses unless it improves
   the clean primary endpoint with statistical support.
4. Separate IBD, CD, UC, activity, therapy, and control-definition analyses
   wherever sample size allows; otherwise explicitly label the analysis coarse.

**Accept only when:** a graph claim names its exact data contract, graph
construction, stability measure, and independent validation result.

## Milestone 6: drug synergy extension

Drug synergy is a separate validated workstream, not a reinterpretation of IBD
classification output.

**Implement only after Milestone 4:**

1. Define the therapeutic setting, response assay, synergy label (Bliss, Loewe,
   HSA, or ZIP), drug representation, dose matrix, and train/test grouping.
2. Build a drug-pair response contract with cell/host context, microbiome or
   metabolite features, exposure, batch, and temporal provenance.
3. Use pair- and scaffold-aware splits to prevent chemical analogue leakage;
   evaluate calibration and external disease/context transfer.
4. Treat validated cohort-stable host-microbe-metabolite features as optional
   mechanism covariates, not surrogate synergy labels.

**Accept only when:** synergy performance exceeds appropriate single-agent and
chemical-only baselines on an external or held-out scaffold/context set, with
uncertainty and prospective validation candidates.

## Manuscript release gate

Before revision submission, regenerate every figure/table from frozen run IDs,
run a textual number-consistency test, publish the data dictionary and source
code needed for each panel, and replace the old claims with the scope justified
by Milestones 3--5.
