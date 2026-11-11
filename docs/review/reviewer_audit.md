# Reviewer and framework-article audit for the next reconstruction

## Source material reviewed

- `old_manuscript.docx (local read-only reference)` (old manuscript structure, author information, title-page style, and original framework narrative).
- `D:/CodexHome/attachments/15c8c710-1e9a-4447-9811-2cd63706892a/已粘贴的文本.txt` (Reviewer 1–3 reports).
- Current manuscript-facing sources in `manuscript_rebuild/clean/` and the frozen evidence/figure registry.

## What the reviewers accepted

All three reviewers regarded the problem as important and considered the framework direction potentially valuable. Their objections were mainly about validity, boundary definition, like-for-like comparison, inconsistent reporting and overstatement of novelty. Reviewer 1 explicitly described the dual-engine comparison and signal-decomposition goal as valuable. Reviewer 2 accepted the potential contribution but requested known-truth validation and a clearer separation of method validation, cross-cohort prediction and biological interpretation. Reviewer 3 agreed that the confound/portability question is timely and that the authors were unusually transparent about graph underperformance, but required correction of the original transfer, leakage and baseline problems.

The revised evidence package now addresses the most direct implementation objections: the cross-cohort primary analysis uses a frozen 126-feature microbiome universe, source-only closure/CLR and selection, a zero exact-name metabolite intersection is reported, participant grouping is explicit, the HMP2/PRISM direction is evaluated in both directions, PRISM origin-held-out transfer is included, clinical denylist handling is shared, synthetic controls are labelled known-truth stress tests, and split-specific paired bootstrap/topology/capacity controls are archived.

## Required reviewer closure map

### Reviewer 1

1. **Graph structural paradox.** The paper must state that the graph branch is an object of diagnosis rather than a universally superior classifier. The revised narrative should present the measurable contribution as an estimand: under a frozen contract, what is the directional performance difference between a relational representation and strong marginal references? The observed forward gain and reverse/origin failure are central evidence for direction dependence, not a reason to claim conserved graph biology.
2. **Validation boundary and MI selection.** Methods must state the order: canonical feature universe first; source-training prevalence, pseudocount, scaling and MI selection inside each source split; target is transformed only after the source object is frozen. The 126-feature value is the harmonized universe, whereas the model receives at most 100 fold-local selected features.
3. **Baseline equivalence.** Every branch must share the same denylist, source-only transformed input and participant-aware splits. Discussion should acknowledge the practical trade-off: a tabular model can be simpler and more portable; the framework's value is controlled decomposition, not raw accuracy dominance.
4. **Edge assignment.** The revised Methods should distinguish the absolute association matrix, fixed top-k directed edge set, rewired null and the fact that fixed top-k makes conventional degree correlation undefined. Pearson/Spearman topology is an audited representation choice, not evidence of ecological interaction. Avoid implying that |r| or k=8 alone establishes biological meaning.
5. **Figure and dimensionality alignment.** Captions must match the actual panels. State raw HMP2/PRISM dimensions, the 126 shared microbiome universe, three HMP2 zero-mass exclusions and zero exact-name metabolite overlap in the input contract section.

### Reviewer 2

1. **Known-truth validation and external application.** The synthetic null, batch-only, signal-only, reversing-batch and topology-shift experiments should be introduced as method-validation gates, not biological validation. PRISM origin-held-out transfer is the additional biological generalization test available within the current evidence package. Do not claim that two public cohorts establish broad clinical generalization.
2. **Small cohorts and heterogeneous endpoint.** Use participant/split-level uncertainty and source-split variability. Explicitly define the binary endpoint as IBD versus non-IBD and state that the experiment does not resolve CD/UC, activity, treatment or longitudinal heterogeneity. This is a scope boundary in Discussion, not a reason to inflate the result.
3. **Numerical consistency.** All numbers must come from the frozen registry. Each result paragraph should cite one figure/table and use the same direction/model/estimand label. Never mix old full-prior or metabolome-inclusive values with the No-Prior microbiome-only primary analysis.
4. **Narrative separation.** The manuscript should have three explicit phases: (i) contract and method validity, (ii) directional transfer and representation comparison, (iii) interpretation boundary and candidate biological outputs. The Abstract and Introduction should state that the paper is a framework/evidence-standard article with a worked IBD application.
5. **Graph contribution.** Reframe graph attention as a tested relational representation. The conclusion should be that relational gain is conditional and asymmetric, while the framework exposes when a graph branch fails to transport.

### Reviewer 3

1. **Scale and shared feature space.** Preserve the explicit closure/CLR equations, the common 126-feature space, the zero-mass exclusion rule and the zero metabolite intersection. This directly answers the original abundance-scale and zero-node criticism.
2. **PRISM pooled-origin confound.** The PRISM discovery/validation origin audit must appear as a named Methods estimand and a dedicated Results section. Report discovery-to-validation and validation-to-discovery separately; do not describe PRISM as a single homogeneous cohort.
3. **Like-for-like baselines.** State that the denylist and source-only transformation are common to every comparator. If a legacy result used a different input boundary, it must be excluded from the primary table and identified only in provenance notes.
4. **Global biological priors.** The primary model must be called No-Prior. Any legacy whole-cohort prior is a sensitivity/negative-control artifact, never evidence for the primary claim. Methods and Discussion should state that the No-Prior branch removes this target-leaking asymmetry.
5. **Novelty and literature.** The Introduction must cite SIAMCAT and prior cross-cohort IBD multi-omics work as direct precedents. The residual contribution is a reproducible audit contract plus a directional relational-versus-tabular estimand with known-truth, topology and capacity gates. It is a framework contribution, not a new universal GAT architecture.
6. **Statistical uncertainty and effective sample size.** Report participant-cluster bootstrap intervals, source-split variability, AUROC and average precision, and explain that ten-split summaries average split-specific quantities and interval endpoints rather than forming a pooled cohort-level CI. Report participants and class counts where available.
7. **Reproducibility.** The final package should expose manifest/hash/model-configuration/Slurm provenance and a figure/table source registry. The manuscript can state that repository deposition is pending if the GitHub repository has not yet been created.

## Framework-article writing plan

### Abstract

Use a structured framework pattern in four compact moves:

1. **Problem.** Cross-cohort microbiome classifiers can conflate disease signal with scale, participant repetition, clinical leakage, cohort origin and relational topology.
2. **Framework.** Define a source-only contract: shared feature universe, source-fitted transformations/selection, participant-grouped splits, directional transfer, known-truth controls, topology nulls, capacity control and paired uncertainty.
3. **Worked application.** State the IBD cohorts and the most decision-relevant results: 126 shared microbiome features, zero exact-name metabolites, leakage collapse after grouping, conditional HMP2-to-PRISM relational gain, and reverse/origin failure.
4. **Interpretation.** Claim an auditable evidence standard and a conditional relational estimand; do not claim a universally superior GAT or broad biological generalization.

Avoid leading with a list of model components or a performance leaderboard. The abstract should answer: what problem does the framework make measurable, what evidence gates does it impose, what did the worked case reveal, and what can a reader safely conclude?

### Introduction

Use the old manuscript's logical progression but sharpen it as a framework article:

1. IBD multi-omics produces reproducible but cohort-sensitive signals; within-cohort accuracy is not a portability claim.
2. Existing practice often combines heterogeneous feature spaces, repeated participants, metadata leakage and graph priors; SIAMCAT and cross-cohort studies establish the need for external/confound-aware evaluation.
3. A relational branch is scientifically plausible but can encode cohort-specific topology; a tabular branch is a necessary marginal reference. The question is not whether GAT is generally better, but under which contract and transfer direction relational information adds value.
4. State the framework design principles: frozen contract, source-only fitting, participant-aware estimands, known-truth validation and explicit interpretation gates.
5. State contributions as framework outputs: an auditable contract/manifest; a directional relational-versus-tabular estimand; a set of validity, topology and capacity controls; and a worked IBD case that exposes conditional, asymmetric graph utility.
6. End with a short roadmap separating method validation, biological transfer and interpretation boundaries.

Do not call the work a new GAT architecture. Do not claim that candidate taxa are validated biomarkers or that synthetic controls prove biology.

### Methods first section and Figure 1

Figure 1 belongs immediately after the `Methods` heading and the first short `Study design and cohorts` paragraph, before detailed preprocessing. Its caption should define the whole contract and the evidence gates. It should not appear in Introduction. Results figures should be inserted after the paragraph that first reports the corresponding result: Fig. 2 after the input-contract result; Fig. 3 after leakage/known-truth controls; Fig. 4 after origin-held-out transfer; Fig. 5 after bidirectional primary transfer; Fig. 6 after topology ablations; Fig. 7 after capacity control.

## Formula requirements for the next Methods rewrite

The Methods should expose the framework estimands, not only describe implementation:

- closure on the shared microbiome space and CLR transformation;
- source-only transformation operator and fold-local selector;
- directional transfer estimand, with source model fitted on S and target labels used only for final evaluation;
- participant-cluster bootstrap metric and paired difference;
- leakage optimism gap (naive AUROC minus grouped AUROC);
- topology construction and rewired null;
- edge-set Jaccard;
- capacity-adjusted relational gain;
- explicit statement that fixed top-k degree correlation is undefined.

Each equation needs a sentence defining every symbol and the corresponding figure/table artifact. Equations should be native Word equations or clean Unicode math with consistent subscripts/superscripts; raw LaTeX delimiters must never be visible in the DOCX.

## Reference expansion plan (minimum 40)

The new manuscript should include at least 40 cited references, distributed across these clusters:

1. HMP2/IBDMDB, PRISM and IBD multi-omics data (at least 5).
2. Cross-study microbiome ML and confound auditing, including SIAMCAT and related benchmark/adjustment work (at least 8).
3. Compositional data, CLR/ALR and microbiome normalization (at least 5).
4. Graph attention/GATv2, graph representation learning and microbiome network learning (at least 8).
5. Domain shift, external validation, leakage and participant-aware/nested resampling (at least 6).
6. Bootstrap, AUROC/average precision, calibration and reporting guidelines (at least 4).
7. IBD biological context, dysbiosis, metabolomics and CD/UC heterogeneity (at least 5).

Every reference must be cited in the text, have verified author/title/journal/year/DOI metadata, and be linked to a concrete claim. Do not add uncited padding references. The existing six references are a starting core, not a sufficient bibliography.

## Decision

The current evidence package can support a framework-style reconstruction, but only if the paper's central claim is narrowed to an auditable, directional representation-comparison framework. It cannot support a claim that GATv2 is a universally portable or biologically validated IBD classifier. The next reconstruction should be evaluated against the reviewer map above before generating the final DOCX.
