# Reviewer comment coverage audit

Date: 2026-08-16

This document maps every supplied reviewer concern to evidence in the rebuilt
framework. `Resolved` means executable evidence exists; it does not mean that
the manuscript text has already been revised. `Partial` and `Open` items remain
release blockers when they affect a central claim.

| ID | Reviewer concern | Assessment | Current status | Evidence and remaining action |
|---|---|---|---|---|
| R1.1 / R2.5 | GAT external collapse and fusion choosing tabular branch | Fully reasonable | Partial | Legacy GAT is no longer a primary predictor and graph claims are blocked. A prior-free graph sensitivity analysis is not yet run; alternatively remove GAT as a central contribution. |
| R1.2 | MI feature selection may use full-cohort labels | Reasonable audit request; latest legacy code did not confirm the suspected leak | Resolved in rebuilt path | Prevalence, pseudocount, CLR, scaling and MI are fitted per source-training split; planted leakage tests pass. Manuscript still needs exact boundary language and fold feature-list provenance. |
| R1.3 | Generic baselines rival or beat the framework | Fully reasonable | Substantially addressed | Clean shared-space logistic-L2, elastic net, RF and source-validation early-stopped XGBoost were rerun. The Stage 17 domain-interaction sensitivity candidate has HMP2-to-PRISM paired advantage versus elastic net in 4/5 splits, but no consistent PRISM-to-HMP2 advantage and different eligible source-site sets. MLP remains open; no general performance-superiority language is allowed. |
| R1.4 | Pearson 0.15 / top-k 8 graph edges lack justification | Fully reasonable | Open | No new graph claim is allowed. If graph work continues, compare Pearson, Spearman and compositional association methods with density, overlap and stability audits. |
| R1.5 | Figure-panel references do not agree | Fully reasonable | Audited; correction open | Read-only DOCX audit confirms placeholder anchors and caption/prose number drift. Requires figure regeneration from frozen run IDs and an automated figure-number manifest. |
| R1.6 | Exact post-harmonization dimensions are missing | Fully reasonable | Partial | External microbiome dimension is now 126 shared taxa; exact shared metabolites are 0. Within-cohort and post-prevalence dimensions must be exported per run and added to Methods. |
| R2.1a | Method lacks known-confound validation | Fully reasonable and central | Substantially addressed; real-data evidence remains directional | Calibrated v2 controls pass null 10/10, batch-only 10/10, signal-only 8/10, topology-shift 10/10 and reversing-batch 7/10. Repeated-subject splitting and planted clinical leakage pass 10/10. The shared-effect/domain-deviation candidate passes 8/10 reversing-batch seeds and, on real data, has HMP2-to-PRISM paired benefit but no consistent reverse-direction benefit. Its strict all-HMP2-site analysis fails at Emory, so it remains a labelled sensitivity analysis. |
| R2.1b | Add another independent biological dataset | Reasonable for strong generalization claims | Open / data-dependent | Not solved by HMP2 and pooled PRISM. Identify and contract a third IBD cohort or narrow the claim to two-cohort benchmarking. |
| R2.2a | Sample size is small for GNN/fusion complexity | Fully reasonable | Partially addressed | Complex GNN is not promoted; simple CPU baselines and cluster uncertainty are primary. Formal learning-curve/model-stability analysis remains open. |
| R2.2b | IBD endpoint mixes CD, UC, activity, therapy and controls | Fully reasonable | Substantially addressed within available metadata | Separate clean CD- and UC-versus-control transfers, with PRISM origins kept separate and target-cluster uncertainty, are complete. Mean AUROC ranges from 0.680 to 0.787 for CD and 0.621 to 0.725 for UC, so the coarse endpoint cannot stand in for both subtypes. Activity and therapy fields are absent in both contracts; those questions remain impossible without new versioned metadata and causal claims remain prohibited. |
| R2.3 / R3.6 | Text, figures and tables contain conflicting numbers | Fully reasonable and publication-blocking | Audited; correction open | DOCX audit confirms figure-anchor drift, legacy invalid transfer numbers and tiny unqualified delta-AUC labels. New results have immutable CSV/JSON sources, but the DOCX must be regenerated and numerically tested. |
| R2.4 | Manuscript mixes method validation, prediction, audit, graph and biology | Fully reasonable | Audited; correction open | DOCX audit confirms mixed framing and target-adaptation placement. Revision structure should separate known-truth validation, clean transfer benchmark and exploratory biology. |
| R3.1a | HMP2 percentages and PRISM proportions were not harmonized | Fully reasonable and publication-blocking | Resolved in rebuilt path | Both cohorts are closed in the declared shared space before training-only CLR. Legacy transfer numbers remain withdrawn. |
| R3.1b | Most target graph nodes were absent and silently zero-filled | Fully reasonable | Resolved for primary benchmark | Models require 126 identically ordered shared taxa and raise on missing fitted features. Graph benchmark remains blocked. |
| R3.1c | Shared metabolite count may be zero | Fully reasonable | Confirmed; transfer blocked | Exact shared metabolite names are 0. No metabolome transfer claim is allowed until HMDB/KEGG/InChIKey mapping is completed. |
| R3.2 | PRISM discovery and Dutch validation cohorts were pooled without audit | Fully reasonable | Substantially addressed | Origin is reconstructed (155 discovery, 65 validation). Five-split discovery-to-validation and validation-to-discovery transfers plus target bootstrap are complete. Country, sequencing-run and metabolomics-batch provenance remain incomplete, so causal attribution to US/NL or platform is not allowed. |
| R3.3a | PRISM baselines contain Age/Fecal Calprotectin | Fully reasonable and reproduced | Resolved in rebuilt path | The common pipeline hard-rejects clinical, outcome and identifier fields; clean baseline runs replaced the contaminated legacy comparisons. |
| R3.3b | Baselines are not tuned equivalently | Fully reasonable | Substantially addressed | Logistic-L2, elastic net, RF and XGBoost all use the same source-training matrix and source-validation selection; XGBoost has source-validation early stopping. Fully nested tuning remains open. |
| R3.3c | CLR, batch-corrected and SIAMCAT-style references are missing | Fully reasonable | Partial | CLR logistic-L2, elastic net, RF and XGBoost are present. ComBat/MMUPHin and SIAMCAT-style leave-one-study-out references are not yet implemented. |
| R3.4 | Whole-HMP2 supervised dysbiosis prior leaks labels and target information | Fully reasonable and publication-blocking | Resolved for primary benchmark | The rebuilt primary benchmark is prior-free. A supervised prior may return only if estimated inside each training split and reported as secondary sensitivity analysis. |
| R3.5 | Novelty is overstated and direct precedents are under-discussed | Fully reasonable | Open in manuscript | Reframe as a rigorous portability/confound audit benchmark; update SIAMCAT, Ning et al. and confound-adjustment literature with verified citations. |
| R3.6a | No confidence intervals, paired tests or PR-AUC | Fully reasonable | Substantially addressed | PR-AUC, Brier, target-cluster bootstrap CIs and paired model-difference bootstrap intervals are complete. A manuscript-level multiplicity policy and a pooled source-split inference plan remain open. |
| R3.6b | Repeated HMP2 samples overstate effective sample size | Fully reasonable | Substantially addressed | Participant-group splits and participant-cluster bootstrap are implemented. Known-truth repeats show random-split AUROC 0.945 versus grouped 0.464. Every manuscript table still needs sample, participant and class counts. |
| R3.O1 | Single split versus repeated-split SD is inconsistent | Fully reasonable | Resolved analytically; manuscript open | Five source splits with fixed targets are used; SD is source-training variability and bootstrap is target-sampling uncertainty. The distinction must be written into Methods/results. |
| R3.O2 | Environment, license, inputs and figure-source results are incomplete | Substantially reasonable | Partial | New environment files, Git history, tests and frozen legacy checksums exist. Exact lockfile, license, acquisition recipes for restricted inputs, and figure-to-source manifest remain open. |
| R3.O3 | PRISM citation metadata and uncited references are wrong | Fully reasonable | Open | DOCX audit retains this issue; validate the bibliography against primary sources and remove uncited entries during manuscript revision. |

## Coverage conclusion

The current executable plan addresses the most severe computational validity
defects, but it does **not** yet resolve every reasonable reviewer request. The
remaining publication blockers are: reporting of the completed conservative
CD/UC results plus narrowed activity/therapy claims; third-cohort evidence or
narrowed claims; graph removal/reframing or stability
validation; metabolite ID mapping or explicit omission; fully nested tuning and
additional reference methods if a performance claim remains; manuscript
number/figure/citation consistency; and a complete release package. Known-truth
controls, PRISM-origin directions, fair baselines and paired inference now exist,
but the domain-interaction result is a direction-specific sensitivity result,
not evidence of universally portable superiority.

## Execution order

1. Validate the synthetic test with positive/negative/topology controls.
2. Add repeated-subject and planted-leakage recovery scenarios.
3. Implement and compare domain-aware candidates only after controls pass.
4. Run PRISM discovery-to-validation and validation-to-discovery baselines.
5. Complete fair comparator and paired-inference tables.
6. Report the completed conservative CD/UC analyses; narrow activity and therapy claims.
7. Decide whether to remove the graph branch or execute the full stability matrix.
8. Regenerate manuscript tables/figures and validate references from frozen runs.
