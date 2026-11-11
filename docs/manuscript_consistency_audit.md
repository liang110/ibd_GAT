# Manuscript consistency audit

Source: `old_manuscript.docx (local read-only reference)` (read-only extraction)

Date: 2026-08-16

## Confirmed inconsistencies

| Location in extracted manuscript | Observation | Required correction |
|---|---|---|
| Figure anchors | The document contains implementation-style anchors `[Figure1_new]`, `[Figure2]`, `[Figure3]`, `[Figure4]`, `[Figure6_BIB_style_latest_clean]`, `[Figure6]`, and `[IBDFigure7]`. | Replace anchors with final numbered figures generated from frozen runs. |
| Baseline figure | Section 3.5 says the benchmark is visualized in Figure 5, while the nearby anchor is `[Figure6_BIB_style_latest_clean]` and its caption calls it Figure 5. | Create a figure manifest mapping final figure number, caption, source run, source CSV and rendered asset. |
| Portable feature figure | Section 3.6 refers to Figure 7A/Table 2A; the nearby anchor `[Figure6]` has a Figure 6 caption. | Reconcile prose, anchor and caption before any submission. |
| Cross-cohort values | Legacy text reports HMP2-to-PRISM XGBoost 0.8023 and PRISM-to-HMP2 0.6993, which belong to the invalid unharmonized pipeline. | Remove or label as legacy exploratory evidence; regenerate using shared-space run IDs. |
| Table 1 interpretation | Delta-AUC values -0.0009, -0.0022, -0.0030 and +0.0135 receive labels such as “Similar performance” and “Modest main gain” without uncertainty. | Replace Table 1 with long-form metrics and paired/bootstrap intervals; omit labels whose uncertainty crosses zero. |
| Target adaptation | Section 3.7 reports target-cohort adaptation as part of the framework story. | Move to a separate explicitly label-assisted deployment experiment or remove from the source-only portability paper. |
| Claim hierarchy | The manuscript alternates among graph method, portability audit, benchmark, biology and adaptation. | Reorganize as known-truth validation; primary clean transfer; origin/confound audit; exploratory graph hypotheses; limitations. |
| References | Reviewers identified a wrong PRISM author list and uncited references. | Validate every bibliography entry against primary sources and run a cite-key/in-text citation audit. |

## Current publication rule

No number, figure, table or claim in the DOCX is eligible for reuse until it is
linked to a committed run ID, exact source table and immutable plotting command.
The existing DOCX remains a legacy submission record, not a source of current
results.

## Required release artifacts

1. `figure_manifest.csv`: figure/table identifier, source run ID, CSV/JSON,
   command, Git revision and manuscript section.
2. A single generated `results_long.csv` with metrics, class/sample/participant
   counts and uncertainty metadata.
3. Textual numeric-consistency test that compares manuscript values with the
   frozen source table.
4. Verified reference database and an in-text citation coverage report.
