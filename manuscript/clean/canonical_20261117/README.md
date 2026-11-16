# Manuscript reconstruction source

This directory is the manuscript-facing source for the HMP2/IBDMDB–PRISM framework paper. It is generated from the Nov-2026 evidence ledger under a single 126-feature microbiome contract. The main scientific claim is a directional, conditional relational gain; the package does not claim universal GAT superiority, causal microbial interactions, a shared metabolome or clinical utility.

The source order is Abstract, Introduction, Methods, Results, Discussion and Conclusion. Figure 1 belongs to Methods; Figures 2–7 are placed at their first evidence point in Results. Equations are converted to native Word OMML by artifacts/build_submission_docx.py.

figure_table_registry.csv is the machine-readable source registry. docs/reconstruction_gap_matrix_20261112.md records the plan-to-evidence audit and docs/reconstruction_todo_20261112.md records experiments that are useful but not required for the bounded two-cohort claim.


## Evidence closure (2026-11-17)

This source tree is paired with `artifacts/canonical_evidence_20261117/`. It uses the reconciled contract hash `a18364d4...`, the 18,460-row XGBoost ledger, the ten-split topology freeze and the 47-reference citation closure. Legacy source CSVs listed in `manuscript_rebuild/source_data/source_data_registry_20261117.csv` are quarantine-only.
