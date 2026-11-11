# Legacy Freeze and Migration Manifest

Frozen source directory:
`${IBD_LEGACY_ROOT}`

Freeze date: 2026-08-15. The source directory is retained in place and is not
modified by this migration. New work is prohibited there.

Migrated evidence:

- `src/scripts/` and `scripts/` into `legacy/code/`.
- Four final Slurm designs into `legacy/workflows/`.
- `environment.yml` and `requirements.txt` into `environment/`.
- Entire legacy `results/` and `tables/` trees into `artifacts/legacy_*`.

Not copied:

- raw/derived input matrices, including `manual_abundance.csv`,
  `metabolites.csv`, `HMP2_metabolomics.biom`, PRISM matrices, spreadsheets,
  and Nextflow work directories;
- obsolete CKD documentation, broad source modules, archived trial jobs, and
  binary model checkpoints.

These exclusions avoid duplicating large data and prevent misleading legacy
material from becoming the new implementation. Future data access must be
declared through a data catalog with checksums and immutable source paths.

Legacy results must be cited as `legacy_unharmonized`. They cannot be used to
support a final cross-cohort multi-omics or portable-GAT claim until the new
shared-space and leakage-audited pipeline reproduces them under the revised
analysis contract.

