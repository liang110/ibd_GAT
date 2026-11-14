# Formal topology audit

The canonical audit tables record 140 edge lists: 20 real Spearman graphs, 20 identity graphs and 100 degree-preserving rewired graphs. `formal_edge_constraint_audit.csv` records node, edge, self-loop, duplicate-edge, degree and identity-overlap constraints. `formal_gateb_hash_audit.csv` compares each exported edge hash and selected-feature hash with the persisted formal run records. All 140 records match.

The individual edge CSVs are stored in `releases/formal_topology_edge_export_20261114.zip` rather than as 140 separate repository files. The archive also contains the canonical audit tables and manifests. The original 24 local reconstruction mismatches were resolved by exporting from the persisted split checkpoints in the original run environment; the row-level resolution is recorded in `mismatch_resolution.csv`.
