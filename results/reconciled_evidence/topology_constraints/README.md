# Topology constraint audit

`edge_constraint_audit.csv` is a historical compact edge audit for exported association/identity edge lists. It records zero self-loops and zero Jaccard overlap with the explicit identity self-loop graph. It is not a substitute for the formal Gate B rewired edge export.

`edge_hash_audit.csv` is the reconstruction ledger for the formal identity-aware run. It contains the source split and rewiring seed for each record. The current evidence package contains 24 hash mismatches (HMP2-to-PRISM split 9 and PRISM-to-HMP2 splits 6–8). The validator fails closed; these records must be regenerated or reconciled before claiming complete edge-level auditability.

`topology_audit_status.json` records this status in machine-readable form. The compact metrics remain sufficient for summary-table and figure-source reconstruction, but raw cohort matrices, participant identifiers, checkpoints and scheduler logs are excluded.
