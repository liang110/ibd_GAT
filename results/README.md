# Study results

The [display index](display_sources.csv) links each figure and table to the data and programs used for its reported quantities.

| Directory | Contents |
| --- | --- |
| [comparisons](comparisons/) | Performance by target, source fit and scoring unit. Joint intervals resample target participants across ten fixed fits. |
| [source_data](source_data/) | Cohort denominators, shared taxa, model scores, source-pair AUROC contributions and participant-weighted Brier comparisons. |
| [controls](controls/) | Known-truth rules and each seed's decision. |
| [topology](topology/) | Inputs and fit-specific scores for the full-batch graph intervention and the node-identity sensitivity analysis. |

In the topology data, `graph_id` joins `edges.csv` to `graphs.csv`. Source and target are zero-based, directed node indices. Use the direction and split to resolve node indices in `features.csv`. Join to `edge_ablation_metrics.csv` on direction, split, graph condition and rewiring seed. Each direction and split has one Spearman graph, five rewired graphs and one self-loop graph.

The [verification commands](../repro/README.md) check graph identities against the recorded training results, compare edge constraints and rebuild taxon-aligned overlaps. The [training configurations](../repro/audits/training_configuration.json) distinguish the primary classifiers and the full-batch interventions.
