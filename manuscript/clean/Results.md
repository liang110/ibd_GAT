# Results

## An auditable transfer contract was fixed before model fitting

The study framework (Fig. 1) fixed the input, feature and evaluation contract before fitting any classifier. The HMP2/IBDMDB input contained 1,638 metagenomic samples from 130 participants and 578 raw microbial features; PRISM contained 220 samples from 220 participants and 209 raw microbial features. Median library sums were approximately 100 in HMP2 and 1 in PRISM, documenting a percentage-like versus proportion-like scale mismatch rather than a biological abundance comparison (Fig. 2a; Table 1). The loader closed both matrices within the shared feature space before source-fitted transformation.

Canonical harmonization produced 126 shared microbiome features. Three HMP2 samples had zero mass in the declared shared space and were excluded before CLR transformation; no PRISM sample was excluded. The exact-name metabolite intersection was zero. Primary cross-cohort analyses therefore used microbiome features only. The 126 features were the frozen contract universe; source-fold mutual-information selection retained up to 100 model features within each split (Fig. 2b,c; Supplementary Table 1). The same clinical and identifier denylist was applied to every model branch.

## Participant-aware and known-truth controls identified failure modes

Across ten repeated leakage controls, naive sample splits shared a mean of 43 participants between training and evaluation sets and produced a mean AUROC of 0.945. Participant-grouped splits had zero participant overlap and a mean AUROC of 0.464 (Fig. 3a; Table 2). These results quantify optimism caused by treating repeated samples from one participant as independent observations; they are a split-design control, not an estimate of external clinical performance.

The planted clinical leakage validator rejected all ten leakage cases. In known-truth semi-synthetic controls, null and batch-only scenarios produced mean AUROCs of 0.512 and 0.495, respectively. Signal-only scenarios recovered 5.1 of six injected disease features and produced AUROC 0.738; reversing-batch and topology-shift scenarios recovered 5.2 of six features and produced AUROCs 0.737 and 0.753 (Fig. 3b–d; Supplementary Table 2). Across the effect-size sweep, target AUROC was 0.644, disease-feature recovery averaged 3.3 of six and all acceptance criteria were met in 20% of settings. The domain-stable selector increased recovery to 4.1 of six and the all-criteria pass rate to 60%. These stress tests define the operating boundary of the selector and do not establish biological generalization.

## PRISM origin produced asymmetric transfer

PRISM discovery contained 155 samples (121 IBD and 34 controls) and validation contained 65 samples (43 IBD and 22 controls). In the five-repeat origin-held-out tabular benchmark, discovery-to-validation AUROC was 0.789 for logistic regression and 0.824 for random forest; validation-to-discovery AUROC was 0.881 and 0.888, respectively (Fig. 4a; Supplementary Table 3). Corresponding PR-AUC values were 0.899 and 0.914 for discovery-to-validation and 0.963 and 0.968 for validation-to-discovery.

Across ten origin-held-out source repeats, No-Prior Spearman GAT reached AUROC 0.718 ± 0.038 for discovery-to-validation and 0.752 ± 0.167 for validation-to-discovery (Fig. 4b). Across the ten source splits, the mean split-specific participant-cluster bootstrap ΔAUROC was −0.051 versus elastic net (mean split-specific interval endpoints −0.223 to 0.126) and −0.110 versus random forest (−0.235 to 0.012) for discovery-to-validation. In validation-to-discovery, corresponding means were −0.099 (−0.182 to −0.013) and −0.134 (−0.200 to −0.067). These endpoints summarize split-specific bootstrap intervals rather than a pooled cohort-level interval. The results show that PRISM performance depends on origin and transfer direction.

## Fair within-contract transfer separated relational and tabular performance

All models used the same 126-feature contract, source-only preprocessing and participant-grouped source splits. In HMP2 within-cohort evaluation, the legacy Full-Model GAT achieved AUROC 0.583 ± 0.142, whereas No-Prior GAT achieved 0.826 ± 0.055. No-Prior was therefore selected as the primary relational analysis (Supplementary Table 4). No-Context GAT achieved 0.629 ± 0.103, while the corresponding fusion achieved 0.791; these values do not support a general claim that the context branch improves portability.

## Directional cross-cohort transfer showed a conditional relational gain

In the ten-split confirmation, No-Prior Spearman GAT achieved AUROC 0.846 ± 0.035 for HMP2-to-PRISM transfer, compared with 0.715 ± 0.054 for elastic net and 0.768 ± 0.039 for random forest. GAT PR-AUC was 0.941, versus 0.875 and 0.907, respectively (Fig. 5b; Table 3). Across source splits, the mean paired ΔAUROC was 0.120 versus elastic net (mean split-specific interval endpoints 0.034–0.205) and 0.080 versus random forest (0.018–0.141). These are mean split-specific bootstrap quantities, not pooled intervals.

For PRISM-to-HMP2, No-Prior Spearman GAT achieved AUROC 0.612 ± 0.035, below elastic net at 0.711 ± 0.023 and random forest at 0.680 ± 0.011. Mean paired ΔAUROC was −0.100 versus elastic net (mean endpoints −0.229 to 0.026) and −0.067 versus random forest (−0.183 to 0.048). The summarized intervals included zero. The result is a direction-specific relational gain rather than a symmetric GAT advantage (Fig. 5c; Table 3).

## Topology and capacity controls localized the transfer difference

In the five-split topology comparison, HMP2-to-PRISM AUROC was 0.573 ± 0.091 for no-edge/identity controls, 0.587 ± 0.130 for rewired topology, 0.658 ± 0.162 for Pearson topology and 0.851 ± 0.007 for Spearman topology. For PRISM-to-HMP2, the corresponding values were approximately 0.539, 0.547, 0.513 and 0.607 (Fig. 6a; Table 4). The topology comparison therefore separated edge construction from the relational model itself.

The approximate capacity-matched no-message-passing MLP achieved AUROC 0.751 ± 0.061 and 0.708 ± 0.019 in the two directions, compared with Spearman GAT values 0.846 ± 0.035 and 0.612 ± 0.035 (Fig. 7; Table 4). The forward difference was therefore not fully explained by generic nonlinear capacity, whereas the reverse direction favoured the comparator.

Source-split edge export showed moderate Pearson–Spearman overlap: mean Jaccard 0.325 for HMP2-to-PRISM and 0.415 for PRISM-to-HMP2. Spearman–rewired overlap was 0.054 and 0.043, respectively (Fig. 6c; Supplementary Table 5). Because every retained node was assigned eight outgoing edges, out-degree vectors were constant and correlation was undefined; the audit artifact’s value 1.0 is reported only as a construction-level placeholder.

## Evidence supports a bounded framework-level conclusion

Across input audits, leakage controls, origin-held-out uncertainty, topology controls and capacity controls, relational gain was conditional on transfer direction and cohort structure. The evidence supports a source-only framework for evaluating relational and tabular microbiome representations and a conditional graph contribution in HMP2-to-PRISM transfer. It does not support a universal GAT advantage, a cross-cohort shared-metabolome claim or biomarker-level interpretation of graph-ranked features. Graph-ranked outputs are candidate, hypothesis-generating features.
