**Figure 1. Source-only information-flow contract.** HMP2/IBDMDB and PRISM are projected to a frozen 126-feature microbial registry; source-training operations produce the same input boundary for No-Prior GAT, elastic net and random forest; target labels enter only after frozen predictions.

**Figure 2. Input scale, feature contract and exclusion audit.** (a) Library-sum ranges and medians; (b) raw microbial features to the shared 126-feature registry and zero exact-name metabolite intersection; (c) source-only transformation boundary; (d) sample and zero-mass exclusion counts.

**Figure 3. Participant leakage and known-truth controls.** (a) Current-contract seed-level sample-level versus participant-grouped AUROC (0.960 versus 0.552182; 10 seeds); (b) participant overlap for the same seeds (mean 43 versus 0 in every grouped split); (c) null, batch-only, signal-only, reversing-batch and topology-shift scenario results. Synthetic panels validate the procedure and are not biological cohort validation.

**Figure 4. PRISM origin-held-out transfer.** (a) Discovery/validation composition and source-target directions; (b) ten-split AUROC by model; (c) split-specific paired GAT minus reference differences; (d) average precision. Origin labels are discovery/validation because the local manifest does not provide a validated geographic field.

**Figure 5. Bidirectional transfer under the frozen registry.** (a) ten source-repeat AUROC points; (b) paired ΔAUROC summaries; (c) model and direction means; (d) average precision.

**Figure 6. Pre-convolution input-edge configuration sensitivity and edge diagnostics.** (a) ten source-split Real, Rewired and Identity AUROC conditions; (b) direction-by-configuration means; (c) taxon-ID-aligned directed edge-set Jaccard; (d) fixed-k degree diagnostic, which is undefined because out-degree variance is zero. The edge audit does not infer convolution-internal self-loop handling.

**Figure 7. Approximate no-message comparator and probability error.** No-Prior Spearman GAT and no-message MLP use the same source-only transformed input boundary. Participant-level predictions support the paired target-participant bootstrap; Brier differences remain learner-specific and are not standalone calibration evidence.

