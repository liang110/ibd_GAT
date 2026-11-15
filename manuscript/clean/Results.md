# Results

## Input and analysis contract

The study began with a frozen input contract. HMP2/IBDMDB contributed 1,638 manifest samples from 130 participants and 578 raw microbial features; PRISM contributed 220 samples from 220 participants and 209 raw microbial features. Median library sums were approximately 100 in HMP2 and 1 in PRISM, documenting a percentage-like versus proportion-like scale mismatch (Fig. 2a; Table 1). The matrices were closed inside the shared feature space before source-fitted transformation.

Canonical harmonization produced 126 shared microbial features. Twelve HMP2 samples had zero mass in the declared shared space and were excluded before CLR transformation; no PRISM sample was excluded. The included HMP2 panel therefore contained 1,626 samples from 130 participants, while all 220 PRISM samples remained included. The exact-name metabolite intersection was zero. Exact, normalized and substring sample- and participant-ID overlap between the HMP2 and PRISM manifests was zero in both the pre-exclusion and included panels. The 126-feature registry is the frozen candidate universe, whereas source-fold mutual-information selection retained at most 100 coordinates for a fitted model (Fig. 2b–d; Supplementary Table 1). The contract is accompanied by a machine-readable manifest and source-table registry.

## Leakage and known-truth controls

Ten repeated leakage controls quantified the effect of repeated participants. Naive sample splits shared a mean of 43 participants between training and evaluation sets and produced a mean AUROC of 0.960. Participant-grouped splits had zero overlap and a mean AUROC of 0.552182 (Fig. 3a,b; Table 2). This comparison is a split-design control and is not the external HMP2 performance estimate.

The planted clinical-field validator rejected all 10/10 leakage cases. In the frozen `v2_prevalent_truth_features` known-truth benchmark, null and batch-only scenarios produced mean AUROCs of 0.512 and 0.495. Signal-only scenarios recovered 5.1 of six injected disease features with mean AUROC 0.738 and an 80% acceptance rate; reversing-batch and topology-shift scenarios recovered 5.2 of six with mean AUROCs 0.737 and 0.753, with acceptance rates of 70% and 100%, respectively (Fig. 3c; Supplementary Table 2). These experiments validate pipeline behaviour under declared perturbations rather than IBD biology.

## PRISM origin transfer

PRISM origin strata contained 155 discovery samples (121 IBD and 34 controls) and 65 validation samples (43 IBD and 22 controls). In the ten-split origin closure, discovery-to-validation AUROC was 0.718 ± 0.038 for No-Prior GAT, 0.770 ± 0.035 for elastic net and 0.827 ± 0.008 for random forest. Validation-to-discovery AUROC was 0.751 ± 0.167, 0.850 ± 0.024 and 0.885 ± 0.011, respectively (Fig. 4a,b; Supplementary Table 3). The mean split-specific paired GAT minus elastic-net difference was −0.052 for discovery-to-validation and −0.098 for validation-to-discovery; the corresponding GAT minus random-forest differences were −0.110 and −0.134. These values are source-split summaries of participant-cluster bootstrap estimates, not pooled intervals. The local data contract records discovery and validation labels but not a validated city or country field; no geographic causal attribution is made.

## Bidirectional transfer under the frozen registry

All models used the same 126-feature registry, source-only preprocessing and participant-grouped source splits. In ten source repeats, No-Prior Spearman GAT achieved AUROC 0.8565 ± 0.0063 for HMP2-to-PRISM transfer, compared with 0.7137 ± 0.0644 for elastic net and 0.7572 ± 0.0372 for random forest. Average precision was 0.9435 ± 0.0044, 0.8677 ± 0.0338 and 0.9032 ± 0.0163, respectively; Brier scores were 0.1409 ± 0.0051, 0.2881 ± 0.0438 and 0.1688 ± 0.0056 (Fig. 5a,c,d; Table 3). The B=5,000 participant-cluster paired bootstrap gave mean GAT-minus-elastic-net ΔAUROC +0.1425 with a split-aggregated 95% interval of [0.0605, 0.2266], and mean GAT-minus-random-forest ΔAUROC +0.0994 with [0.0364, 0.1640] (Fig. 5b).

For PRISM-to-HMP2, No-Prior Spearman GAT achieved AUROC 0.6138 ± 0.0370, below elastic net at 0.6660 ± 0.0256 and random forest at 0.6973 ± 0.0085. Average precision was 0.7973 ± 0.0122, 0.8361 ± 0.0135 and 0.8424 ± 0.0082, respectively; Brier scores were 0.1923 ± 0.0011, 0.2206 ± 0.0144 and 0.1984 ± 0.0047. The paired bootstrap intervals were [−0.1761, 0.0712] for GAT versus elastic net and [−0.1755, 0.0062] for GAT versus random forest. Thus, the forward advantage is stable under target-participant resampling, whereas reverse transfer does not support universal GAT superiority.

## Contract-aligned topology stability and capacity comparison

The final topology export used ten source splits in each direction and three graph constructions: Pearson, Spearman and destination-shuffled. Mean edge-set Jaccard values across split pairs were 0.0603, 0.0601 and 0.0528 for HMP2-to-PRISM, and 0.0441, 0.0469 and 0.0399 for PRISM-to-HMP2, respectively (Fig. 6a; Supplementary Table 5). Spearman edge records and edge counts matched the transfer metrics for all 60 rows checked. Each retained node had eight outgoing edges, so the degree vector had zero variance and conventional degree correlation was reported as NA under the prespecified policy. These diagnostics establish graph-construction provenance and stability; they are not evidence of biological interaction conservation.

The capacity-matched no-message-passing MLP used the same source-only transformed input boundary and paired source splits. Its AUROC was 0.7665 ± 0.0572 for HMP2-to-PRISM and 0.6794 ± 0.0217 for PRISM-to-HMP2, with average precision 0.8945 ± 0.0346 and 0.8368 ± 0.0127, respectively (Fig. 7a,b; Supplementary Table 6). The participant-level predictions supported a B=5,000 paired bootstrap: GAT-minus-MLP ΔAUROC was +0.0900 with a split-aggregated interval of [0.0177, 0.1631] forward and −0.0656 with [−0.1853, 0.0503] reverse. The forward contrast therefore exceeds a capacity-matched nonlinear reference under the frozen contract, while the reverse contrast again favours the reference on AUROC.

## Bounded framework-level interpretation

The evidence supports a source-only, directional framework for comparing relational and tabular microbiome representations. In this contract, relational gain is positive for HMP2-to-PRISM, is not reproduced by PRISM origin-held-out transfer, and is negative for PRISM-to-HMP2. The results do not establish universal GAT superiority, geographic effects, causal microbial interactions, graph-ranked biomarkers or clinical utility.



