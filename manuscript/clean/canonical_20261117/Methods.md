# Methods

## Study design, cohorts and analysis units

We designed the study as a source-only, participant-aware comparison of marginal and relational microbiome representations for directional transfer. The primary cohorts were HMP2/IBDMDB and PRISM [1,2]. The endpoint was IBD versus non-IBD control; subtype, activity and therapy effects were not estimated. HMP2/IBDMDB contained 1,638 metagenomic samples from 130 participants and 578 raw microbial features. PRISM contained 220 samples from 220 participants and 209 raw microbial features. PRISM was retained as discovery (155 samples: 121 IBD and 34 controls) and validation (65 samples: 43 IBD and 22 controls) origin strata. The local manifest records these origin labels but does not contain a validated city or country field; the analysis therefore reports discovery/validation rather than assigning geographic labels. Participant identifiers, repeated-sample structure, class counts, origin labels and the input manifest were frozen before model fitting.

The framework is a linked set of estimands and evidence gates rather than a new neural layer. Figure 1 gives the information-flow contract: the shared microbial registry and metadata boundary are fixed first, source-training operations are fitted second, predictions are generated on a locked target third, and uncertainty and validity gates are evaluated last.

## Figure 1: information-flow contract

The figure shows the exact order used by every main analysis. HMP2 and PRISM microbial matrices are projected to the label-free registry U. Participant-grouped source splits determine the training participants. Closure, prevalence filtering, source-training pseudocount and CLR/scaling are fitted within each split; mutual-information selection and source-association topology are also split-local. The resulting matrix is sent to the same elastic-net, random-forest, no-message-passing MLP, XGBoost and No-Prior GATv2 comparison objects [27–36]. Target labels are introduced only after frozen predictions are written. Origin strata, target participant bootstrap, topology controls and capacity controls are evaluation layers and do not feed back into fitting.

## Shared feature registry and exclusion rules

Microbial identifiers were normalized with a conservative canonicalizer, and identifiers mapping to the same canonical taxon were summed. Clinical, outcome, medication, antibiotic, steroid, age, faecal-calprotectin and sample/participant identifier fields were removed by one shared denylist before any model branch was built. The label-free cross-cohort registry is the deterministic intersection of the post-canonicalization registries:

\[
U = C_{HMP2} \cap C_{PRISM}, \quad |U| = 126.
\]

Here U is the frozen candidate universe, not the number of coordinates presented to every split-local model. Twelve HMP2 samples had zero total mass in U and were excluded before transformation; no PRISM sample was excluded. The exact-name metabolite intersection was zero, so the primary cross-cohort estimand is microbiome-only. Any multi-omics fields remain outside the cross-cohort model contract.

## Compositional transformation and source-fold fitting

For sample x with non-negative coordinates in U, closure was evaluated after zero-mass exclusion:

\[
h_j(x) = \frac{x_j}{\sum_{\ell \in U}x_{\ell}}, \quad j \in U.
\]

For source s, split r and split-local coordinate set A_r, the pseudocount was δ_r = 0.5 times the smallest positive source-training value. The centred log-ratio transformation was:

\[
CLR_j(x;A_r,δ_r) = \log(h_j(x)+δ_r) - \frac{1}{|A_r|}\sum_{\ell \in A_r}\log(h_{\ell}(x)+δ_r).
\]

The fixed order was canonicalization and intersection, zero-mass exclusion, closure, source-training prevalence filtering, source-training pseudocount, source-training scaling and source-training mutual-information selection. Prevalence below 0.10 was evaluated on source-training participants only. At most 100 features with the largest source-training mutual information with the endpoint were retained; F_r therefore satisfies |F_r| ≤ 100 while U remains fixed at 126. The selected feature order, fitted transform and split manifest were archived.

For each directed source s to target t and source split r, the fitted object and deployment operator were:

\[
T_{s,r} = Fit(X^{train}_{s,r},y^{train}_{s,r}), \quad \hat y_t = f_{m,s,r}(T_{s,r}(X_t)).
\]

The target label, target prevalence, target ranking, target threshold and target fusion weight were absent from T_{s,r}; no domain-adversarial adaptation layer was added [37]. Constants and non-finite association values were excluded from topology ranking; ties were resolved by canonical taxon identifier.

## Model comparison objects

The tabular references received the same T_{s,r}(X) matrix and denylist as the graph branch. Elastic net used C=1, l1 ratio 0.5 and class balancing; random forest used 300 trees, minimum leaf size 2 and class balancing. XGBoost used 300 histogram trees with depth 2/3, learning rate 0.03/0.1, row and column subsampling 0.8, no early stopping and no class weighting. No-Prior GATv2 used two graph-attention layers with hidden width 16, two heads in the first layer and one in the second, dropout 0.1, global mean pooling and a two-layer classifier head. Adam optimization used learning rate 10^-3, weight decay 10^-3, 40 epochs and batch size 32. Topology-intervention and identity-aware sensitivity runs are reported as separate configurations: they use the same loss, optimizer and epoch budget but full-source-panel training with shuffling disabled; identity-aware runs additionally use an 8-dimensional learnable node-identity representation. Absolute differences between these configurations are therefore not interpreted as a pure topology effect.

Association matrices were computed from source-training profiles. Each retained node received its eight strongest non-self absolute associations under the declared rule:

\[
E_{s,r} = \{(i,j): j \in Top_k(|R_{s,r,i,j}|; \text{stable taxon-ID ties}), j \ne i\}, \quad k=8.
\]

Pearson and Spearman association supplied topology sensitivity; the topology-stability export also included a destination-shuffled null. No-edge and identity graphs were retained as secondary controls in the source package but were not mixed with the final contract-aligned topology summary. Edge construction is a representation choice and is not interpreted as a biological interaction.

The capacity comparator was an approximate no-message-passing multilayer perceptron using the same T_{s,r}(X) input and a similar hidden-parameter budget. It was paired by source split and target participant set, and participant-level predictions were archived for the B=5,000 comparison.

## Directional estimands and uncertainty

For metric M, target evaluation layer h_t, target weights w_t, source s, target t and source split r, the primary paired contrast was:

\[
\Delta M^{w_t}_{s\to t,r}(h_t;f,g) = M^{w_t}_{s\to t,r}(h_t;f) - M^{w_t}_{s\to t,r}(h_t;g).
\]

AUROC and average precision are higher-is-better metrics; Brier score is lower-is-better. Differences are always reported as GAT minus the named reference. Source-repeat standard deviation describes variation across fitted source splits. It is not a target confidence interval.

For each source split, B=5,000 bootstrap draws resampled target participants and retained all samples from a selected participant. The same participant indices were applied to all models and references before the paired difference was calculated:

\[
\Delta M^{(b)}_{s\to t,r} = M^{(b)}_{s\to t,r,GAT} - M^{(b)}_{s\to t,r,reference}.
\]

The 2.5th and 97.5th percentiles of the paired draw distribution define the split-specific 95% interval [39]. AUROC follows the standard ranking interpretation [40,41], and average precision is reported because the cohorts are IBD-skewed [42]. Ten-split summaries are arithmetic means of split-specific point estimates and interval endpoints; they are not pooled cohort-level intervals. The baseline GAT-versus-elastic-net/random-forest analyses completed 40 comparisons with zero invalid draws, and the capacity-matched GAT-versus-MLP analysis completed 20 comparisons with zero invalid draws. Source-side model selection is kept separate from target uncertainty [38].

## Validity gates and graph diagnostics

The leakage gate compares naive sample splitting with participant-grouped splitting and tests planted clinical fields. The known-truth gate evaluates null, batch-only, signal-only, reversing-batch and topology-shift scenarios with declared injected features and acceptance rules. The origin gate keeps PRISM discovery and validation separate.

The split-matched optimism gap was defined as:

\[
\Omega_{leak,r} = AUROC_{naive,r} - AUROC_{grouped,r}.
\]

For edge sets E_a and E_b, Jaccard overlap was computed by canonical taxon ID and directed edge convention:

\[
J(E_a,E_b) = \frac{|E_a \cap E_b|}{|E_a \cup E_b|}.
\]

An empty union is reported as undefined. Because fixed k=8 assigns the same out-degree to every retained node, Var(d)=0 and conventional out-degree correlation is undefined. The contract-aligned topology export used ten source splits and verified the exported Spearman edge records against the transfer metrics for all 60 direction/split/model rows.

The capacity contrast was:

\[
\Delta^{M,w_t}_{cap,s\to t,r}(h_t) = M^{w_t}_{s\to t,r}(h_t;No\text{-}Prior\ graph) - M^{w_t}_{s\to t,r}(h_t;no\text{-}message\ comparator).
\]

## Reproducibility and reporting

The public package provides the production preprocessing and model code, compact source tables, figure/table registries, executable audit entry points and tests. Raw HMP2/IBDMDB and PRISM files are acquired from their originating repositories and are not redistributed. Reproducibility is presented in three tiers: local code validation, evidence-table audit and optional full rerun after public-input acquisition. Reporting follows TRIPOD, STROBE and PROBAST principles [43–45]. Model interpretation and multiplicity are secondary descriptive analyses [46,47].


