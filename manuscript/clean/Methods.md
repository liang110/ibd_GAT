# Methods

## Study design and cohorts

We designed the study as a source-only, participant-aware comparison of marginal and relational microbiome representations for directional transfer. The primary cohorts were HMP2/IBDMDB and PRISM [1,2]. The binary endpoint was IBD versus non-IBD control; the analysis was not designed to resolve Crohn disease, ulcerative colitis, disease activity or therapy effects. HMP2/IBDMDB contained 1,638 metagenomic samples from 130 participants and 578 raw microbial features, and PRISM contained 220 samples from 220 participants and 209 raw microbial features. PRISM was retained as two origin strata: discovery (155 samples; 121 IBD and 34 controls) and validation (65 samples; 43 IBD and 22 controls). Participant identifiers, repeated-sample structure, class counts, origin labels and input checksums were frozen before model fitting.

The methodological contribution is a linked set of estimands and evidence gates rather than a new neural layer. Figure 1 gives the complete information-flow contract. The first stage fixes the shared feature universe and the metadata boundary; the second fits transformations and model choices only in source-training participants; the third measures directional transfer and uncertainty; the final gates test leakage, origin, topology and generic capacity. The framework is designed so that a graph branch can fail without invalidating the benchmark, and so that an apparent graph gain can be attributed only to the declared representation and transfer direction.

## Frozen input contract and feature harmonization

Microbial identifiers were normalized with a conservative canonicalizer, and identifiers mapping to the same canonical taxon were summed. Clinical, outcome, medication, antibiotic, steroid, age, faecal-calprotectin and sample/participant identifier fields were removed by one shared denylist before any model branch was built. The cross-cohort universe was the deterministic intersection of the two microbial registries:

\[
C = C_HMP2 ∩ C_PRISM,  |C| = 126.
\]

Here C_HMP2 and C_PRISM are the post-canonicalization microbial feature registries. This universe is label-free and defines the same node and tabular coordinate system in both cohorts. Three HMP2 samples had zero total mass in C and were excluded before transformation; no PRISM sample was excluded. Exact-name metabolite intersection was zero, so the primary cross-cohort estimand is microbiome-only. Multi-omics fields were retained only for within-cohort sensitivity records and never used to support a cross-cohort metabolome claim. The contract, feature order, denylist, input hashes and zero-mass decisions were written to the Gate0 manifest (Fig. 2; Supplementary Table 1).

## Compositional transformation and source-only fitting

Both cohorts were closed within C before model fitting. For sample x with non-negative abundance coordinates, closure was defined as:

\[
C(x)ⱼ = xⱼ / Σₖ₌₁ᵖ xₖ.
\]

The pseudocount δ was one half of the smallest positive source-training value. The centred log-ratio transformation was:

\[
CLR(x)ⱼ = log(C(x)ⱼ + δ) − (1/p)Σₖ₌₁ᵖ log(C(x)ₖ + δ).
\]

The order of operations was fixed: canonicalization and intersection, zero-mass exclusion, closure, source-training prevalence filtering, source-training pseudocount, source-training scaling and source-training mutual-information selection. Prevalence below 0.10 was evaluated in source-training participants only. At most 100 features with the largest source-training mutual information with the endpoint were retained for a model split; the 126-feature value therefore denotes the frozen universe, not the number of features presented to every fitted model.

For source cohort s and target cohort t, all fitted quantities were represented by a source object Tₛ; [train] denotes the source-training participants. The deployment boundary was:

\[
Tₛ = Fit(Xₛᵗʳᵃⁱⁿ, yₛᵗʳᵃⁱⁿ),  ŷₜ = f(Tₛ(Xₜ)).
\]

No target label, target prevalence, target ranking, target threshold or target fusion weight entered Tₛ. We did not add a domain-adversarial adaptation layer; the target remained a locked evaluation domain [37]. The retained feature list, feature-order hash, fitted transform, split identifiers and predictions were archived for every source split. This explicit operator is the main protection against the lookahead concern raised by Reviewers 1 and 3.

## Directional transfer estimands and participant-aware evaluation

We evaluated HMP2→PRISM and PRISM→HMP2 transfer with target samples locked after source fitting. A separate PRISM origin estimand evaluated discovery→validation and validation→discovery transfer, rather than pooling the two origins. Within-cohort source splits and all repeated leakage controls used participant-grouped stratification whenever identifiers were available. The primary confirmation and paired uncertainty analyses used ten source repeats; topology controls used five source repeats.

For metric M, model f, reference g, source s and target t, the primary relational contrast was:

\[
ΔMₛ→ₜ(f,g) = Mₛ→ₜ(f) − Mₛ→ₜ(g).
\]

This contrast is intentionally directional. It does not assert a population-wide interaction or a universal model ranking. AUROC, average precision and Brier score were reported with source-repeat variability and participant-cluster uncertainty. Target labels were used only after predictions were frozen for evaluation.

## Representation branches and graph construction

The tabular branch contained prevalence, L2 logistic regression, elastic net, random forest and XGBoost references, all receiving the same Tₛ(X) matrix and the same denylist. Elastic net used C=1, l₁ ratio 0.5 and class balancing; random forest used 300 trees, minimum leaf size 2 and class balancing. XGBoost was retained as a source-validation reference, not as evidence of a new architecture. These references address the reviewers’ request for like-for-like baselines and make the practical engineering trade-off visible: a simpler tabular model may be preferable when the graph does not transport.

The relational branch was the No-Prior GATv2 model. The legacy whole-cohort biological prior was excluded from the primary analysis because it was estimated outside split boundaries and could encode target information. Each sample was represented as a graph whose nodes were the fold-local selected microbial features and whose scalar node values were the transformed abundances. Association matrices were computed only from source-training profiles. For Pearson and Spearman topologies, each node retained its eight strongest non-self absolute associations:

\[
Eₖ = {(i,j) : j ∈ Topₖ(|Rᵢⱼ|), j ≠ i},  k = 8.
\]

R is the source-training association matrix. No-edge and identity graphs were controls. The rewired null preserved each node’s outgoing edge count while randomly reassigning destinations. The edge rule is a representation choice, not proof of a biological interaction; Pearson and Spearman were compared because compositional dependence can make linear association sensitive to scale and closure [13–21].

The model used two GATv2Conv layers with hidden width 16, two heads in the first layer and one in the second, dropout 0.1, global mean pooling and a two-layer classifier head. Parameters were optimized with Adam (learning rate 10⁻³, weight decay 10⁻³) for 40 epochs with batch size 32. GATv2 denotes the graph-attention architecture [33,34], not a manuscript version label. The implementation details and topology hashes were stored with each run.

## Audit gates and capacity control

The framework contains four validity gates. The leakage gate compares naive sample splitting with participant-grouped splitting and tests a planted clinical leakage validator. The known-truth gate evaluates null, batch-only, signal-only, reversing-batch and topology-shift semi-synthetic scenarios. The origin gate keeps PRISM discovery and validation separate. The topology gate reports the five topology controls, edge-set overlap and the fixed-k construction diagnostic.

The participant-overlap optimism gap was defined as:

\[
Ω_leak = AUROC_naive − AUROC_grouped.
\]

For edge sets Eₐ and E_b, overlap was defined as:

\[
J(Eₐ,E_b) = |Eₐ ∩ E_b| / |Eₐ ∪ E_b|.
\]

A no-message-passing multilayer perceptron used the same Tₛ(X) input boundary and approximately matched hidden capacity. The capacity-adjusted relational contrast was:

\[
Δ_cap = M(No-Prior\ graph) − M(capacity\ comparator).
\]

Because fixed top-k assigns eight outgoing edges to every retained node, the out-degree vector d has zero variance:

\[
Var(d) = 0  ⇒  conventional\ corr(d,d′)\ is\ undefined.
\]

The value 1.0 recorded by an audit artifact is retained only as an engineering placeholder and is not interpreted as ecological conservation.

## Participant-cluster bootstrap and uncertainty

For each source split, bootstrap draws resampled participants rather than individual samples and retained all samples from a participant within a draw. For metric M and bootstrap draw b, the paired difference was:

\[
ΔM⁽ᵇ⁾ = M⁽ᵇ⁾_GAT − M⁽ᵇ⁾_reference.
\]

The 2.5th and 97.5th percentiles defined a split-specific 95% interval. Ten-split summaries are arithmetic means of split-specific point estimates and interval endpoints; they are not pooled cohort-level confidence intervals. This separation distinguishes source-training variability from target participant resampling, as requested by Reviewer 3. AUROC interpretation was accompanied by average precision because the cohorts were IBD-skewed [39–42].

## Reproducibility and reporting

Each run recorded the Git revision, input and output SHA-256 hashes, feature-contract hash, preprocessing manifest, split identifiers, random seed, topology hash, model parameters, predictions, metrics and Slurm job record. Figure and table source registries map every manuscript object to one frozen artifact. Reporting follows the logic of TRIPOD and STROBE for transparent prediction and observational reporting [43–45]; model interpretation and multiplicity are treated as secondary, descriptive analyses [46,47]. The final release package will add the environment lockfile, license, acquisition instructions and checkpoint policy before submission.
