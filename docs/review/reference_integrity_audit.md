# Reference integrity review

## A. Input Match Check

The clean manuscript sections and the 47-entry reference list were available. Source PDFs were not bundled with the manuscript package, so this is a claim-to-reference alignment screen rather than a source-text certification.

## B. Review Scope Determination

Paragraph-level audit of Introduction, Methods, Results, Discussion and Supplementary text, with a citation-index coverage check.

## C. Main Reference Integrity Findings

- The 47 references are numbered consecutively and all 47 are cited at least once in the clean manuscript.
- Background citations are separated from method citations: IBD multi-omics and cohort context use [1–12], compositional and batch-aware methods use [13–22], microbial profiling uses [23–26], model baselines and graph methods use [27–37], and uncertainty/reporting guidance uses [38–47].
- Claims about graph topology are explicitly limited to representation diagnostics; the text does not cite a graph paper as evidence of a biological interaction.
- The framework cites domain-adversarial learning [37] only to state that this adaptation layer was not fitted because the target domain remained locked.

## D. Major Integrity Risks

No major citation mismatch was detected by the available text-level audit. Direct source checking remains required before submission for any claim that is upgraded from association to mechanism or clinical utility.

## E. Moderate and Minor Integrity Issues

- Several references are method or reporting standards rather than IBD-specific studies; this is appropriate for the framework sections but should not be presented as biological evidence.
- The paper uses candidate-feature language and does not claim that the cited microbiome studies validate the graph-ranked taxa.
- DOI and bibliographic metadata should receive one final Crossref/PubMed check at submission time.

## F. Correction Priority Plan

1. Keep the current bounded wording in Introduction and Discussion.
2. Verify the 2024 differential-abundance benchmark and SVA DOI metadata against the publisher records when the final reference manager is populated.
3. Preserve the current one-to-one numbering when converting the DOCX to the Bioinformatics submission template.

## G. Integrity Logic Explanation

The central claim is methodological: directional transferability is an estimand defined by a frozen input contract and source-only fitting. The citations support the motivation, compositional constraints, model families and reporting standards; the performance numbers are supported by the current experiments and source-data registry rather than by external references. This separation prevents external papers from being used as substitute evidence for the present cross-cohort results.

## H. What Additional Information Would Improve Accuracy

Publisher PDFs or a finalized BibTeX/EndNote library would permit sentence-level verification of every DOI, page range and method-specific claim. The current DOCX is internally consistent but is not a substitute for that final bibliographic check.
