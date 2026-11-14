# Known-truth acceptance audit

`per_seed_acceptance.csv` is the compact truth-panel ledger. Positive scenarios use the declared AND rule: target AUROC >= 0.70 **and** recovery >= 3 of 6 disease-truth features. Negative scenarios use the declared AUROC interval 0.35–0.65. The acceptance rate is the number of passing seeds divided by 10; no mean is used to decide a seed.
