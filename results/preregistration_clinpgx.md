# Pre-registration: pharmacogene deficit in liver organoids (ClinPGx / PharmGKB)
Honest note: the ClinPGx genes.tsv file (release 2026-09-05) was downloaded minutes before this commit; only its header was read. It has not been joined to any deficit data.
Gene set: ClinPGx genes with "Is VIP" = Yes OR "Has CPIC Dosing Guideline" = Yes (core pharmacogenes), matched by Ensembl ID.
Deficit vectors: per-organ consensus deficits delta_g = z(GTEx organ) - z(organoid) from src/deficit_bodymap.deficits (GTEx reference).
H1 (applied relevance for drug testing): in liver organoids, core pharmacogenes have higher deficit than other genes (one-sided Mann-Whitney p < 0.01).
H2 (specificity): the liver pharmacogene median deficit exceeds that of each other organ's pharmacogene median deficit (permutation of organ labels on the per-organ standardized deficit, p < 0.05).
Fail = reported negative.
