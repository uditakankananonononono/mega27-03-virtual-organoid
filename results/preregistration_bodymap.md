# Pre-registration: reference robustness - Illumina Body Map (Expression Atlas E-MTAB-513) vs GTEx (committed before download)
Question: are organoid deficit rankings an artefact of the GTEx reference?
Method: recompute per-organ consensus deficit vectors (as src/deficit_decoupler.py) with the E-MTAB-513 baseline TPMs (liver, kidney, lung, brain, colon) in place of GTEx, on shared genes.
H: per-organ Spearman rho(deficit_GTEx, deficit_BodyMap) >= 0.5 in all 5 organs, and top-100 deficit overlap exceeds chance (hypergeometric p < 0.001) in all 5.
Secondary: the Open Targets diagonal-dominance test (a73bd26) re-run with BodyMap top-100 lists still has permutation p < 0.05.
Any organ failing = reported as reference-sensitive.
