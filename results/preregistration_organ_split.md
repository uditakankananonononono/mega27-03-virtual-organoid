# Pre-registration: out-of-sample test of the organ fidelity split (written before any new data was downloaded)
Date: 2026-09-25 06:23 IST. Commit precedes download of any ArrayExpress / Expression Atlas data.
Hypothesis H1: on organoid bulk RNA-seq profiles NOT in the GEO scan, liver, kidney and lung organoids rank their own GTEx tissue worse than colon/small-intestine and brain organoids.
Scoring: identical to SALL_PARONLY (rank correlation over HPA parenchymal-only genes, 54 GTEx tissues); outcome = rank of labelled tissue (1 = best); one profile per accession (mean log-CPM over all samples, organoid samples only when annotated).
Test: one-sided Mann-Whitney U, unfaithful-organ ranks > faithful-organ ranks, alpha 0.05. Secondary: top-1 rate per group.
Falsified if the one-sided p >= 0.05 or the median rank of the unfaithful group is not higher.
Sources: EBI BioStudies/ArrayExpress processed files and EBI Expression Atlas processed counts; human only; any accession already in the GEO scan (E-GEOD mirrors of scanned GSEs) is excluded.
