# GEO organoid-fidelity scan: nulls and method checks (exploratory)
Source: results/geo_fidelity_methods.json, results/geo_fidelity_methods_ranks.csv (151 of 152 scored GSEs; GSE261999 profile failed the 150 MB decompression guard).
- The matched GTEx tissue is the top-1 match in 27.2% (S3000), 33.1% (all-gene Spearman) and 25.2% (GTEx-centred Pearson) of series.
  Label-permutation null (10,000 draws): 6.0%, 8.0% and 5.7%. p < 1e-4 for all three. Analytic chance is 6.2%.
- Robust across all three methods: "Cells - Cultured fibroblasts" is the most frequent best match (34, 38 and 52 of 151 series).
- Method-sensitive: "Pancreas" as best match (35 / not in top 5 / 18). Treated as a method artefact of the variance-selected gene set, not a finding.
- Per organ (top-1 rate): brain 0.75 (n=20); colon/intestine 0.40 (48); liver 0.17 (23); lung 0.09 (11); breast 0.00 (7); pancreas 0.00 (8).
Caveats: organ labels come from series titles (not yet audited against summaries). Series profiles average all samples, which may include
non-organoid controls, 2D cultures or cell lines. Counts vs TPM inputs are mixed.
Candidate (named, falsifiable, not yet claimed): "culture-fibroblast attractor". A large share of organoid transcriptomes are closer to
GTEx cultured fibroblasts than to their organ of origin. It should hold after restricting to organoid-only samples and audited labels.
