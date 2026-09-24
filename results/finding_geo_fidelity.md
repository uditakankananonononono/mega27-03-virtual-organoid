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

## Label audit and clean-subset re-test (results/geo_label_audit.csv, results/geo_fidelity_clean.json)
- The title-derived organ is confirmed by the GSE summary in 144/152 series. The clean subset (confirmed, single organ, "organoid" in summary) has 114 series; 113 have cached profiles.
- Tissue recovery on the clean subset: top-1 30% / 33% / 27% (S3000 / SALL / CENT).
- Fibroblast best match on the clean subset: 23% [15, 31] / 26% [18, 34] / 36% [27, 45] (bootstrap 95% CI). The one-of-54-columns chance level is 1.9%. The attractor survives the label audit.
- Alternative explanation, not yet excluded: GTEx has only two in-vitro references (cultured fibroblasts, EBV lymphocytes), so a generic culture/proliferation signature could drive the match rather than a mesenchymal identity. Next test: remove proliferation/ECM gene sets, or compare with non-organoid 2D cultures of the same organ.

## Tool: `vorganoid fidelity` (src/vorganoid/fidelity.py, tests/test_fidelity.py)
Real-data run on GSE278954 raw counts (results/tool_fidelity_GSE278954.json): intended tissue Liver ranks 50/54; fibroblast rank 20.
Known issue found while building the tool: HTSeq summary rows ("__no_feature" etc.) were included in CPM denominators in the scan
pipeline (src/geo_fidelity.py). The tool excludes them. The effect on rank correlations is a per-sample scale factor, expected to be small; a rescan is pending.

## HTSeq-row fix and culture-confound test (results/geo_fidelity_methods.json, results/geo_fidelity_clean.json)
- HTSeq "__" rows removed and profiles re-cached. Top-1 recovery rates are unchanged (27.2% / 33.1% / 25.2%), so the bug had no material effect.
- Removed 853 Hallmark proliferation/MYC/EMT genes (MSigDB v2023.2: E2F, G2M, mitotic spindle, MYC v1/v2, EMT).
  Clean-subset fibroblast best match: 25% [17, 33] (SALL_NOCULT) and 42% [34, 51] (CENT_NOCULT). The attractor is NOT explained by these gene sets.
- Remaining alternative, not excluded: cell-type purity. GTEx tissues contain blood, immune, vascular and stromal cells, while cultured fibroblasts
  and organoids are purified cultures. The match may reflect missing non-parenchymal signal rather than mesenchymal identity. A deconvolution-
  or epithelium-only-gene test is the next falsification step. Status: named candidate, not a claimed discovery.

## Purity test with Human Protein Atlas single-cell types (results/geo_fidelity_clean.json, results/geo_mcnemar.json)
- Gene sets from the HPA single-cell type enhancement API: 7,595 genes enhanced in any immune, blood, vascular or stromal type (NONPAR, removed), and 9,699 genes enhanced only in parenchymal types (PARONLY, kept alone).
- Pre-stated prediction (paper v3): with epithelium/parenchyma-only genes, the fibroblast best-match share in the clean subset stays above 15%.
  Result: SALL_PARONLY 28% [20, 36] and CENT_PARONLY 22% [15, 29]. The prediction held, but CENT's lower bound sits exactly at 15%. The purity explanation is weakened, not excluded.
- Method result: removing non-parenchymal genes raises organ-of-origin top-1 recovery on all 151 series. SALL goes 33% -> 42% (McNemar exact p = 0.0044) and CENT 25% -> 38% with PARONLY (p = 0.00031).
  Removing the Hallmark culture genes does not change recovery (p = 1). Interpretation: organoid-vs-tissue fidelity scores are depressed by stromal/immune genes that organoids lack by design. A purity-restricted score is fairer.
- Tool: `vorganoid fidelity --purity HPA_TSV` restricts scoring to HPA parenchymal-only genes. Real run on GSE278954 (results/tool_fidelity_GSE278954_purity.json): Liver rank improves 50 -> 36 of 54, and the fibroblast match becomes top-1 in this series.

## What drives the attractor? g:Profiler enrichment (results/attractor_genes.csv, results/attractor_gprofiler.csv)
- 29 clean-subset series whose SALL best match is fibroblasts. Top 300 genes by d_g = mean_s z(p_sg)(z(G_fib) - z(G_organ)), tested against all scored genes with g:SCS correction.
- Two components:
  (1) Mitotic cell cycle, top term p_adj = 7.8e-56 (66 genes; TOP2A, MKI67, CCNB1, CDK1, BIRC5).
  (2) Immune/complement genes that the tissue has and the organoid lacks: adaptive immune response p_adj = 9.5e-15; complement cascade 5.9e-10 (C1QA, LY86, CTSG).
- Interpretation: the attractor combines culture proliferation with missing non-parenchymal signal. Removing either gene family alone (Hallmark culture sets; HPA non-parenchymal) leaves it in place, so neither alone explains it. It stays a named candidate with a mechanism hypothesis, not a claimed discovery.

## Protein-coding-only robustness (HGNC locus groups; results/geo_fidelity_clean.json, results/geo_mcnemar.json)
- Restricted to 19,254 HGNC protein-coding genes (_PC). The fibroblast best-match share on the clean subset is 44% [35, 53] (SALL_PC) and 28% [20, 36] (CENT_PC). The attractor is not a non-coding/pseudogene artefact.
- Recovery changes: SALL 33% -> 28% (McNemar p = 0.23), CENT 25% -> 30% (p = 0.016).
- Ensembl REST (about 50 s per 500 IDs) and BioMart (timeout) were tried and abandoned for this lookup. HGNC was used instead. Neither Ensembl route is counted as a tool.

## Sample-level audit with GEOparse (results/geo_sample_meta.csv, results/geo_fidelity_strict.json)
GEOparse pulled GSM metadata for all 114 clean series. Only 45 series have every sample annotated as organoid (frac_organoid=1); 59 have no organoid-annotated sample (e.g. GSE343459: fibroblasts isolated from skin organoids; GSE287925: LNCaP 2D cells), despite organoid keywords at series level.
On the strict 45, the "culture-fibroblast attractor" shrinks: fibroblast-best fraction 0.09-0.29 across 11 method variants (vs 0.28-0.62 in the no-organoid-sample series; Fisher p<0.05 in 10/11 variants). Top-1 tissue recovery on the strict subset is 0.31-0.56.
Verdict (negative/partial): a large part of the attractor signal came from series-level label contamination, not from organoids. The attractor remains a candidate for a minority (~10-25%) of strict organoid series. Not claimed as a discovery.

## Independent annotation of attractor driver genes: Enrichr + Reactome AnalysisService (src/attractor_enrichr_reactome.py)
Top-200 genes by d (results/attractor_genes.csv). Note d>0 arises two ways: sample-high and fibroblast-high (proliferation), or sample-low and tissue-high (missing tissue programs).
- Enrichr MSigDB Hallmark: G2-M checkpoint q=2.5e-33, E2F targets q=3.2e-32, mitotic spindle q=6.2e-17. CellMarker_2024: cycling/MKI67+ progenitor signatures (q=4e-60). PanglaoDB: immune types (gamma-delta T q=1.6e-10, plasma cells q=4e-6).
- Reactome AnalysisService: classical antibody-mediated complement activation and Cell Cycle, Mitotic (both FDR 1.4e-14).
This confirms the g:Profiler result with two independent services: the attractor is proliferation plus absent immune/complement programs. Contrast set (ranks 801-1000, still d>0) is weak (Hallmark q>=0.07; Reactome phagocytosis FDR 1.2e-4). Given the GEOparse strict-subset result above, this describes what drives the fibroblast match, not an organoid biology discovery.
