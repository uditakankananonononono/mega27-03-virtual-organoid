# Which organoids are faithful? Strict-subset organ analysis (src/strict_organ_deficits.py, src/strict_organ_parenchymal.py)

Strict series = all GSMs annotated organoid (45; results/geo_sample_meta.csv). Recovery = GTEx rank of labelled organ, SALL_PARONLY (results/strict_organ_recovery.csv).

- Colon (19 series) top-1 0.84, brain cortex (9) 0.78; liver (6) 0.17, kidney cortex (4) 0.00, lung (4) 0.00.
- Colon+brain 23/28 top-1 vs liver+kidney+lung 1/14; Fisher p=3.9e-6 (results/strict_organ_summary.json). Organ groups were chosen after seeing the table, so this p is descriptive, not confirmatory.
- What the unfaithful organoids lack (top-100 consensus deficit genes, STRING enrichment + Enrichr PanglaoDB; results/strict_organ_deficit_enrichment.csv):
  - Liver: hepatocyte program (PanglaoDB Hepatocytes q=1.7e-31; STRING organic acid metabolism FDR 4.6e-12, PPAR signaling 3e-5).
  - Kidney: proximal tubule program (PanglaoDB q=5.7e-22; ALDOB, SLC17A3) plus complement C1q.
  - Lung: mostly stroma (fibroblasts, pericytes) and immune genes, with NAPSA (AT2) among top deficits.
  - Faithful colon/brain organoids lack mainly stroma, vascular and immune genes (colon: peritubular myoid / smooth muscle; brain: oligodendrocytes, microglia, complement C1q).
- Mechanism test (negative/inconclusive): HPA parenchymal-only share of deficit genes, unfaithful vs faithful organs, OR=1.45, p=0.061 (results/strict_organ_parenchymal.json). Liver (0.68) and kidney (0.55) deficits are parenchymal; lung (0.15) is not.
Named candidate, not claimed: "metabolic-parenchyma gap" - liver and kidney organoids fail to match their tissue because they lack mature metabolic epithelium (hepatocyte, proximal tubule), whereas intestinal and cortical organoids match despite lacking stroma. Falsifiable: in a new set of strict liver/kidney organoid series, restoring or scoring only hepatocyte/proximal-tubule genes should explain most of the rank gap. Caveats: small n per organ; organ confounded with lab and protocol (PSC-derived vs adult stem cell), not tested.

## Protocol confound (src/strict_protocol.py, results/strict_protocol.csv, results/strict_protocol.json)
Each strict series was classified from series + GSM text (GEOparse) by keyword counts as PSC-derived (16), adult/tissue-derived (22) or unclear (7). Classification is regex-based and not manually validated.
- Derivation alone: top-1 0.44 (PSC) vs 0.73 (adult), Fisher p=0.099.
- Organ effect within derivation: PSC faithful-organ 7/11 vs unfaithful 0/5 (p=0.034); adult 16/16 vs 0/5 (p=4.9e-5). CMH (Haldane 0.5) pooled OR=46, p=5.8e-5.
- Unfaithful organs fail regardless of derivation (0/5 PSC, 0/5 adult).
Verdict: the organ split is not explained by PSC vs adult derivation. Brain organoids are all PSC-derived, so derivation cannot be separated from organ for brain. The "metabolic-parenchyma gap" stays a named candidate (post hoc grouping, small n).

## UniProt secretome check (src/strict_deficit_uniprot.py, results/strict_deficit_uniprot.csv)
Share of deficit genes whose reviewed human UniProt entry has a signal peptide or 'Secreted' location, vs 2,000 random expressed genes (16.6%):
brain 0.44, colon 0.41, kidney 0.36, liver 0.30, lung 0.57; all enriched (OR 2.2-6.8, q<=0.001).
Verdict (negative for the organ split): every organoid type lacks secreted/extracellular proteins relative to tissue, faithful or not. The secretome deficit is universal and does not explain why liver/kidney/lung fail. Note: gene_exact queries can return extra entries (n_mapped slightly >100 for some organs).

## Pre-registered out-of-sample test on ArrayExpress (results/preregistration_organ_split.md committed 23fa32e before download; src/ae_scan.py, src/ae_replication.py, results/ae_replication.csv/.json)
BioStudies/ArrayExpress search -> 18 curated human non-cancer organoid E-MTAB accessions with processed files; 11 scored (7 unusable: mouse gene IDs E-MTAB-9181/14831/11273; only Illumina probe IDs E-MTAB-4591; malformed table E-MTAB-12548; single-cell/ATAC only E-MTAB-10876/15659). E-MTAB-17066 used its first 3 per-sample featureCounts files (download budget).
- Faithful organs (intestine 5, brain 2): median rank 1, top-1 6/7 (brain E-MTAB-10037 ranked 26, best match cultured fibroblasts).
- Unfaithful organs: median rank 19.5, top-1 1/4. Pre-registered one-sided Mann-Whitney p=0.031 -> H1 formally supported.
Honest reading: support is weak. The unfaithful side has no lung; both liver datasets are cholangiocyte (biliary) organoids on Illumina microarrays (ranks 39, 37; best match stomach), so the liver signal mixes cell type and platform. Both kidney datasets match kidney well (tubuloids rank 1; PSC kidney organoids rank 2), which contradicts the kidney part of the candidate. The "metabolic-parenchyma gap" remains a candidate; the kidney component is not replicated.

## GSEApy prerank GSEA of the full deficit ranking (src/strict_gsea.py, results/strict_gsea_hallmark.csv, results/strict_gsea_summary.json)
Hallmark v2023.2, 500 permutations, per organ consensus delta over all genes.
- Liver top deficits: xenobiotic (NES 2.29), coagulation, fatty acid, OXPHOS, bile acid metabolism (all FDR ~0). Kidney: complement/coagulation/xenobiotic. Lung: interferon, EMT (stroma/immune).
- Mean NES over 6 metabolic hallmarks: liver 1.96, kidney 1.95, brain 1.65, colon 0.74, lung 0.46. Liver+kidney vs colon+brain one-sided MWU over NES values p=0.004 (hallmarks overlap, so values are not independent; descriptive).
Verdict: supports the metabolic-deficit description for liver and kidney organoids, but it is not specific: brain organoids also lack OXPHOS/fatty-acid programs (NES 2.63/1.77) yet match their tissue. So a metabolic deficit alone does not decide fidelity; the candidate stays unproven.
