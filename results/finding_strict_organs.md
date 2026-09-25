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

## Open Targets disease relevance (pre-registered a73bd26; src/deficit_opentargets.py, results/deficit_opentargets.json)
Top-500 Open Targets targets per organ disorder (MONDO liver/kidney/lung/brain/colonic). 5x5 overlap of deficit lists vs disease sets.
- Pre-registered primary: diagonal dominance D = 2.6 (null mean -0.01), permutation p = 0.0005 -> PASS. Organoid deficits are organ-specifically disease-relevant as a set.
- Per organ (Fisher, own vs other lists): lung OR 4.5 p=0.017; colon OR 4.5 p=0.017; kidney OR 2.0 p=0.13; brain p=0.45; liver OR 1.0 p=0.62.
- Negative for the metabolic-parenchyma idea: liver, the core metabolic organ, shows no own-disease enrichment (2/100 hits). The pass is driven by lung and colon. Counts are small (2-8 hits per cell).

## decoupler TF / pathway activity (pre-registered a72ccc4; src/deficit_decoupler.py, results/deficit_decoupler*.{json,csv})
decoupler 2.1.4 ulm with CollecTRI (TFs) and PROGENy top-500 (pathways) from OmniPath, on per-organ deficit vectors (4,261 genes shared across all strict series).
- Pre-registered positive control FAILED: only 2/12 organ master TFs in the top 10% of deficient TFs (binomial p=0.34). Liver CEBPA (rank 6/397) and HNF1A (24) are deficient; colon CDX1/CDX2, lung NKX2-1/FOXA2, kidney PAX8/HNF4A are not. Brain master TFs are not covered on the shared gene set.
- Exploratory (not pre-registered): every organ shows the same PROGENy pattern - PI3K higher in organoids than tissue (score -4.8 to -8.7) and JAK-STAT/p53 lower in organoids (+2.7 to +5.8). This is a shared culture signature, which fits the non-specific GSEA result, not an organ-specific gap.

## Reference robustness: Illumina Body Map (pre-registered 3d898d6; src/deficit_bodymap.py, results/deficit_bodymap.json)
Deficits recomputed with Expression Atlas E-MTAB-513 baseline TPMs instead of GTEx.
- Per-organ Spearman GTEx vs Body Map deficits: liver 0.72, lung 0.72, kidney 0.67, brain 0.61, colon 0.57; top-100 overlap 52-75 (hypergeometric p < 1e-89) -> PASS in all 5.
- Open Targets diagonal-dominance re-run with Body Map lists: D=2.1, permutation p=0.0026 -> PASS.
Caveat: both deficit vectors share the same organoid profiles, so part of the agreement comes from the organoid side; this checks reference sensitivity only.

## Pharmacogenes in liver organoids (pre-registered afcb277; src/deficit_clinpgx.py, results/deficit_clinpgx.json)
ClinPGx (formerly PharmGKB) genes.tsv, release 2026-09-05.
- Pre-registered gene set (VIP or CPIC) is UNUSABLE: the "Is VIP" column is "Yes" for all 25,041 genes in this release. With that set H1 trivially fails (p=0.97); this run is kept as a negative/data-defect record.
- Documented deviation: CPIC dosing-guideline genes only (33 with Ensembl IDs; 21 on the liver deficit vector).
  - H1: liver organoids under-express CPIC pharmacogenes vs other genes (median deficit 0.32 vs 0.03, Mann-Whitney p=0.0004) -> PASS under deviation. Top: CYP4F2, CYP2D6, SLCO1B1, VKORC1, CYP2C9, TPMT, UGT1A1, CYP2C19.
  - H2 (liver-specific): FAIL (p=0.23; only 6 CPIC genes are on all 5 organ vectors; brain and kidney medians are similar).
Practical reading: liver organoids in these public series are short on CPIC-actionable drug-metabolism genes, which matters for drug testing, but the evidence for liver specificity is absent and the set was changed after pre-registration.

## gnomAD population constraint of missing organoid genes (pre-registered 58a422c; src/deficit_gnomad.py, results/deficit_gnomad.json)
Official gnomAD v2.1.1 gene-level pLoF constraint (LOEUF, `oe_lof_upper`); ambiguous duplicate gene symbols excluded. Liver 99/100 mapped, brain 94/100. Median LOEUF liver 1.256 vs brain 1.107; pre-registered one-sided Mann-Whitney U=4900, p=0.263 -> **FAIL**. The descriptive kidney-vs-brain comparison also fails (p=0.228). Lower LOEUF means more LoF intolerance; this test gives no clear evidence that missing liver genes are less LoF constrained than missing brain genes. Brain's 1.107 median also cautions against assuming its missing genes are mostly essential developmental factors; the top deficits include blood and immune signals. Organs and their deficits were chosen post hoc before this test; population haploinsufficiency is not a causal test of organoid fidelity. The v2.1.1 constraint scale must not be compared to v4 cutoffs.

## Ensembl coding-composition check (pre-registered 8ea5c50; src/deficit_ensembl_biotype.py, results/deficit_ensembl_biotype.json)
Ensembl REST `POST /lookup/id` mapped all 500 fixed top-deficit IDs on GRCh38. Protein-coding fractions: brain 100/100, colon 100/100, liver 99/100, kidney 93/100, lung 92/100. Pre-registered liver > brain one-sided Fisher exact p=1.0 -> **FAIL**. This is a low-information check: almost all top deficits in both groups are coding and the input GTEx expression universe already favors annotated genes. It does not support a distinct coding-gene composition explanation for the liver gap; it does not test whether particular mature hepatocyte genes matter.

## Ontology label audit (pre-registered fecf4f7; OLS4, results/ols_anatomy_audit.json)
OLS4 returned unique exact-label ontology entries for cortex of kidney (UBERON:0001225) vs kidney (UBERON:0002113), and cholangiocyte (CL:1000488) vs hepatocyte (CL:0000182). The literal query `proximal tubule cell` returned no exact result; this alias is unresolved under the pre-registered query (exploratory rephrasing finds `epithelial cell of proximal tubule`, CL:0002306, but it is not counted as pre-registered resolution). H1 passes only the distinct-label audit, not a part-of proof or expression analysis. This confirms that the liver cholangiocyte-to-whole-liver and kidney organoid-to-adult-cortex comparisons should not be read as matched-cell-type validations. It does not change the failed kidney component or prove the metabolic-parenchyma candidate. Individual query URLs and response labels are recorded in results/ols_anatomy.json.

## Independent single-cell kidney-organoid stress test (pre-registered 3330723; GSE108291; src/sc_kidney_markers.py)
Two deposited human kidney organoid 10x matrices, not the mouse kidney reference, were analyzed as one accession with Scanpy sparse QC. The main `org` matrix has 2,211,840 barcodes including empty droplets, of which 9,638 pass >=500 UMIs and >=200 genes; `org4` has 1,421/1,421 passing. In `org`, ALDOB occurs in 6 QC cells, SLC17A3 in 4, with zero coexpressing both; `org4` has 2 and 1 respectively, again zero both. The pre-registered >=1% two-marker H1 fails in both runs; the <20% H2 passes. **This is a negative for this strict two-marker coexpression operationalization, not evidence that proximal tubule cells are absent.** Sparse scRNA-seq dropout, developmental state, protocol, and marker choice can all explain zeros. In addition, QC thresholds were fixed after pre-registration, before reading marker coexpression, and are explicitly disclosed in results/sc_kidney_markers.json. This new cohort does not overturn the earlier whole-tissue kidney fidelity success; it suggests bulk match and mature proximal tubule coexpression answer different questions. A single accession GSE108291 is added to the ledger, not two datasets. Original HCA project: https://explore.data.humancellatlas.org/projects/7b947aa2-43a7-4082-afff-222a3e3a4635 .

## HPO clinical phenotype annotation (pre-registered 16771b5; results/deficit_hpo.json)
Official HPO 2026-09-01 gene-to-phenotype file maps only 35/100 liver and 28/100 brain top deficit genes to **any** direct HPO phenotype, below the pre-registered 50/100 coverage minimum. Therefore the primary liver-vs-brain test is **UNINTERPRETABLE** and no Fisher p or clinical claim is reported. The fixed direct metabolism/homeostasis term HP:0001939 occurs in 3 liver genes (AGXT, CYP2D6, DMGDH), 1 kidney (HMGCS2), and none in brain/colon/lung. Without ancestor propagation and with missing clinical annotations, the zeroes are not biological absences. This resource helps trace annotations, not diagnose patients.

## InterPro fixed-target audit (pre-registered cdb1051, invalid primary test)
Twelve proteins were mapped to unique reviewed human UniProt entries and InterPro domain entries (results/deficit_interpro.json). The pre-registration incorrectly called six liver drug-metabolism genes "top-100 deficits". **Only CYP2D6 actually belongs to the fixed liver top 100**, whereas all six chosen brain genes are in its top 100. The purported H1 contrast is therefore INVALID; it cannot count as a biological success. The annotation still confirms CYP2D6/2C9/2C19 have P450 domains, UGT1A1 a UDP-glucuronosyltransferase family, SLCO1B1 an organic anion transporter domain, and VKORC1 a vitamin K epoxide reductase domain. But the target choice does not fairly compare the groups. This mistake is kept visibly rather than silently editing the pre-registration or claiming the naive 5/6 keyword pass.
