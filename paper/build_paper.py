"""Build the item-3 paper from committed data/results. Run from repo root."""
import json, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, "paper")
from paperkit import Paper

# per-donor attenuation table, recomputed from raw per-organoid data and saved as a results file
m = pd.read_csv("data/raw/orgasegment/dis_merged_A0.csv"); m = m[m["forskolin_concentration_µM"] > 0]; m["ls"] = np.log(m.swelling)
q = np.quantile(m.A0, [0.25, 0.75]); m["sz"] = np.where(m.A0 < q[0], "small", np.where(m.A0 > q[1], "large", "mid"))
rows = []
for don, g in m.groupby("donor"):
    a = g[g.condition == "VX445_VX661_VX770"]; c = g[g.condition == "DMSO"]
    if len(a) < 30 or len(c) < 30: continue
    es = a[a.sz == "small"].ls.mean() - c[c.sz == "small"].ls.mean(); el = a[a.sz == "large"].ls.mean() - c[c.sz == "large"].ls.mean()
    rows.append([don, g.genotype.iloc[0], len(a), len(c), round(es, 3), round(el, 3), round(el - es, 3), round(a.ls.mean() - c.ls.mean(), 3)])
R = pd.DataFrame(rows, columns=["donor", "genotype", "n_trikafta", "n_dmso", "eff_small", "eff_large", "attenuation", "overall"])
R.to_csv("results/per_donor_attenuation.csv", index=False)
rho, pv = spearmanr(R.overall, R.attenuation)
S1 = json.load(open("results/seg_eval.json")); S2 = json.load(open("results/seg_eval_tuned_256.json")); S3 = json.load(open("results/seg_eval_512.json")); S4 = json.load(open("results/seg_eval_tuned_512.json"))

P = Paper("Small-Organoid Attenuation of CFTR-Modulator Swelling in Patient-Derived Intestinal Organoids: "
          "a Single-Organoid Re-Analysis with a Biophysical Swelling Twin",
          "MEGA-PROGRAM-27, Item 3 - Udita Phookan (program owner); computational work by an AI research agent. Draft of 24 September 2026.")
P.h("Abstract")
P.p("The forskolin-induced swelling (FIS) assay on patient-derived intestinal organoids is used to predict which people with cystic fibrosis "
    "(CF) respond to CFTR modulators. Standard pipelines sum organoid area per well, so any dependence of the response on organoid size is "
    "averaged away. We re-analysed the public single-organoid OrgaSegment data (17 CF donors, per-organoid area before and after "
    "forskolin) with a biophysical swelling model and a DMSO-controlled, donor-matched design.")
P.p(f"Finding (candidate). Modulator-induced swelling is attenuated in the smallest organoids and saturates above a size threshold. For "
    f"elexacaftor/tezacaftor/ivacaftor (Trikafta) the drug-minus-DMSO log-swelling effect rises from 0.336 [0.262, 0.405] in the smallest "
    f"size octile to a plateau near 0.55-0.64; top minus bottom 0.251 [0.184, 0.323] (donor bootstrap). Per donor, the effect in the "
    f"largest size quartile exceeds that in the smallest in {int((R.attenuation > 0).sum())} of {len(R)} donors, and the attenuation does "
    f"not track overall response (Spearman rho = {rho:.2f}, p = {pv:.2f}). Geometry predicts the opposite sign: under a uniform "
    "surface flux, small organoids should swell more. Area noise also biases the slope negative. So the effect is unlikely to be a "
    "measurement artefact.")
P.p(f"Negatives. Size-adjusted and size-filtered readouts do not improve per-donor Trikafta-vs-DMSO separation, so the finding does not "
    f"change theratyping calls at this assay's well counts. Replication on the public FIS time series is impossible because it is "
    f"well-level. Our best U-Net segmentation (512 px, tuned) reaches mAP@0.5 = {S4['eval_mAP50']:.3f} +/- {S4['eval_sd']:.3f} on the OrgaSegment eval split "
    f"(published 0.76 +/- 0.12): below the state of the art, although the gap is within one standard deviation on 12 images.")
P.p("Fidelity. Scoring 151 organoid GEO series against GTEx, cultured fibroblasts are the most common best match. A GEOparse sample-level audit shows this "
    "'culture-fibroblast attractor' comes mostly from non-organoid samples in organoid-titled series: on 45 strictly organoid series it falls to 9-29%. "
    "We retract it as a general organoid property.")

P.h("1. Introduction")
P.p("CFTR moves chloride and bicarbonate across the apical membrane of epithelial cells; water follows, and in a closed organoid the lumen "
    "swells. Forskolin raises cAMP and opens CFTR, so the swelling of rectal organoids measures residual and drug-rescued CFTR function in "
    "each patient. The assay is used to support access to modulators for people with rare CFTR genotypes.")
P.p("Organoids in a well vary in size by more than an order of magnitude. If the swelling response depends on size, the well-level readout "
    "depends on the size mix in each well, which varies with passage, seeding density and culture time. We asked whether the single-organoid "
    "response depends on starting size, and whether correcting for it improves drug-response calls.")

P.h("2. Biophysical swelling twin")
P.p("Treat an organoid as a near-spherical shell with volume V and projected area A, V proportional to A^(3/2). Let net secretion scale as a "
    "power of volume:")
P.equation("dV/dt = k V^alpha")
P.p("For alpha != 1 the solution over an assay window t is")
P.equation("V_1^(1-alpha) = V_0^(1-alpha) + (1 - alpha) k t")
P.p("and the observed area fold change is s = A_1/A_0 = (V_1/V_0)^(2/3). For small responses,")
P.equation("log s ≈ (2/3) k t V_0^(alpha - 1)")
P.p("so the slope of log log s on log A_0 carries the sign of alpha - 1:")
P.equation("d log(log s) / d log A_0 = (3/2)(alpha - 1)")
P.p("Surface-limited secretion (constant flux per unit apical area) gives alpha = 2/3 and a negative slope: small organoids swell more in "
    "relative terms. Size-invariant secretion gives alpha = 1 and zero slope. The within-well regression we fit is")
P.equation("log s_ow = a_w + b ( log A0_ow - mean_w log A0 ) + e_ow")
P.p("with a well intercept a_w that absorbs donor, plate and well effects. Measurement noise in A0 enters both sides with opposite sign "
    "(s = A1/A0), which biases b downward:")
P.equation("E[ b-hat ] = b - Var(u) / ( Var(log A0) + Var(u) ),   u = noise in log A0")
P.p("so a positive observed b is conservative. The donor-matched drug effect per size bin k is")
P.equation("Delta_k = mean_{o in drug, bin k} log s_o - mean_{o in DMSO, bin k} log s_o")
P.p("and the attenuation statistic is Delta_top - Delta_bottom, with 95% intervals from a donor-level bootstrap:")
P.equation("CI_95 = [Q_0.025, Q_0.975] of { Delta*_top - Delta*_bottom }, donors resampled with replacement")
P.p("Segmentation quality is scored by average precision at IoU 0.5, as in OrgaSegment:")
P.equation("AP_50 = TP / (TP + FP + FN),   match if IoU(P, G) = |P ∩ G| / |P ∪ G| ≥ 0.5")

P.h("3. Data")
P.table(["#", "dataset", "source / accession", "content", "use"], [
    [1, "OrgaSegment DIS single-organoid measurements", "Zenodo 10610438 (Lefferts et al. 2024)", "per-organoid A0/A1, 17 donors, 4 conditions", "size-dependence analysis"],
    [2, "OrgaSegment FIS database", "Zenodo 10610438", "well-level area, 7 time points, 868 rows", "replication attempt (not possible)"],
    [3, "OrgaSegment annotated images", "Zenodo 10278229", "train/val/eval images with instance masks", "U-Net segmentation benchmark"],
], "Dataset manifest (distinct, accession-level). Honest count: 3. The program target of 120+ was not reached.")

P.h("4. Results")
P.h("4.1 Within-well size slopes", 2)
P.table(["condition", "organoids", "wells", "slope b", "95% CI", "donors b>0"], [
    ["DMSO", 3810, 94, 0.011, "[-0.002, 0.028]", "12/17"], ["VX770", 1496, 30, 0.226, "[0.185, 0.269]", "5/5"],
    ["VX661+VX770", 4493, 112, 0.042, "[0.028, 0.056]", "12/16"], ["VX445+VX661+VX770", 4984, 121, 0.116, "[0.090, 0.142]", "14/17"],
], "Within-well slope of log fold change on centred log A0 (forskolin > 0). Well-cluster bootstrap, 500 resamples.")
P.p("The slope is near zero without CFTR rescue and positive whenever a modulator is present. The sign is opposite to the surface-limited "
    "prediction. A power-law reading (alpha > 1) was rejected: restricted to organoids above 2,000 px the slope turns negative (Trikafta "
    "-0.117), and binned means rise then plateau. Simulated segmentation noise on a size-invariant truth gave slopes of only -0.002 to -0.055.")
P.h("4.2 Small-organoid attenuation", 2)
P.figure("results/figures/fig_size_attenuation.png", "Drug-minus-DMSO mean log swelling per starting-size octile, donor-matched, with donor-bootstrap 95% CI.")
P.table(["modulator", "donors", "smallest octile effect", "plateau", "top minus bottom [95% CI]"], [
    ["Trikafta", 17, "0.336 [0.262, 0.405]", "0.55-0.64", "0.251 [0.184, 0.323]"],
    ["VX661+VX770", 16, 0.136, "0.28-0.31", "0.169 [0.081, 0.259]"],
    ["VX770", 5, -0.153, "0.31-0.42", "0.462 [0.335, 0.572]"],
], "Size-bin drug effects.")
P.h("4.3 Per-donor consistency", 2)
P.table(["donor", "genotype", "eff small", "eff large", "attenuation", "overall"],
        R.sort_values("overall")[["donor", "genotype", "eff_small", "eff_large", "attenuation", "overall"]].values.tolist(),
        "Trikafta minus DMSO log-swelling effect in the smallest and largest starting-size quartiles (results/per_donor_attenuation.csv).")
P.p("Association between attenuation and overall response uses Spearman's rank correlation:")
P.equation("rho = 1 - 6 sum_d (r_d - r'_d)^2 / ( n (n^2 - 1) ),   r, r' = ranks of attenuation and overall effect over n donors")
P.p(f"Attenuation is positive in {int((R.attenuation > 0).sum())} of {len(R)} donors and unrelated to overall response "
    f"(rho = {rho:.2f}, p = {pv:.2f}). Two of three donors with two nonsense alleles (X/X) and one F508del/R117H donor are exceptions.")
P.h("4.4 Does size correction help theratyping? No", 2)
P.p("Locked before running: per-donor standardised Trikafta-vs-DMSO separation with raw, size-adjusted and size-band well readouts. Median "
    "raw 4.73, adjusted 4.35, band 4.39; adjusted better in 8 of 17 donors. Excluding organoids below 1,069 px (a post-hoc threshold): median "
    "4.73 -> 4.75, better in 9 of 17, Wilcoxon p = 0.68. The size effect is real but does not change calls.")
P.h("4.5 Segmentation benchmark", 2)
P.table(["model", "eval mAP@0.5", "sd", "source file"], [
    ["U-Net v1, 256 px, default post-processing", round(S1["mAP50"], 3), round(S1["sd"], 3), "results/seg_eval.json"],
    ["U-Net v1 + val-tuned post-processing + flip TTA", round(S2["eval_mAP50"], 3), round(S2["eval_sd"], 3), "results/seg_eval_tuned_256.json"],
    ["U-Net v2, 512 px, 60 epochs, default post-processing", round(S3["mAP50"], 3), round(S3["sd"], 3), "results/seg_eval_512.json"],
    ["U-Net v2 + val-tuned post-processing + flip TTA", round(S4["eval_mAP50"], 3), round(S4["eval_sd"], 3), "results/seg_eval_tuned_512.json"],
    ["OrgaSegment (published, Mask R-CNN)", 0.76, 0.12, "Lefferts et al. 2024"],
], "Instance segmentation on the OrgaSegment eval split (12 images).")
P.p("Caveat: the validation images used for tuning were also in the training set, so the tuning signal is optimistic; the eval split was "
    "scored once. The 512-px model (resumed from checkpoint after an out-of-memory kill at epoch 18) raises the tuned score from "
    f"{S2['eval_mAP50']:.3f} to {S4['eval_mAP50']:.3f}, still below the published mean. With 12 eval images and sd near 0.13, "
    "neither a win nor a loss against OrgaSegment is statistically resolved; we report it as below state of the art.")

P.h("4.6 Organoid-to-tissue fidelity across 151 GEO series", 2)
GM = json.load(open("results/geo_fidelity_methods.json")); GC = json.load(open("results/geo_fidelity_clean.json"))
P.p(f"To test organoid fidelity at scale, we searched GEO through NCBI E-utilities for human organoid bulk RNA-seq series (2,838 hits; 1,545 with "
    f"text count matrices). For {GM['n_series']} series with a usable gene-level matrix and a title naming a single organ, each sample was converted to "
    "log counts per million and averaged. The series profile was then correlated with every GTEx v8 tissue median profile (54 tissues).")
P.equation("CPM_gs = 10^6 x_gs / sum_g x_gs,   p_g = (1/S) sum_s log(1 + CPM_gs)")
P.equation("r_t^cent = <p - mu_p - (m - mu_m), G_t - m - mean(G_t - m)> / (||.|| ||.||),   m_g = mean_t G_gt")
P.p("The permutation p-value for the top-1 recovery rate reshuffles organ labels across series, B = 10,000 times:")
P.equation("p_perm = (1 + #{b : top1_b >= top1_obs}) / (1 + B)")
P.p("Uncertainty on shares uses a series-level percentile bootstrap:")
P.equation("CI_95 = [q_0.025, q_0.975] of { f(D*_b) },  D*_b drawn with replacement from the series set, b = 1..2000")
rowsm = [[m, GM["variants"][m]["top1"], GM["variants"][m]["perm_top1_mean"], GM["variants"][m]["perm_top1_p"], GC[m]["top1"], GC[m]["fibroblast_best_frac"],
          f"[{GC[m]['fibroblast_best_ci95'][0]:.2f}, {GC[m]['fibroblast_best_ci95'][1]:.2f}]"] for m in GM["variants"]]
P.table(["method", "top-1 (all)", "null top-1", "p_perm", "top-1 (clean)", "fibroblast best (clean)", "95% CI"],
        [[r[0]] + [round(x, 3) if isinstance(x, float) and x > 0.001 else (f"{x:.0e}" if isinstance(x, float) else x) for x in r[1:]] for r in rowsm],
        "Tissue recovery and fibroblast best-match share. S3000: Spearman on the 3,000 most tissue-variable genes; SALL: all genes; CENT: GTEx-centred Pearson; "
        f"_NOCULT: {GM['n_culture_genes_removed']} MSigDB Hallmark proliferation/MYC/EMT genes removed. Clean subset: {GC['n_clean']} series whose summary confirms a single organ and "
        "mentions organoids (113 with profiles). Files: results/geo_fidelity_methods.json, results/geo_fidelity_clean.json.")
P.figure("results/figures/fig_geo_fidelity.png", "Per-organ share of clean-subset series whose organ of origin (blue) or GTEx cultured fibroblasts (orange) is the top-1 match. Organs with n >= 3.")
P.p("Organ of origin is recovered four to five times above chance. Brain organoids are the most faithful, and lung, breast and pancreas the least. "
    "Across every method, the most frequent best match is GTEx cultured fibroblasts, which we name the culture-fibroblast attractor. It survives the label audit and the removal of "
    "proliferation and EMT genes. It is not yet a discovery: GTEx tissues contain immune, vascular and stromal cells, while organoids and cultured "
    "fibroblasts are purified cultures, so cell-type purity is an unexcluded explanation. Falsifiable prediction: restricted to epithelium-specific "
    "genes, the fibroblast best-match share in the clean subset will stay above 15%. The tool vorganoid fidelity scores any count matrix this way.")

P.h("4.7 Purity-corrected fidelity (Human Protein Atlas)", 2)
MC = json.load(open("results/geo_mcnemar.json"))
P.p(f"GTEx tissues contain immune, blood, vascular and stromal cells that organoids lack by design. Using the Human Protein Atlas single-cell type "
    f"atlas, we removed {GM['n_nonparenchymal_genes_removed']:,} genes enhanced in any such cell type (NONPAR), or kept only the "
    f"{GM['n_parenchymal_only_genes']:,} genes enhanced solely in parenchymal types (PARONLY). Paired top-1 recovery per series is compared with the exact McNemar test:")
P.equation("p_McN = 2 * sum_{k=0}^{min(b,c)} C(b+c, k) 0.5^(b+c),   b = base-only hits, c = variant-only hits")
P.table(["comparison", "top-1 base", "top-1 variant", "b", "c", "exact p"],
        [[k.replace("_vs_", " vs "), round(v["top1_base"], 3), round(v["top1_variant"], 3), v["table"][0][1], v["table"][1][0], f"{v['p_exact']:.2g}"] for k, v in MC.items()],
        "Organ-of-origin recovery on 151 series, base method vs gene-restricted variant (results/geo_mcnemar.json).")
P.p(f"Purity restriction raises recovery (all-gene Spearman {MC['SALL_vs_SALL_NONPAR']['top1_base']:.0%} -> {MC['SALL_vs_SALL_NONPAR']['top1_variant']:.0%}, "
    f"p = {MC['SALL_vs_SALL_NONPAR']['p_exact']:.2g}; centred Pearson {MC['CENT_vs_CENT_PARONLY']['top1_base']:.0%} -> {MC['CENT_vs_CENT_PARONLY']['top1_variant']:.0%}, "
    f"p = {MC['CENT_vs_CENT_PARONLY']['p_exact']:.2g}). Removing culture genes does not change it. Standard whole-transcriptome fidelity scores are therefore biased low for organoids. "
    f"The prediction stated before this test (fibroblast best-match share above 15% on parenchymal-only genes) held: "
    f"{GC['SALL_PARONLY']['fibroblast_best_frac']:.0%} [{GC['SALL_PARONLY']['fibroblast_best_ci95'][0]:.2f}, {GC['SALL_PARONLY']['fibroblast_best_ci95'][1]:.2f}] and "
    f"{GC['CENT_PARONLY']['fibroblast_best_frac']:.0%} [{GC['CENT_PARONLY']['fibroblast_best_ci95'][0]:.2f}, {GC['CENT_PARONLY']['fibroblast_best_ci95'][1]:.2f}]. "
    "The lower bound of the second interval touches the threshold, so purity is weakened but not excluded as the explanation. "
    "vorganoid fidelity --purity applies this correction to any count matrix.")

P.h("4.8 What drives the attractor, and is it robust to gene class?", 2)
AG = _pd_gp = pd.read_csv("results/attractor_gprofiler.csv"); AG = AG[AG.term_size < 2000]
P.p("For the clean-subset series whose best match is cultured fibroblasts, each gene is scored by how much it pulls the profile toward "
    "fibroblasts and away from the organ of origin, averaged over series:")
P.equation("d_g = (1/|S|) sum_{s in S} z(p_sg) ( z(G_g,fib) - z(G_g,organ(s)) ),   z = rank-based standard score within a profile")
P.p("The top 300 genes were tested with g:Profiler (GO:BP, Reactome, KEGG; g:SCS correction; background = all scored genes).")
sel = pd.concat([AG.head(4), AG[AG.name.str.contains("immune|complement", case=False)].head(3)])
P.table(["source", "term", "adjusted p", "genes"], [[r.source, r.name, f"{r.p_adj:.1e}", int(r.intersection)] for r in sel.itertuples()],
        "Enrichment of attractor-driving genes (results/attractor_gprofiler.csv; terms with more than 2,000 genes omitted).")
P.p("Two components appear: mitotic cell-cycle genes (TOP2A, MKI67, CDK1) that cultures express and tissues do not, and immune and complement "
    "genes that tissues contain and organoids lack. Removing either family alone leaves the attractor in place (Sections 4.6-4.7). "
    f"Restricting to {GM['n_protein_coding_genes']:,} HGNC protein-coding genes also leaves it in place (fibroblast best match "
    f"{GC['SALL_PC']['fibroblast_best_frac']:.0%} [{GC['SALL_PC']['fibroblast_best_ci95'][0]:.2f}, {GC['SALL_PC']['fibroblast_best_ci95'][1]:.2f}] and "
    f"{GC['CENT_PC']['fibroblast_best_frac']:.0%} [{GC['CENT_PC']['fibroblast_best_ci95'][0]:.2f}, {GC['CENT_PC']['fibroblast_best_ci95'][1]:.2f}]), so it is not a non-coding artefact. "
    "Our working hypothesis is that the attractor is the joint signature of proliferation and purity. It remains a candidate: a direct test needs "
    "organoid and matched primary tissue from the same donors, profiled together.")


P.h("4.9 Sample-level audit and independent annotation: the attractor is mostly label contamination", 2)
ST = json.load(open("results/geo_fidelity_strict.json")); SM = pd.read_csv("results/geo_sample_meta.csv")
EN = pd.read_csv("results/attractor_enrichr.csv"); RE = pd.read_csv("results/attractor_reactome.csv")
meths = [k for k, v in ST.items() if isinstance(v, dict)]
fs = [ST[m]["strict"]["fibroblast_best_frac"] for m in meths]; fn = [ST[m]["no_organoid"]["fibroblast_best_frac"] for m in meths]
nsig = sum(ST[m]["fisher_p_strict_vs_none"] < 0.05 for m in meths)
P.p(f"Series-level text says a study is about organoids; it does not say every profiled sample is an organoid. We therefore pulled the sample (GSM) "
    f"metadata of all {len(SM)} clean series with GEOparse and classified each sample as organoid or not from its title, source and characteristics. "
    f"Only {ST['n_strict']} series have every sample annotated as organoid; {ST['n_no_organoid_samples']} have no organoid-annotated sample at all "
    "(for example GSE343459 profiles fibroblasts isolated from skin organoids, and GSE287925 profiles LNCaP cells in 2D). We compared the fibroblast-best "
    "fraction between the strict and the no-organoid series with Fisher's exact test:")
P.equation("p = sum_{x : P(x) <= P(a)} C(K, x) C(N - K, n - x) / C(N, n)")
P.table(["method", "strict n", "strict top-1", "strict fibroblast-best [95% CI]", "no-organoid fibroblast-best", "Fisher p"],
        [[m, ST[m]["strict"]["n"], f"{ST[m]['strict']['top1']:.2f}", f"{ST[m]['strict']['fibroblast_best_frac']:.2f} [{ST[m]['strict']['fibroblast_ci95'][0]:.2f}, {ST[m]['strict']['fibroblast_ci95'][1]:.2f}]",
          f"{ST[m]['no_organoid']['fibroblast_best_frac']:.2f}", f"{ST[m]['fisher_p_strict_vs_none']:.2g}"] for m in meths],
        "Attractor on strict organoid series vs series without organoid samples (results/geo_fidelity_strict.json).")
P.p(f"On the strict subset the attractor shrinks to {min(fs):.2f}-{max(fs):.2f} of series (vs {min(fn):.2f}-{max(fn):.2f}); the difference is significant in "
    f"{nsig} of {len(meths)} method variants. Much of the attractor therefore came from non-organoid samples inside organoid-titled series, not from organoids. "
    "We retract it as a general property of organoids; it persists only in a minority of strict series.")
up = EN[EN.direction == "top200"]
h = up[up.library == "MSigDB_Hallmark_2020"].iloc[0]; c = up[up.library == "CellMarker_2024"].iloc[0]; r0 = RE[RE.direction == "top200"].iloc[0]
P.p(f"Independent annotation of the top 200 driver genes with Enrichr and the Reactome AnalysisService reproduces the g:Profiler picture: Hallmark {h.term} "
    f"(q = {h.q_bh:.1e}), CellMarker '{c.term}' (q = {c.q_bh:.0e}), and Reactome '{r0['name']}' (FDR = {r0.fdr:.1e}) together with mitotic cell cycle. "
    "The driver genes describe what makes a profile look like cultured fibroblasts (proliferation, missing immune and complement programs); they are not an organoid discovery.")

P.h("4.10 Which organoids are faithful? Organ, deficits and protocol on the strict subset", 2)
REC = pd.read_csv("results/strict_organ_recovery.csv"); SUM = json.load(open("results/strict_organ_summary.json"))
ENR = pd.read_csv("results/strict_organ_deficit_enrichment.csv"); PAR = json.load(open("results/strict_organ_parenchymal.json")); PRO = json.load(open("results/strict_protocol.json"))
P.table(["labelled organ", "strict series", "top-1 recovery", "median rank of own tissue"],
        [[r.organ, int(r.n), f"{r.top1:.2f}", f"{r.median_rank:g}"] for r in REC.sort_values("n", ascending=False).itertuples()],
        "Organ recovery on the 45 strictly organoid series (SALL_PARONLY; results/strict_organ_recovery.csv).")
fb = SUM["fisher_colon_brain_vs_liver_kidney_lung"]
P.p(f"Colon and cortical organoids recover their own tissue ({fb['table'][0][0]} of {sum(fb['table'][0])} top-1); liver, kidney and lung organoids almost never do "
    f"({fb['table'][1][0]} of {sum(fb['table'][1])}; Fisher p = {fb['p']:.1e}). The grouping was chosen after seeing the table, so this p is descriptive. To see what the "
    "unfaithful organoids lack we scored, for each series, the deficit of each gene relative to its labelled tissue and averaged within organ:")
P.equation("delta_g = (1/|S_o|) sum_{s in S_o} [ z(G_g,o) - z(p_sg) ]")
def top(o, src):
    x = ENR[(ENR.organ == o) & (ENR.source == src)].sort_values("fdr")
    return f"{x.iloc[0].term} (q = {x.iloc[0].fdr:.1e})" if len(x) else "none"
P.p(f"The 100 largest deficits were annotated with STRING functional enrichment and Enrichr PanglaoDB. Liver organoids lack the hepatocyte program "
    f"({top('Liver', 'Enrichr:PanglaoDB')}; STRING {top('Liver', 'STRING:Process')}); kidney organoids lack proximal tubule "
    f"({top('Kidney - Cortex', 'Enrichr:PanglaoDB')}); lung organoids lack mainly stroma ({top('Lung', 'Enrichr:PanglaoDB')}). Faithful colon and cortical organoids "
    f"lack mostly stromal, vascular and immune cells ({top('Colon - Transverse', 'Enrichr:PanglaoDB')}; {top('Brain - Cortex', 'Enrichr:PanglaoDB')}).")
fp = PAR["fisher_unfaithful_vs_faithful"]
P.p(f"A direct mechanism test is inconclusive: the share of HPA parenchymal-only genes among the deficits is only modestly higher in unfaithful organs "
    f"(OR = {fp['odds_ratio']:.2f}, p = {fp['p']:.3f}), because lung deficits are stromal (share {PAR['Lung']['par_share_of_annotated']:.2f}) while liver "
    f"({PAR['Liver']['par_share_of_annotated']:.2f}) and kidney ({PAR['Kidney - Cortex']['par_share_of_annotated']:.2f}) deficits are parenchymal.")
cm = PRO["cmh_organ_effect_given_derivation"]; ps = PRO["organ_effect_within_PSC"]; ad = PRO["organ_effect_within_adult"]; fd = PRO["fisher_psc_vs_adult_top1"]
P.p(f"Protocol. Classifying series from GEO text as PSC-derived or adult-tissue-derived (keyword rules, not manually validated), derivation alone is weakly "
    f"associated with recovery (Fisher p = {fd['p']:.2f}). The organ effect holds within each stratum (PSC {ps['table'][0][0]}/{sum(ps['table'][0])} vs "
    f"{ps['table'][1][0]}/{sum(ps['table'][1])}, p = {ps['fisher_p']:.3f}; adult {ad['table'][0][0]}/{sum(ad['table'][0])} vs {ad['table'][1][0]}/{sum(ad['table'][1])}, "
    f"p = {ad['fisher_p']:.1e}) and in a Cochran-Mantel-Haenszel test:")
P.equation("OR_MH = sum_k (a_k d_k / n_k) / sum_k (b_k c_k / n_k)")
P.p(f"(pooled OR = {cm['pooled_OR_haldane']:.0f} with a 0.5 continuity correction, p = {cm['p']:.1e}). We name the candidate the metabolic-parenchyma gap: liver and "
    "kidney organoids miss their tissue because they lack mature metabolic epithelium, while intestinal and cortical organoids match despite missing stroma. "
    "It is falsified if, in new strictly organoid liver or kidney series, hepatocyte and proximal-tubule genes do not account for most of the rank gap. "
    "Brain organoids are all PSC-derived, so organ and protocol cannot be separated there; per-organ n is small.")

UP = pd.read_csv("results/strict_deficit_uniprot.csv")
P.p(f"Secretome check (UniProt). The share of deficit genes whose reviewed UniProt entry has a signal peptide or a 'Secreted' location is "
    f"{UP.share.min():.2f}-{UP.share.max():.2f} per organ against {UP.bg_share.iloc[0]:.2f} for 2,000 random expressed genes (all q <= {UP.q_bh.max():.3f}). "
    "Every organoid type lacks secreted proteins, faithful or not, so this does not explain the organ split (negative).")
P.h("4.11 Pre-registered out-of-sample test on ArrayExpress", 2)
AE = pd.read_csv("results/ae_replication.csv"); AJ = json.load(open("results/ae_replication.json")); ok = AE[AE.status == "ok"]
P.p("Before downloading any new data we committed a pre-registration (results/preregistration_organ_split.md): on organoid bulk profiles outside the GEO "
    "scan, liver, kidney and lung organoids should rank their own tissue worse than intestinal and brain organoids, tested by a one-sided Mann-Whitney U at 0.05:")
P.equation("U = sum_{i in unfaithful} sum_{j in faithful} [ 1(r_i > r_j) + 0.5 * 1(r_i = r_j) ]")
P.p(f"A BioStudies/ArrayExpress search gave {len(AE)} curated human non-cancer organoid studies with processed files; {len(ok)} could be scored "
    "(the others had mouse gene IDs, probe IDs only, malformed tables, or single-cell data only).")
P.table(["accession", "organ", "rank of own tissue", "best GTEx match"], [[r.accession, r.organ, int(r.rank_match), r.best_tissue] for r in ok.itertuples()],
        "Out-of-sample organoid profiles (results/ae_replication.csv).")
P.p(f"Median rank {AJ['median_rank_unfaithful']:g} vs {AJ['median_rank_faithful']:g}; top-1 {AJ['top1_unfaithful']:.2f} vs {AJ['top1_faithful']:.2f}; one-sided p = {AJ['mwu_one_sided_p']:.3f}. "
    "The pre-registered test passes, but the support is weak: there is no lung dataset, both liver datasets are cholangiocyte organoids on microarrays "
    "(so cell type and platform are confounded), and both kidney datasets match kidney well (ranks 1 and 2), which contradicts the kidney part of the candidate. "
    "The metabolic-parenchyma gap remains a candidate; its kidney component did not replicate.")
P.h("5. Negative results (kept by design)")
for t in ["The initial power-law (alpha > 1) interpretation was withdrawn after large-organoid checks contradicted it.",
          "Size-adjusted and size-filtered readouts do not improve donor-level modulator discrimination.",
          "Replication on the public FIS time series is impossible: it is well-level.",
          "Our segmentation is below the published state of the art.",
          "Fidelity: 'Pancreas' as best match (35/151 series with S3000) disappears with other methods, so it is a method artefact.",
          "Fidelity: a 'liver disease' series (GSE278954) ranks Liver 50th of 54, and lung and breast organoids are rarely matched to their organ.",
          "Fidelity: organ labels come from text, and series profiles average all samples, including any non-organoid controls.",
          "Fidelity: the HPA parenchymal share of deficit genes does not clearly separate faithful from unfaithful organs (p = 0.061); lung organoids break the pattern.",
          "Fidelity: all organoid types lack secreted proteins (UniProt), so the secretome does not explain the organ split.",
          "Replication: the kidney part of the organ split failed out of sample (two ArrayExpress kidney datasets rank kidney 1st and 2nd).",
          "Fidelity: the culture-fibroblast attractor is RETRACTED as a general organoid property: on 45 series whose samples are all organoids it shrinks to 9-29% of series (Section 4.9)."]:
    P.p("- " + t)
P.h("6. Discussion")
P.p("Small organoids respond less to CFTR modulators, relative to DMSO, than large ones, and the response saturates above a size threshold. "
    "Because geometry predicts the reverse, the likeliest explanations are biological: small organoids may lack a mature lumen, contain "
    "fewer differentiated secretory cells, or have a closed tight-junction seal later. Falsifiable prediction: in time-lapse single-organoid "
    "data from an independent lab, the Trikafta-minus-DMSO effect in the top starting-size quartile will exceed that in the bottom quartile "
    "by at least 0.1 log units in most donors. Checking this needs a second single-organoid FIS dataset, which we did not find in public form.")
P.p("Practical note. Well-level FIS readouts depend on the size mix. Labs comparing drug response across passages or sites should report "
    "the starting-size distribution.")
P.p("Novelty. We did not find this effect reported for intestinal FIS in the literature we searched (OrgaSegment; the FIS_analysis pipeline, "
    "which sums area per well). Related: nasal 2D-derived organoids show swelling variation concentrated in large structures "
    "(bioRxiv 2021.07.20.453105). This was not a systematic review.")
P.h("7. Tools used")
TL = pd.read_csv("results/tools_ledger.csv")
TL["gate"] = TL.counts_for_gate.astype(str).map({"True": "counts", "False": "infra (excluded)"})
P.table(["tool", "kind", "where used", "gate"], TL[["tool", "kind", "where_used", "gate"]].values.tolist(), f"Tools genuinely used (results/tools_ledger.csv): {len(TL)} entries, {int((TL.gate == 'counts').sum())} counting toward the gate after excluding infrastructure. Target of 40 not reached.")
P.h("References")
for r in ["Lefferts JW, et al. OrgaSegment: deep-learning based organoid segmentation to quantify CFTR dependent fluid secretion. Commun Biol 2024. https://www.nature.com/articles/s42003-024-05966-4",
          "Boj SF, Vonk AM, et al. Forskolin-induced swelling in intestinal organoids: an in vitro assay for assessing drug response in cystic fibrosis patients. J Vis Exp 2017. https://pmc.ncbi.nlm.nih.gov/articles/PMC5408767/",
          "GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 2020. https://gtexportal.org",
          "Karlsson M, et al. A single-cell type transcriptomics map of human tissues. Sci Adv 2021. https://www.proteinatlas.org",
          "Szklarczyk D, et al. The STRING database in 2023. Nucleic Acids Res 2023.",
          "The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Res 2023.",
          "Moreno P, et al. Expression Atlas and ArrayExpress/BioStudies at EMBL-EBI. Nucleic Acids Res 2022.",
          "Kuleshov MV, et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Res 2016.",
          "Milacic M, et al. The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res 2024.",
          "Gumienny R. GEOparse: Python library to access Gene Expression Omnibus. https://github.com/guma44/GEOparse",
          "Kolberg L, et al. g:Profiler - interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update). Nucleic Acids Res 2023. https://biit.cs.ut.ee/gprofiler",
          "Seal RL, et al. Genenames.org: the HGNC resources in 2023. Nucleic Acids Res 2023. https://www.genenames.org",
          "Liberzon A, et al. The Molecular Signatures Database Hallmark gene set collection. Cell Syst 2015. https://www.gsea-msigdb.org",
          "Barrett T, et al. NCBI GEO: archive for functional genomics data sets. Nucleic Acids Res 2013. https://www.ncbi.nlm.nih.gov/geo/",
          "Botelho HM. FIS_analysis. https://github.com/hmbotelho/FIS_analysis",
          "Measuring cystic fibrosis drug responses in organoids derived from 2D differentiated nasal epithelia. bioRxiv 2021. https://www.biorxiv.org/content/10.1101/2021.07.20.453105v2"]:
    P.p(r)
P.save("paper/mega27-03-virtual-organoid-paper.docx")
print("ok", P.eq, "eq", P.tab, "tab", P.fig, "fig")
