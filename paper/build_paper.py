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

P.h("5. Negative results (kept by design)")
for t in ["The initial power-law (alpha > 1) interpretation was withdrawn after large-organoid checks contradicted it.",
          "Size-adjusted and size-filtered readouts do not improve donor-level modulator discrimination.",
          "Replication on the public FIS time series is impossible: it is well-level.",
          "Our segmentation is below the published state of the art."]:
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
P.table(["tool", "role"], [["Python 3.10", "language"], ["NumPy", "arrays, bootstrap"], ["pandas", "data"], ["SciPy", "statistics, labelling"],
    ["scikit-image", "watershed"], ["PyTorch", "U-Net"], ["Pillow", "image IO"], ["matplotlib", "figures"], ["pytest", "tests"],
    ["python-docx", "paper"], ["git / GitHub", "version control"]], "Tools table (honest count: 11). Target of 40 not reached.")
P.h("References")
for r in ["Lefferts JW, et al. OrgaSegment: deep-learning based organoid segmentation to quantify CFTR dependent fluid secretion. Commun Biol 2024. https://www.nature.com/articles/s42003-024-05966-4",
          "Boj SF, Vonk AM, et al. Forskolin-induced swelling in intestinal organoids: an in vitro assay for assessing drug response in cystic fibrosis patients. J Vis Exp 2017. https://pmc.ncbi.nlm.nih.gov/articles/PMC5408767/",
          "Botelho HM. FIS_analysis. https://github.com/hmbotelho/FIS_analysis",
          "Measuring cystic fibrosis drug responses in organoids derived from 2D differentiated nasal epithelia. bioRxiv 2021. https://www.biorxiv.org/content/10.1101/2021.07.20.453105v2"]:
    P.p(r)
P.save("paper/mega27-03-virtual-organoid-paper.docx")
print("ok", P.eq, "eq", P.tab, "tab", P.fig, "fig")
