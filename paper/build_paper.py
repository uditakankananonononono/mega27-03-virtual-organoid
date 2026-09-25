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
          "an Exploratory Single-Organoid Re-Analysis and Negative Validation Tests",
          "MEGA-PROGRAM-27, Item 3 - Udita Phookan (program owner); computational work by an AI research agent. Working draft of 25 September 2026.")
P.h("Abstract")
P.p("Novelty status. This working project tests whether baseline organoid size moderates CFTR-modulator swelling, using a donor-matched analysis and a swelling-geometry model. This is a testable candidate, NOT an established new biological discovery. Its stronger pre-registered donor x plate x dose check failed; the tested segmentation pipeline does not match the published benchmark; no diagnostic tool has been validated. The contribution is a transparent, reproducible analysis with negative controls, not an award or clinical claim.")
P.p("The forskolin-induced swelling (FIS) assay on patient-derived intestinal organoids is used to predict which people with cystic fibrosis "
    "(CF) respond to CFTR modulators. Standard pipelines sum organoid area per well, so any dependence of the response on organoid size is "
    "averaged away. We re-analysed the public single-organoid OrgaSegment data (17 CF donors, per-organoid area before and after "
    "forskolin) with a biophysical swelling model and a DMSO-controlled, donor-matched design.")
P.p(f"Finding (candidate). In this accession, estimated modulator-induced swelling is lower in the smallest tracked organoids and approaches a plateau across larger size bins; this is a post hoc descriptive curve, not a validated biological threshold. For "
    f"elexacaftor/tezacaftor/ivacaftor (Trikafta) the drug-minus-DMSO log-swelling effect rises from 0.336 [0.262, 0.405] in the smallest "
    f"size octile to a plateau near 0.55-0.64; top minus bottom 0.251 [0.184, 0.323] (donor bootstrap). Per donor, the effect in the "
    f"largest size quartile exceeds that in the smallest in {int((R.attenuation > 0).sum())} of {len(R)} donors, and the attenuation does "
    f"not track overall response (Spearman rho = {rho:.2f}, p = {pv:.2f}). Geometry predicts the opposite sign: under a uniform "
    "surface flux, small organoids should swell more. Area noise also biases the slope negative. So the effect is unlikely to be a "
    "simple uniform-area-noise artefact under the tested model; track-selection bias and source-specific confounding remain unbounded. A stronger matched-block test fails, so this is not donor-general evidence.")
P.p(f"Negatives. Size-adjusted and size-filtered readouts do not improve per-donor Trikafta-vs-DMSO separation, so the finding does not "
    f"change theratyping calls at this assay's well counts. A pre-registered donor x plate x dose test is 9/12 positive (p=0.073), failing H1. Replication on the public FIS time series is impossible because it is "
    f"well-level. Our best U-Net segmentation (512 px, tuned) reaches mAP@0.5 = {S4['eval_mAP50']:.3f} +/- {S4['eval_sd']:.3f} on the OrgaSegment eval split "
    f"(published 0.76 +/- 0.12): below the reported reference. A separate train-only, validation-only pipeline also fails (0.73867 by the published AP scorer). Neither score establishes a match, superiority, or clinical use.")
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
P.p("For a strictly positive response (s > 1), the idealized model implies the local slope of log(log s) on log A_0 carries the sign of alpha - 1:")
P.equation("d log(log s) / d log A_0 = (3/2)(alpha - 1)")
P.p("Surface-limited secretion (constant flux per unit apical area) gives alpha = 2/3 and a negative slope: small organoids swell more in "
    "relative terms. Size-invariant secretion gives alpha = 1 and zero slope in this idealized model. "
    "We do NOT fit log(log s): it is undefined when s <= 1, as occurs in 1,520 of 14,783 positive-forskolin observations (10.3%). Instead our empirical within-well regression is")
P.equation("log s_ow = a_w + b ( log A0_ow - mean_w log A0 ) + e_ow")
P.p("A positive empirical b has the same directional prediction only for positive responses under the idealized model; it does not estimate alpha and can reflect non-secretion mechanisms, especially for s <= 1. With a well intercept a_w absorbing donor, plate and well effects, measurement noise in A0 enters both sides with opposite sign "
    "(s = A1/A0), which biases b downward under independent zero-mean log-area errors u at baseline "
    "and endpoint, with latent log baseline size x_true:")
P.equation("E[ b-hat ] ≈ (b Var(x_true) - Var(u)) / (Var(x_true) + Var(u)),   x = log A0 = x_true + u")
P.p("Thus a positive observed b is conservative only under this independent-error model; correlated or size-dependent segmentation and track selection can reverse the bias. The donor-matched drug effect per size bin k is")
P.equation("Delta_k = mean_{o in drug, bin k} log s_o - mean_{o in DMSO, bin k} log s_o")
P.p("and the attenuation statistic is Delta_top - Delta_bottom, with 95% intervals from a donor-level bootstrap:")
P.equation("CI_95 = [Q_0.025, Q_0.975] of { Delta*_top - Delta*_bottom }, donors resampled with replacement")
P.p("Segmentation quality is scored by average precision at IoU 0.5, as in OrgaSegment:")
P.equation("AP_50 = TP / (TP + FP + FN),   match if IoU(P, G) = |P ∩ G| / |P union G| >= 0.5")

P.h("3. Data")
P.table(["#", "dataset", "source / accession", "content", "use"], [
    [1, "OrgaSegment DIS single-organoid measurements", "Zenodo 10610438 (Lefferts et al. 2024)", "per-organoid A0/A1, 17 donors, 4 conditions", "size-dependence analysis"],
    [2, "OrgaSegment FIS database", "Zenodo 10610438", "well-level area, 7 time points, 868 rows", "replication attempt (not possible)"],
    [3, "OrgaSegment annotated images", "Zenodo 10278229", "train/val/eval images with instance masks", "U-Net segmentation benchmark"],
    [4, "Drevinek intestinal FIS cohort", "Zenodo 4771466", "20 patient IDs, 54 plates, well-level 0-60 minute area series", "triple-vs-double treatment ranking"],
    [5, "Borek-Dohalska intestinal FIS cohort", "Zenodo 15754800", "12 patient IDs, 33 plates, well-level areas; possible overlap with 4771466", "negative VTI-vs-ETI ranking test"],
], "Dataset manifest. Four distinct Zenodo accessions are shown: the DIS and FIS files are two assay tables under the SAME accession 10610438, not two datasets. The ledger has 168 primary-used accession entries across imaging, FIS and transcriptomic analyses; separate Drevinek-group releases may overlap patients and are not individual-object replication.")

P.h("4. Results")
P.h("4.1 Within-well size slopes", 2)
P.table(["condition", "organoids", "wells", "slope b", "95% CI", "donors b>0"], [
    ["DMSO", 3810, 94, 0.011, "[-0.002, 0.028]", "12/17"], ["VX770", 1496, 30, 0.226, "[0.185, 0.269]", "5/5"],
    ["VX661+VX770", 4493, 112, 0.042, "[0.028, 0.056]", "12/16"], ["VX445+VX661+VX770", 4984, 121, 0.116, "[0.090, 0.142]", "14/17"],
], "Within-well slope of log fold change on centred log A0 (forskolin > 0). Well-cluster bootstrap, 500 resamples.")
P.p("The empirical slope is near zero in DMSO and positive in each modulator group. Its sign contrasts with a simple surface-limited "
    "secretion prediction only for strictly positive responses; it does not estimate alpha for the full data, which include shrinking objects. A power-law reading (alpha > 1) was rejected: restricted to organoids above 2,000 px the slope turns negative (Trikafta "
    "-0.117), and binned means rise then plateau. A specific simulation of independent, size-invariant segmentation noise on a flat true response gave slopes of -0.002 to -0.055; size-dependent segmentation and differential track retention were not ruled out.")
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
WR = json.load(open("results/withincohort_replication.json"))
P.h("4.3a Internal experiment-date replication, not a new cohort", 2)
P.p(f"We fixed the earlier small (<722 px) and large (>=1400 px) cutoffs, then split each donor's first "
    f"experiment date from later dates in the SAME OrgaSegment DIS accession. Requiring at least five organoids "
    f"per DMSO/Trikafta by size cell, {WR['portions']['later']['n_positive']}/{WR['portions']['later']['n_donors']} "
    f"later-date donors have a positive donor-matched large-minus-small Trikafta effect: median "
    f"{WR['portions']['later']['median_attenuation']:.3f} log-swelling units, exact one-sided sign "
    f"p={WR['portions']['later']['sign_p_one_sided']:.6f}. On first-date plates, "
    f"{WR['portions']['earliest']['n_positive']}/{WR['portions']['earliest']['n_donors']} eligible donors are positive "
    f"(median {WR['portions']['earliest']['median_attenuation']:.3f}, p={WR['portions']['earliest']['sign_p_one_sided']:.3f}). "
    "The later-date H1 passes its pre-registration, but thresholds were chosen using the whole accession, donors "
    "overlap, and there is no external individual-organoid cohort. This supports repeatability over experiment "
    "dates, not independent discovery, clinical validity, or the lumen-maturation mechanism "
    "(results/withincohort_replication.json).")

BL = json.load(open("results/blocked_size.json"))
P.h("4.3b Matched plate-dose stress test: fails", 2)
P.p(f"The stronger pre-registered stress test required >=3 organoids per Trikafta/DMSO by small/large "
    f"cell within each donor x plate x forskolin concentration block. Only {BL['n_eligible_blocks']} "
    f"blocks and {BL['n_eligible_donors']} donors qualify; {BL['n_positive_donors']}/"
    f"{BL['n_eligible_donors']} donor median block effects are positive, with exact one-sided "
    f"sign p={BL['sign_p_one_sided']:.3f}. H1 FAILS (registered alpha 0.05). "
    f"{BL['n_excluded_donors']} donors have no eligible block. The pooled block median "
    f"({BL['median_block_effect']:.3f}) is descriptive and not a donor-independent test. "
    "This challenges the pooled later-date 12/12 result: dose/plate control with narrower cells alters "
    "three donor signs. It does not show the size effect is absent, but the evidence does not establish "
    "a new validated discovery (results/blocked_size.json).")

EXT = json.load(open("results/external_fis_demo.json"))
P.h("4.3c External single-plate intestinal FIS assay", 2)
P.p(f"A separate public intestinal organoid demonstration experiment from the FIS_image_analysis repository "
    f"provides object-level tracked areas in {len(EXT['paths'])} CSVs from ONE plate, one class-II CFTR genotype, "
    f"and two wells per treatment at each of eight matched forskolin doses. The registered analysis paired "
    f"{EXT['n_tracked_pairs']:,} unique time-0/time-60-min tracks and compared combined VX-809/VX-770 "
    "with forskolin-only vehicle. Fixed-analysis pooled baseline quartile boundaries were "
    f"{EXT['quartile_cutoffs_area_micronsq']['small_lte']:.0f} and "
    f"{EXT['quartile_cutoffs_area_micronsq']['large_gte']:.0f} square microns, distinct units/cutoffs from OrgaSegment. "
    f"The drug-minus-vehicle large-minus-small log-swelling effect is positive in "
    f"{EXT['n_positive_doses']}/{EXT['n_eligible_doses']} dose strata "
    f"(one-sided sign p={EXT['sign_p_one_sided']:.4f}; H1 passes); the lowest dose is negative. "
    "This is assay-source support under an older double modulator, not independent donor replication: "
    "eight doses share a plate and donor, so sign-test strata are correlated and its p is descriptive. "
    "The original donor x plate x dose stress test still fails. Tracking attrition may depend on size. "
    "No Trikafta generalization, new established discovery, or clinical prediction follows "
    "(results/external_fis_demo.json; results/finding_external_fis_demo.md).")
ER = json.load(open("results/external_fis_robustness.json"))
EM = json.load(open("results/external_fis_missingtracks.json"))
P.p(f"Well-as-unit sensitivity on the same demo plate (registered before analysis): all eight doses have "
    f"four eligible wells; {ER['R1']['n_positive']}/8 dose effects are positive, with a "
    f"median of {ER['R1']['median_effect']:+.3f}. A 2,000-draw well bootstrap *conditional on "
    f"this donor/plate* yields [{ER['R1']['well_bootstrap_ci95'][0]:+.3f}, "
    f"{ER['R1']['well_bootstrap_ci95'][1]:+.3f}], not an across-donor interval. "
    "Tracking attrition is material: among unique labeled baseline objects, large-object retention "
    f"is {ER['R2']['by_condition_size']['fsk']['large']['rate']:.1%} for vehicle and "
    f"{ER['R2']['by_condition_size']['fsk_770_809']['large']['rate']:.1%} for combo, "
    "while labeled small-object retention is 100% in both. Including unlabeled baseline objects "
    f"lowers small retention to {EM['M1']['fsk']['small']['paired_fraction_all_baseline']:.1%} and "
    f"{EM['M1']['fsk_770_809']['small']['paired_fraction_all_baseline']:.1%}. "
    "Assigning all missing outcomes to adverse observed same-cell 5th/95th percentiles gives "
    f"{EM['M2']['n_adverse_positive']}/8 positive doses, median "
    f"{EM['M2']['adverse_median']:+.3f}, but this arbitrary bound cannot rule out stronger "
    "missingness bias. Source-level evidence is conditional and must not be called donor-level "
    "validation (results/external_fis_robustness.json; results/external_fis_missingtracks.json).")
FIJI = json.load(open("results/external_fis_fiji.json"))
P.p(f"Pipeline sensitivity on the SAME images/plate/donor: the repository's alternative Fiji/ImageJ "
    f"tracking output yielded {FIJI['n_tracked_pairs']:,} paired tracks and a positive combination-minus-vehicle "
    f"size contrast at {FIJI['H1_fiji']['n_positive']}/{FIJI['H1_fiji']['n_eligible']} doses "
    f"(median {FIJI['H1_fiji']['median_effect']:+.3f}; pre-registered coverage/sign gate passes). "
    "Fiji-specific quartile cutoffs and segmentation/tracking differ from CellProfiler, so effect magnitudes "
    "are not interchangeable. In Fiji, labeled large-object retention was "
    f"{FIJI['well_and_retention_descriptive']['retention_by_condition_size']['fsk']['large']['rate']:.1%} "
    f"vehicle and {FIJI['well_and_retention_descriptive']['retention_by_condition_size']['fsk_770_809']['large']['rate']:.1%} "
    "combo, unlike CellProfiler's 91.5%/83.3%. This pipeline-dependent selection limits any effect-size "
    "claim. These two pipelines are NOT two independent biological datasets "
    "(results/external_fis_fiji.json).")
P.table(["forskolin (uM)", "paired wells / condition", "large-minus-small drug effect (log units)"],
    [[f"{d['dose_uM']:g}", f"{d['n_wells']['fsk']}/{d['n_wells']['fsk_770_809']}",
      f"{d['effect_log']:+.3f}" if d['effect_log'] is not None else "not eligible"] for d in EXT['doses']],
    "Pre-registered external single-plate demo assay; doses are correlated and not independent donors.")

DV = json.load(open("results/drevinek_fis.json"))
P.h("4.3d Distinct multi-patient well-level FIS cohort: published drug ranking reproduced", 2)
P.p(f"A genuinely distinct dataset (Drevinek et al., Zenodo 4771466) contains well-level area time series "
    f"from {DV['patients_source']} CF patient IDs on {DV['plates_source']} patient-date plates. "
    "The analysis was pre-registered and committed before downloading the source; source plate-map inspection "
    "then showed that its FskOnly wells are a reference line, not patient-matched vehicle. That provenance "
    "correction was committed before outcome analysis. Of 432 patient-plate-dose blocks, "
    f"{DV['blocks_eligible']} have >=2 wells per treatment and matched 0-60 minute measurements; "
    "one TEZ/IVA block has one well and is excluded. Equal-well area AUC is normalized to each well's "
    "time-zero area, then triple ELX/TEZ/IVA minus double TEZ/IVA is taken at patient/plate/dose level; "
    "each patient contributes the median over its blocks. "
    f"Positive in {DV['primary']['positive_patients']}/{DV['primary']['eligible_patients']} patients, "
    f"median {DV['primary']['median_patient_effect']:+.3f} baseline-area AUC (exact one-sided sign "
    f"p={DV['primary']['exact_one_sided_sign_p']:.2g}), so the pre-registered treatment-ranking H1 passes. "
    "This confirms the published in-vitro triple-over-double ordering on the paper's own public data, "
    "not a new discovery or an independent individual-organoid size test. "
    "Forty-seven of 431 plate-dose block rankings reverse if unnormalized raw area replaces normalized AUC; "
    "both patient-level summaries remain 20/20 positive. No clinical prediction follows "
    "(results/drevinek_fis.json; results/finding_drevinek_fis.md).")
P.table(["patient", "blocks", "triple-double AUC", "patient", "blocks", "triple-double AUC"],
        [[a['patient'], a['n_blocks'], f"{a['effect']:+.3f}", b['patient'], b['n_blocks'], f"{b['effect']:+.3f}"]
         for a,b in zip(DV['patient_results'][:10], DV['patient_results'][10:])],
        "Patient-equal well-level treatment-ranking test; no individual organoid areas in this accession.")

VT = json.load(open("results/vti_eti.json"))
P.page_break()
P.h("4.3e Newer VTI-versus-ETI patient ranking fails registered test", 2)
P.p(f"A further Drevinek-group public release (Zenodo 15754800) has well-level area series for "
    f"{VT['patients_source']} anonymized CF patient IDs and {VT['plates_source']} plates. "
    "Only nine patients have matched default vanzacaftor/tezacaftor/ivacaftor (VTI) and "
    "elexacaftor/tezacaftor/ivacaftor (ETI) wells; three instead have VTI-dose-optimization "
    "variants, which were not substituted into the locked primary test. Of 264 plate-dose "
    f"blocks, {VT['blocks_eligible']} qualify. Patient-equal median VTI-minus-ETI normalized "
    f"area AUC is positive in {VT['primary']['positive_patients']}/{VT['primary']['eligible_patients']} "
    f"patients (median {VT['primary']['median_patient_effect']:+.3f}, one-sided exact sign "
    f"p={VT['primary']['exact_one_sided_sign_p']:.3f}); H1 fails. Only 112/208 matched "
    "plate-dose blocks are positive. Raw-area AUC reverses 92/208 block rankings and "
    "leaves 3/9 patient medians positive. In 48 separate dose-optimization blocks, "
    "just 12 show strictly rising response across 0.02, 0.2 and 2 µM VTI, a descriptive "
    "check, not an alternative pass. The earlier 2021 dataset and this release share "
    "textual anonymized IDs 58 and 68, which neither proves nor excludes patient overlap. "
    "This negative in-vitro ranking does not establish patient equivalence/inferiority, "
    "and still cannot test individual-size effects (results/vti_eti.json; "
    "results/finding_vti_eti.md).")
P.table(["patient", "blocks", "VTI-ETI AUC"],
        [[d['patient'], d['n_blocks'], f"{d['effect']:+.3f}"] for d in VT['patient_results']],
        "Pre-registered patient-equal VTI-versus-ETI effects; negative IDs and missing default-VTI patients are retained in the full result file.")

P.h("4.4 Does size correction help theratyping? No", 2)
P.p("Locked before running: per-donor standardised Trikafta-vs-DMSO separation with raw, size-adjusted and size-band well readouts. Median "
    "raw 4.73, adjusted 4.35, band 4.39; adjusted better in 8 of 17 donors. Excluding organoids below 1,069 px (a post-hoc threshold): median "
    "4.73 -> 4.75, better in 9 of 17, Wilcoxon p = 0.68. The size trend does not improve calls, and the stricter matched-block H1 fails.")
P.page_break()
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

SA = json.load(open("results/seg_integrity_audit.json"))
P.h("4.5a Segmentation split-integrity correction", 2)
P.p(f"Training combined {SA['split_counts']['train']} train and {SA['split_counts']['val']} validation images. "
    "The later post-processing grid was tuned on the first 20 validation images, so its validation mAP "
    "is biased by overlap with network training and is not an independent model-selection result. "
    f"The disjoint {SA['split_counts']['eval']}-image eval split was held out, and its mAP@0.5 "
    f"remains {SA['heldout_eval_mAP50']:.3f} versus published {SA['published_mAP50']:.2f}. "
    f"Bootstrap across our eval images gives CI [{SA['eval_image_bootstrap_ci95'][0]:.3f}, "
    f"{SA['eval_image_bootstrap_ci95'][1]:.3f}] for our mean, but leader per-image predictions "
    "and guaranteed identical AP implementation are unavailable. No benchmark break "
    "(results/seg_integrity_audit.json).")

P.h("4.5b Clean train-only and val-only segmentation attempt: negative", 2)
TS = json.load(open("results/trainonly_seg/sealed_eval.json")); VSEL = json.load(open("results/trainonly_seg/val_selection.json"))
P.p(f"A separate protocol was frozen before fitting (results/preregistration_trainonly_seg.md): 184 training images only, "
    f"35 disjoint validation images, 12 evaluation images untouched until final selection. Sixty epochs were fit; "
    f"epoch {VSEL['selected_epoch']} minimized validation cross-entropy. A fixed grid of {VSEL['grid_size']} "
    f"postprocessing choices used the first 20 validation images only, selecting {VSEL['selected_params']} "
    f"at validation AP50 {VSEL['selected_val_mAP50']:.3f}. This validation score is a selection statistic, not an independent benchmark.")
P.p(f"The single sealed final evaluation scored {TS['mean_author']:.5f} AP50 under the authors' released Cellpose-style scorer "
    f"and {TS['mean_ours']:.5f} under our unique-mask-count scorer, below the published 0.76. The only material "
    "matcher difference arises from one released GT mask with a skipped label ID: the original scorer uses the maximum ID "
    "as object count. This small numerical difference is not a model gain. Author baseline per-image predictions were not run in "
    "this 2-CPU, no-TensorFlow/GPU workspace, so we do not assert a paired statistical comparison. Clean selection fixes an "
    "older train/val overlap, but it is not a novel model breakthrough or a benchmark match. Per-image TP/FP/FN and both AP "
    "versions are preserved in results/trainonly_seg/sealed_eval.json; selection frozen before eval at commit 120f944.")
P.table(["image #", "author AP50", "unique-count AP50", "GT masks", "predicted masks"],
        [[i, f"{d['ap_author']:.3f}", f"{d['ap_ours']:.3f}", d['n_gt'], d['n_pred']] for i, d in enumerate(TS['per_image'], 1)],
        "Sealed 12-image evaluation of the clean train-only model; image numbers follow the JSON per_image order, which has full filenames and matching counts. None selected the model or postprocessing.")

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

P.h("4.12 Orthogonal checks of the metabolic-parenchyma candidate", 2)
GS = json.load(open("results/strict_gsea_summary.json"))
P.p("GSEApy preranked GSEA used MSigDB Hallmark v2023.2 on all ranked strict-subset deficit genes, with 500 permutations. Positive normalized enrichment "
    "score (NES) means higher expression in the reference tissue than the organoids. The mean over six metabolic Hallmarks is "
    f"{GS['mean_metabolic_NES']['Liver']:.2f} in liver, {GS['mean_metabolic_NES']['Kidney - Cortex']:.2f} in kidney, "
    f"{GS['mean_metabolic_NES']['Brain - Cortex']:.2f} in brain, and {GS['mean_metabolic_NES']['Colon - Transverse']:.2f} in colon "
    f"(descriptive Mann-Whitney p={GS['mwu_metabolic_NES_liver_kidney_vs_colon_brain']:.3f}; pathways overlap and are not independent). "
    f"Brain nevertheless has oxidative phosphorylation NES {GS['metabolic_NES']['Brain - Cortex']['HALLMARK_OXIDATIVE_PHOSPHORYLATION']:.2f}, "
    "despite matching its reference tissue. A metabolic deficit alone cannot explain fidelity (results/strict_gsea_summary.json).")
OT = json.load(open("results/deficit_opentargets.json"))
P.p("Pre-registered Open Targets disease relevance. We compared each organ's top 100 deficit genes against the top 500 targets of its "
    "corresponding organ-disease term. Diagonal dominance in the five-by-five overlap matrix is:")
P.equation("D = (1/5) sum_o M_oo - (1/20) sum_{o != d} M_od")
P.p(f"D={OT['D']:.1f}, permutation p={OT['perm_p']:.4g}: the primary test passes (results/deficit_opentargets.json). "
    f"But liver-specific disease enrichment fails ({OT['per_organ']['Liver']['own_hits']}/100 own-disease hits, "
    f"p={OT['per_organ']['Liver']['p']:.2f}). The pass is driven more by lung and colon than liver; the few hits per cell limit its use as a discovery claim.")
DC = json.load(open("results/deficit_decoupler.json"))
P.p(f"Pre-registered decoupler univariate linear modeling used OmniPath CollecTRI transcription factor regulons and PROGENy footprints "
    f"on {DC['n_genes']:,} shared genes. Only {int(DC['pooled_hits'])}/{DC['pooled_tested']} tested lineage master factors ranked in "
    f"the top decile of deficient factors (binomial p={DC['binom_p']:.2f}); the positive control FAILS. "
    "Exploratory PROGENy scores show PI3K higher and JAK-STAT and p53 lower in organoids across all five organs, "
    "a shared culture signature rather than an organ-specific mechanism (results/deficit_decoupler.json).")
BM = json.load(open("results/deficit_bodymap.json"))
rr = [v['spearman'] for v in BM['per_organ'].values()]; overlap = [v['top100_overlap'] for v in BM['per_organ'].values()]
P.p(f"Reference robustness. We replaced GTEx with Expression Atlas Illumina Body Map E-MTAB-513. The pre-registered test passed: "
    f"deficit-vector Spearman rho {min(rr):.2f}-{max(rr):.2f} across organs, top-100 overlap {min(overlap)}-{max(overlap)}. "
    f"The Open Targets diagonal-dominance rerun also passes (D={BM['opentargets_rerun']['D']:.1f}, "
    f"permutation p={BM['opentargets_rerun']['perm_p']:.4f}). Both deficit vectors share the organoid measurements, "
    "so this checks reference sensitivity, not independent biological replication (results/deficit_bodymap.json).")
CP = json.load(open("results/deficit_clinpgx.json"))
P.p("ClinPGx pharmacogenes: the pre-registered VIP-or-CPIC gene set is unusable because the source release marks all 25,041 genes "
    "as VIP. We retain that failure. A documented, post-registration CPIC-guideline-only deviation finds liver deficits in "
    f"{CP['liver']['n_pg_tested']} tested genes (median {CP['liver']['median_deficit_pg']:.2f} vs "
    f"{CP['liver']['median_deficit_other']:.2f} in other genes, p={CP['liver']['mwu_p_greater']:.4f}), "
    f"but liver specificity fails (p={CP['H2']['perm_p']:.2f}, only {CP['H2']['n_pg_all_organs']} genes shared across all organ vectors). "
    "This is a drug-testing caution, not evidence of a liver-specific mechanism (results/deficit_clinpgx.json).")
GN = json.load(open("results/deficit_gnomad.json")); ENB = json.load(open("results/deficit_ensembl_biotype.json"))
P.p(f"Pre-registered population-constraint and coding-composition checks were negative. In gnomAD v2.1.1, liver deficit genes have median "
    f"LOEUF {GN['per_organ']['Liver']['median_loeuf']:.3f} vs brain {GN['per_organ']['Brain - Cortex']['median_loeuf']:.3f} "
    f"(one-sided p={GN['H1']['p_one_sided']:.3f}; lower LOEUF is more loss-of-function constrained). "
    f"Ensembl REST mapped all fixed top-100 genes per organ: liver {ENB['per_organ']['Liver']['n_protein_coding']}/100 and brain "
    f"{ENB['per_organ']['Brain - Cortex']['n_protein_coding']}/100 are protein coding (one-sided Fisher p={ENB['H1']['p_one_sided']:.1f}). "
    "Neither check supplies a mechanism; the gene universe is already coding-enriched (results/deficit_gnomad.json; "
    "results/deficit_ensembl_biotype.json).")
P.p("Synthesis. Organoid deficits are reproducible to a reference swap and linked to organ-disease targets in aggregate, but the "
    "metabolic-parenchyma gap lacks liver specificity and the kidney finding fails external replication. It remains a falsifiable "
    "candidate, not a new discovery or a validated diagnostic tool.")
P.h("4.13 Independent protein and pathway annotation, single-cell challenge, and prior art", 2)
WP = json.load(open("results/deficit_wikipathways.json"))
P.p(f"Independent pathway annotation (WikiPathways 2026-09-10, NCBI Gene mapping). A pre-registered fixed-name filter "
    f"selected {WP['n_selected_pathways']} pathways containing 'bile acid', 'fatty acid', 'drug metabolism' or 'xenobiotic metabolism'. "
    f"Their union intersects {WP['per_organ']['Liver']['n_union']}/100 liver versus "
    f"{WP['per_organ']['Brain - Cortex']['n_union']}/100 brain top deficits "
    f"(Fisher OR={WP['H1']['odds_ratio']:.2f}, p={WP['H1']['p_one_sided']:.4f}; test passes). "
    "This supports the liver top-deficit pathway label in a second resource, not causation or organoid specificity; "
    "the comparison follows earlier exploratory results (results/deficit_wikipathways.json).")
GOA = json.load(open("results/deficit_goa.json"))
KE = json.load(open("results/deficit_kegg.json"))
RH = json.load(open("results/deficit_rhea.json"))
P.table(["curated source", "selected definition", "liver top-100", "brain top-100", "one-sided Fisher p"], [
    ["EBI GOA", "3 direct metabolic BP terms", GOA["per_organ"]["Liver"]["n_direct_metabolic_union"], GOA["per_organ"]["Brain - Cortex"]["n_direct_metabolic_union"], f"{GOA['H1']['p_one_sided']:.4f}"],
    ["KEGG REST", "9 fixed-title metabolic pathways", KE["per_organ"]["Liver"]["n_union"], KE["per_organ"]["Brain - Cortex"]["n_union"], f"{KE['H1']['p_one_sided']:.4f}"],
    ["Rhea", "curated Swiss-Prot reaction mapping", RH["per_organ"]["Liver"]["n_rhea_positive"], RH["per_organ"]["Brain - Cortex"]["n_rhea_positive"], f"{RH['H1']['p_one_sided']:.2g}"],
], "Further pre-registered annotation checks on fixed 100-gene deficit lists, not independent experiments. Each p is descriptive in an exploratory sequence; results/deficit_goa.json, deficit_kegg.json, deficit_rhea.json.")
P.p(f"Direct GOA BP terms cover {GOA['per_organ']['Liver']['n_direct_metabolic_union']}/100 liver versus "
    f"{GOA['per_organ']['Brain - Cortex']['n_direct_metabolic_union']}/100 brain top deficits, while any BP annotation "
    f"covers 94 and 90 respectively. KEGG REST selected {KE['n_selected_pathways']} pathways with {KE['n_union_symbols']} "
    "mapped genes. Rhea reactions mark 42 liver versus 10 brain genes; kidney has 41, so the reaction signal "
    "is not specific to liver. GOA, KEGG, Rhea, UniProt and WikiPathways annotations may share literature and curation, "
    "and prior g:Profiler work included KEGG. These positives cannot be counted as independent confirmations. "
    "Neither organoid protein nor enzyme activity was measured, and the original kidney fidelity candidate failed external replication. "
    "No mechanism, discovery, or diagnostic use follows from these annotation contrasts.")
IA = json.load(open("results/deficit_intact.json"))
JA = json.load(open("results/deficit_jaspar.json"))
P.p(f"Additional pre-registered annotation/coverage checks: IntAct PSIQUIC queried the fixed top-20 deficit "
    f"genes per group and measured {IA['per_organ']['Liver']['n_measured']} liver and "
    f"{IA['per_organ']['Brain - Cortex']['n_measured']} brain genes. Median indexed interaction RECORD "
    f"counts were {IA['per_organ']['Liver']['median_count']:.1f} versus "
    f"{IA['per_organ']['Brain - Cortex']['median_count']:.1f} (one-sided rank-test "
    f"p={IA['H1']['p_one_sided']:.4f}). Study density, nonhuman interactors, and blood/immune genes confound "
    "this contrast; no organoid interactome was measured. JASPAR CORE vertebrate matrices cover 2/3 "
    "pre-fixed brain lineage TFs absent from CollecTRI (NEUROD2 and TBR1, not NEUROD6) and kidney PAX2; "
    "this is motif availability only and cannot repair the failed decoupler activity control. "
    "Neither is an independent biological replication (results/deficit_intact.json; results/deficit_jaspar.json).")
HMS = json.load(open("results/deficit_hpa_ms.json"))
P.p(f"Adult-tissue protein check (HPA v25.1 mass spectrometry). Liver protein is detected in "
    f"{HMS['per_organ']['Liver']['n_liver_detected']}/100 liver versus "
    f"{HMS['per_organ']['Brain - Cortex']['n_liver_detected']}/100 brain top deficits "
    f"(Fisher OR={HMS['H1']['odds_ratio']:.2f}, p={HMS['H1']['p_one_sided']:.1e}; pre-registered H1 passes). "
    "Only 26 liver and 22 brain list proteins have measured intensity in both adult liver and cortex, "
    "below the registered 50-per-list minimum, so the quantitative H2 is uninterpretable. "
    "There is no organoid proteome in this check (results/deficit_hpa_ms.json).")
SC = json.load(open("results/sc_kidney_markers.json"))
P.p(f"Single-cell challenge on an independent public human kidney-organoid accession, GSE108291. Scanpy sparse QC retains "
    f"{SC['runs']['org']['n_qc']:,} cells in the main run and {SC['runs']['org4']['n_qc']:,} in a second run; "
    f"no cell coexpresses ALDOB and SLC17A3 (zero in both). The registered 1% two-marker H1 fails; the "
    "<20% H2 passes. This does not show that proximal-tubule cells are absent, because sparse single-cell dropout and "
    "developmental state remain untested. QC thresholds were set after registration but before looking at marker coexpression; "
    "both runs count as one accession (results/sc_kidney_markers.json).")
OLS = json.load(open("results/ols_anatomy_audit.json"))
P.p(f"Ontology audit (EBI OLS4) distinguishes kidney cortex UBERON:0001225 from kidney UBERON:0002113, "
    "and cholangiocyte CL:1000488 from hepatocyte CL:0000182. This is a label check only, not proof of "
    "anatomical hierarchy or expression similarity; the exact 'proximal tubule cell' query was unresolved "
    "(results/ols_anatomy_audit.json).")
HPO = json.load(open("results/deficit_hpo.json"))
P.p(f"HPO direct gene-phenotype annotation could not test liver-vs-brain metabolic phenotype enrichment: "
    f"only {HPO['per_organ']['Liver']['n_any_hpo']}/100 liver and "
    f"{HPO['per_organ']['Brain - Cortex']['n_any_hpo']}/100 brain deficits have any HPO annotation, "
    "below the pre-registered 50/list coverage floor. No clinical inference from the missing annotations "
    "(results/deficit_hpo.json).")
P.p("InterPro domain audit is invalid as a top-deficit comparison. Its registered six liver pharmacogenes "
    "were incorrectly called top-100 deficits: only CYP2D6 actually is one. The protein-domain annotations "
    "remain correct, but the primary contrast cannot be scored (results/deficit_interpro.json).")
EP = json.load(open("results/epmc_priorart.json"))
P.p(f"A fixed three-query Europe PMC title/abstract search returned {EP['n_unique_returned']} distinct records. "
    "None explicitly reported human intestinal CF organoid starting-size moderation of modulator-induced FIS "
    "in its abstract, but absence in a narrow search does not establish novelty. Related work links nasal "
    "organoid lumen to baseline CFTR function and pig pancreatic organoid size to absent lumen "
    "(results/epmc_priorart_review.md).")
P.p("Wider full-text prior art narrows the size claim further: Calucho et al. 2021 tested individual WT nasospheroid starting size against CFTR-dependent FSK shrinking and found no significant association among 138 responders (Spearman r=-0.1216, p=0.1537 for 60-minute fractional response; supplementary slope r=0.04124, p=0.6298). This airway system has opposite polarity/readout and does not directly duplicate donor-matched intestinal drug-by-size interaction, but starting-size analysis of CFTR spheroids is not generally novel. Kim et al. 2020 found no initial-size/growth relation in two colorectal tumor organoid lines and warned that small-object detection can add variation. Full-text audit and source URLs: results/epmc_priorart_review.md and results/external_fis_source_audit.md. No discovery follows from a narrower question alone.")
P.p("Berical et al. 2022 tracked individual human iPSC-derived airway spheroids after forskolin or vehicle and tested CFTR-modulator FIS across patient lines. This is strong prior art for individual-object CFTR swelling assays, not a release of matched individual-object baseline/end size by patient and treatment: its public Source Data workbook provides figure-level experiment values and separate baseline sizes without object IDs linking the panels. Airway tissue and a 20-24-hour response differ from intestinal FIS. Our narrower intestinal drug-by-size result remains unvalidated, not novel merely because this released workbook lacks paired tracks (results/berical_source_suitability.md).")
P.p("Synthesis. The liver top-deficit list has independent pathway and adult-tissue protein support, "
    "while no new organoid proteome, matched tissue, perturbation, or clinical outcomes validate its effect on fidelity. "
    "The kidney external profiles recover their tissue despite the strict single-cell marker test failing. "
    "The candidate remains unproven; do not offer the tool as diagnosis.")

P.p("External replication source audit (results/external_fis_source_audit.md): a 2025 respiratory FIS study "
    "with single-organoid measures excluded baseline areas below 1500 px and used forskolin alone rather "
    "than modulator vs DMSO; OrganoID's open individual time series used pancreatic cancer organoids "
    "and gemcitabine; two nasal CFTR-modulator papers reported either well-summed areas or summary plots "
    "without a verified public matched-object table. Drevinek's 20-patient area series was suitable for "
    "well-level drug-ranking analysis, not individual-size replication. This search is not exhaustive.")

P.h("5. Negative results (kept by design)")
for t in ["The initial power-law (alpha > 1) interpretation was withdrawn after large-organoid checks contradicted it.",
          "Size-adjusted and size-filtered readouts do not improve donor-level modulator discrimination.",
          "Matched donor x plate x dose robustness H1 fails: 9/12 positive donor medians, p=0.073, despite the pooled later-date result.",
          "Small-organoid replication on the public FIS time series is impossible: it is well-level; the 20-patient triple-vs-double ranking passes, but the separate VTI-vs-ETI ranking fails (6/9, p=0.254).",
          "Our segmentation is below the published state of the art.",
          "Fidelity: 'Pancreas' as best match (35/151 series with S3000) disappears with other methods, so it is a method artefact.",
          "Fidelity: a 'liver disease' series (GSE278954) ranks Liver 50th of 54, and lung and breast organoids are rarely matched to their organ.",
          "Fidelity: organ labels come from text, and series profiles average all samples, including any non-organoid controls.",
          "Fidelity: the HPA parenchymal share of deficit genes does not clearly separate faithful from unfaithful organs (p = 0.061); lung organoids break the pattern.",
          "Fidelity: all organoid types lack secreted proteins (UniProt), so the secretome does not explain the organ split.",
          "Replication: the kidney part of the organ split failed out of sample (two ArrayExpress kidney datasets rank kidney 1st and 2nd).",
          "GSEA: metabolic deficits are also present in faithful brain organoids; not specific to poor fidelity.",
          "Open Targets: liver deficits show no liver-disease enrichment despite the aggregate diagonal-dominance pass.",
          "decoupler: lineage master-TF positive control fails (2/12); pathway signals are shared across cultures.",
          "ClinPGx: pre-registered VIP flag is defective; CPIC-only deviation passes liver deficit but fails liver specificity.",
          "gnomAD and Ensembl: no liver-vs-brain differentiation by LOEUF or protein-coding composition.",
          "Single-cell GSE108291: no ALDOB/SLC17A3 coexpression under the registered two-marker definition; dropout is an alternative.",
          "HPO: too few annotated genes for the registered clinical phenotype contrast.",
          "InterPro: primary contrast invalid due to mistaken liver top-100 gene-set membership.",
          "Fidelity: the culture-fibroblast attractor is RETRACTED as a general organoid property: on 45 series whose samples are all organoids it shrinks to 9-29% of series (Section 4.9)."]:
    P.p("- " + t)
P.h("6. Discussion")
P.p("In this accession, drug-minus-DMSO swelling tends to rise from small to larger organoids and then plateau. "
    "The plate-and-dose-matched donor sign test fails at alpha 0.05, so we do not establish a reliable donor-general "
    "effect. Uniform surface-flux geometry predicts the reverse direction, but that simple model is not a full "
    "measurement-bias test. Lumen maturity and cell composition are hypotheses, not measured explanations. "
    "Falsifiable prediction: in an independent, matched, time-lapse single-organoid intestinal FIS cohort, with "
    "size thresholds locked before analysis, a donor-level modulator-minus-DMSO large-minus-small effect will "
    "be positive in most donors, with a sign-test p<0.05. This needs a second suitable dataset and mechanistic "
    "measurements; the public leads audited here do not supply those controls.")
P.p("Practical note. In this dataset, starting size may affect single-organoid readouts, but size correction did not "
    "improve per-donor discrimination. A lab comparing assays across passages may record its starting-size "
    "distribution as a quality-control variable; no clinical adjustment is recommended from these data.")
P.p("Novelty remains unproved. We did not find the exact intestinal drug-by-starting-size interaction in the bounded literature search, but Calucho et al. 2021 tested starting size versus individual CFTR nasospheroid FSK response and found no association. Berical et al. 2022 tracked individual iPSC-derived airway spheroid FIS and tested modulators across multiple patient lines; the released source workbook does not link patient-treatment outcomes to per-object baseline areas. A nasal 2D-derived-organoid study reported swelling variation in large structures (bioRxiv 2021.07.20.453105). The narrower intestinal question is unconfirmed and this was not a systematic review.")
P.h("7. Tools used")
TL = pd.read_csv("results/tools_ledger.csv")
TL["gate"] = TL.counts_for_gate.astype(str).map({"True": "counts", "False": "infra (excluded)"})
P.table(["tool", "kind", "where used", "gate"], TL[["tool", "kind", "where_used", "gate"]].values.tolist(), f"Tools genuinely used (results/tools_ledger.csv): {len(TL)} entries, {int((TL.gate == 'counts').sum())} counting toward the gate after excluding infrastructure. Target of 40 reached; this counts resources genuinely used, not independent biological confirmations.")
P.h("References")
for r in ["Lefferts JW, et al. OrgaSegment: deep-learning based organoid segmentation to quantify CFTR dependent fluid secretion. Commun Biol 2024. https://www.nature.com/articles/s42003-024-05966-4",
          "Boj SF, Vonk AM, et al. Forskolin-induced swelling in intestinal organoids: an in vitro assay for assessing drug response in cystic fibrosis patients. J Vis Exp 2017. https://pmc.ncbi.nlm.nih.gov/articles/PMC5408767/",
          "GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 2020. https://gtexportal.org",
          "Karlsson M, et al. A single-cell type transcriptomics map of human tissues. Sci Adv 2021. https://www.proteinatlas.org",
          "Szklarczyk D, et al. The STRING database in 2023. Nucleic Acids Res 2023.",
          "The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Res 2023.",
          "Moreno P, et al. Expression Atlas and ArrayExpress/BioStudies at EMBL-EBI. Nucleic Acids Res 2022.",
          "Open Targets Platform. https://platform.opentargets.org/",
          "Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinform Adv 2022. https://bioconductor.org/packages/decoupleR/",
          "Turei D, et al. OmniPath: guidelines and gateway for literature-curated signaling pathway resources. Nat Methods 2016. https://omnipathdb.org/",
          "ClinPGx / PharmGKB. https://www.clinpgx.org/",
          "Karczewski KJ, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 2020. https://gnomad.broadinstitute.org/help/constraint",
          "Ensembl REST API. https://rest.ensembl.org/documentation/info/lookup_post",
          "WikiPathways monthly releases. https://data.wikipathways.org/20260910/gmt/",
          "EBI GOA human GAF. https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz",
          "IntAct PSIQUIC. https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/",
          "JASPAR CORE vertebrate motif API. https://jaspar.elixir.no/api/v1/matrix/",
          "KEGG REST API. https://www.kegg.jp/kegg/rest/keggapi.html",
          "Rhea reaction-to-Swiss-Prot mapping. https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_sprot.tsv",
          "Human Protein Atlas protein mass-spectrometry normal tissue download. https://www.proteinatlas.org/humanproteome/tissue/data",
          "Human Phenotype Ontology gene-to-phenotype annotation. https://obophenotype.github.io/human-phenotype-ontology/annotations/genes_to_phenotype/",
          "EBI Ontology Lookup Service 4. https://www.ebi.ac.uk/ols4/api-docs",
          "Human Cell Atlas GSE108291 kidney organoid single-cell project. https://explore.data.humancellatlas.org/projects/7b947aa2-43a7-4082-afff-222a3e3a4635",
          "Demchenko A, et al. A semi-automated algorithm for image analysis of respiratory organoids. PLoS Comput Biol 2025. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013589",
          "Matthews J, et al. OrganoID. PLoS Comput Biol 2022. https://pmc.ncbi.nlm.nih.gov/articles/PMC9645660/",
          "Anderson JD, et al. CFTR function and clinical response to modulators parallel nasal epithelial organoid swelling. AJP Lung 2021. https://pmc.ncbi.nlm.nih.gov/articles/PMC8321858/",
          "Calucho M, et al. Validation of nasospheroids to assay CFTR functionality and modulator responses in cystic fibrosis. Sci Rep 2021. https://www.nature.com/articles/s41598-021-94798-x",
          "Berical A, et al. A multimodal iPSC platform for cystic fibrosis drug testing. Nat Commun 2022. https://www.nature.com/articles/s41467-022-31854-8 . Source-data suitability: results/berical_source_suitability.md",
          "Kim S, et al. Comparison of Cell and Organoid-Level Analysis of Patient-Derived 3D Organoids to Evaluate Tumor Cell Growth Dynamics and Drug Response. SLAS Discov 2020. https://www.slas-discovery.org/article/S2472-5552(22)06605-9/fulltext",
          "Botelho H, Hagemeijer MC, et al. FIS_image_analysis demonstration dataset. https://github.com/hmbotelho/FIS_image_analysis",
          "Drevinek P, et al. Response to elexacaftor/tezacaftor/ivacaftor in intestinal organoids derived from people with cystic fibrosis. J Cyst Fibros 2021. https://doi.org/10.1016/j.jcf.2021.07.006 . Data: https://zenodo.org/records/4771466",
          "Borek-Dohalska L, et al. Effect of vanzacaftor on cystic fibrosis airway epithelial cells compared to elexacaftor. Data: https://zenodo.org/records/15754800",
          "InterPro API. https://interpro-documentation.readthedocs.io/en/latest/api.html",
          "Europe PMC RESTful Web Service. https://europepmc.org/RestfulWebService",
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
# A. Auditable patient-block outcomes, not additional biological cohorts.
# Generated entirely from the frozen full-result JSONs. Keep excluded blocks visible.
P.page_break()
P.h('Appendix A. Patient, plate and dose-level treatment contrasts')
P.p('This appendix makes the two distinct well-level FIS rankings inspectable without downloading a raw table. Every row is one source patient x plate x forskolin-dose block, not a new patient or accession. Effects average technical wells within an arm before the patient median is formed. Normalized area AUC is dimensionless and divided by the observation horizon; values here are differences of treatment-arm AUCs. These results cannot test individual-organoid size and should not be used to choose therapy.')
P.p('A1 uses Zenodo 4771466 and triple ELX/TEZ/IVA minus double TEZ/IVA. A2 uses Zenodo 15754800 and default VTI minus ETI. The source releases come from related teams, and anonymized patient numbers may overlap. The two experiments, drugs and eligibility rules are not pooled. See results/preregistration_drevinek_fis.md and results/preregistration_vti_eti.md for the frozen questions; the source scripts and JSON contain every eligible well, time-point count, exclusion and sensitivity test.')
for label, obj, arm, key in [
    ('A1', DV, 'ELX/TEZ/IVA minus TEZ/IVA', 'drevinek_fis'),
    ('A2', VT, 'VTI minus ETI', 'vti_eti'),
]:
    P.h(f'{label}. {arm}', 2)
    P.p(f"Source {obj['source']}. {obj['patients_source']} IDs and {obj['plates_source']} plates appear in the public source; {obj['blocks_eligible']}/{obj['blocks_total']} blocks meet its preregistered primary criteria. The table preserves positive and negative contrasts. The registered patient-level result is {obj['primary']['positive_patients']}/{obj['primary']['eligible_patients']} positive medians, exact one-sided p={obj['primary']['exact_one_sided_sign_p']:.7g}. Unnormalized raw-area AUC is shown only to reveal sensitivity to baseline area, never as an independent trial.")
    pats = sorted(set(b['patient'] for b in obj['blocks']))
    for ip, patient in enumerate(pats):
        rows = sorted((b for b in obj['blocks'] if b['patient']==patient), key=lambda b:(str(b['date']),b['dose']))
        primary = next((p for p in obj['patient_results'] if p['patient']==patient), None)
        P.h(f'{label}.{ip+1}. Source patient {patient}', 3)
        if primary is None:
            P.p(f'No eligible primary block for this ID under the fixed treatment arms and time/well criteria. {len(rows)} source blocks are retained below as exclusions; this ID contributes no patient median or sign-test observation.')
        else:
            P.p(f"{primary['n_blocks']} eligible of {len(rows)} blocks; median normalized treatment contrast {primary['effect']:+.4f}; median raw-area contrast {primary['raw_effect']:+.2f}. This is one patient-weighted observation regardless of the number of doses or plates. A negative or small block remains in the table rather than being filtered by outcome.")
        data = []
        for b in rows:
            why = '; '.join(b['reasons'])
            data.append([str(b['date']), f"{b['dose']:.5g}", f"{b['effect']:+.4f}" if not why else 'excluded', f"{b['raw_effect']:+.1f}" if not why else why.replace('fewer_than_three_shared_treatment_times_including_zero','missing matched times').replace('fewer_than_two_eligible_wells_per_arm','<2 wells/arm')])
        P.table(['plate date','fsk (uM)','normalized delta','raw delta or exclusion'], data,
                f'{label} source patient {patient}: {len(rows)} plate-dose blocks from results/{key}.json. Raw-area delta is not baseline normalized.')
P.p('Audit boundary. A plate-dose block is not an independent donor. Within each patient the median across eligible blocks gives one observation for the exact sign test, preventing a large plate or many concentrations from dominating the patient count. Raw-area reversals flag sensitivity to baseline organoid density. For excluded A2 patient blocks, a missing default treatment arm is not replaced with a different VTI concentration. This appendix is a reproducibility trace of two already published in-vitro comparisons, not evidence for clinical benefit, equivalence, or the single-object size candidate.')

P.save("paper/mega27-03-virtual-organoid-paper.docx")
print("ok", P.eq, "eq", P.tab, "tab", P.fig, "fig")
