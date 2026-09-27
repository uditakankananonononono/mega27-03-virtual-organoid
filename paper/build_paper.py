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

P = Paper("Is Starting Organoid Size Associated with the CFTR-Modulator Swelling Contrast? "
          "A Size-Aware Computational Reanalysis with a Track-Attrition Audit",
          "MEGA-PROGRAM-27, Item 3 - Udita Phookan (program owner). Working draft of 27 September 2026.")
P.h("Abstract")
P.p("We ask one primary question: is starting organoid size associated with the CFTR-modulator swelling contrast among tracked intestinal organoids, and what does track attrition leave unidentified? The primary endpoint is the donor-level large-minus-small elexacaftor/tezacaftor/ivacaftor minus DMSO log-swelling contrast, and the primary analysis is the pre-registered donor-block mixed model with an exact donor-level sign-flip test, reported together with its companion pre-registered tests. All other analyses in this paper are exploratory and are labeled as such.")
P.p("The primary pre-registered tests, stated first: the sign-only donor test FAILED (9/12 positive, p=0.073), the zero-forskolin specificity test FAILED (3/11, p=0.967), and the matched plate-dose stress test FAILED. Under the separately locked donor-block analysis (31 eligible donor x experiment x dose blocks from 12 donors), the magnitude-aware treatment-by-large coefficient is +0.252815 (positive-direction one-sided p=0.00001151; 4,096-assignment exact donor-level sign flip p=0.002930), and the exploratory donor-matched octile contrast rises from 0.336 [0.262, 0.405] in the smallest starting-size octile toward 0.55-0.64 in larger octiles (top minus bottom 0.251 [0.184, 0.323], donor bootstrap). These are within-accession, magnitude-aware results, not independent replication or a clinical predictor.")
P.p("The methodological contribution is an explicit audit of track attrition. In the original eligible blocks and size cells, 6,964 baseline-detected organoids include 2,000 without a positive measured endpoint; observability varies by size, treatment and image context. Donor-clustered observability models and missing-endpoint scenarios show that the direction of the contrast in the full baseline-detected population is not identified: the scenario envelope crosses zero. The cross-fit inclusion weights address a narrower selection step among already-positive endpoints, not those 2,000 missing outcomes. This is a reusable audit design, not a correction that makes the biological size effect robust.")
P.p("A separate segmentation arm is reported as a secondary exploratory analysis: a train-only U-Net reached author-scored AP50 0.761108 against a rounded published 0.76, but repeated acquisition fields cross the train/evaluation and train/validation filenames, so no independent benchmark-beat claim is made and an acquisition-candidate-group-disjoint exploratory reassignment was subsequently run on the same 231-file source set, including 12 previously scored eval files. This paper is an in vitro computational reanalysis only, not a diagnostic or treatment predictor. The repository preserves all gate outcomes, code and negative controls.")
P.h("1. Introduction")
P.p("CFTR moves chloride and bicarbonate across the apical membrane of epithelial cells; water follows, and in a closed organoid the lumen "
    "swells. Forskolin raises cAMP and opens CFTR, so the swelling of rectal organoids measures residual and drug-rescued CFTR function in "
    "each patient. The assay is used to support access to modulators for people with rare CFTR genotypes.")
P.p("Organoids in a well vary in size by more than an order of magnitude. If the swelling response depends on size, the well-level readout "
    "depends on the size mix in each well, which varies with passage, seeding density and culture time. We asked whether the single-organoid "
    "observed response differs by starting size, and how attrition before endpoint measurement limits that inference. We do not show that size correction improves drug-response calls.")

P.h("1.1 Primary question, primary endpoint and primary analysis", 2)
P.p("Primary question: is starting organoid size associated with the CFTR-modulator swelling contrast among tracked intestinal organoids, and what does track attrition leave unidentified? Primary endpoint: the donor-level large-minus-small Trikafta-minus-DMSO log-swelling contrast. Primary analysis: the pre-registered donor-block mixed model plus the 4,096-assignment exact donor-level sign-flip test, interpreted alongside the companion pre-registered sign-only, zero-forskolin specificity and matched plate-dose stress tests. This framing is a post-verdict amendment of 27 September 2026 (owner-provided verdict archived with message-ID provenance in results/judge_rounds/; amendment queue locked before execution in notes/judge_verdict_01_amendment_queue.md). It changes presentation order and emphasis only: every pre-registered gate, outcome and negative control is unchanged, and the failed primary tests are reported in Section 4 before any magnitude-aware result.")
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

P.p("An additional pre-specified eligibility sensitivity used the same accession and original size bins, varying only minimum organoid count in every treatment-by-size cell. At five per cell, the result weakens to 8/12 positive donors (one-sided p=0.193848). At ten, six of six donors are positive (nominal p=0.015625), but six is below the frozen eight-donor eligibility floor, so that setting is UNINTERPRETABLE, not a rescue. At twenty, only three donors remain. A frozen post hoc donor-composition audit shows all three min-three nonpositive donors (PDIO_01, PDIO_16, PDIO_18) disappear at minimum ten; the six survivors were already positive at min three, and three other positive donors also disappear. This is selective eligibility, not evidence that excluded donors are nonresponders. The frozen three-per-cell H1 still FAILS (results/block_eligibility_sensitivity.json; results/finding_eligibility_composition.md).")
P.table(["minimum objects/cell", "eligible blocks", "donors", "positive", "one-sided p", "verdict"], [[3,31,12,9,".072998","FAIL"],[5,28,12,8,".193848","FAIL"],[10,15,6,6,".015625","n<8"],[20,9,3,3,".125000","n<8"]], "Same-accession eligibility sensitivity preregistered before rerun at b9978eb; the nominal six-donor p at minimum ten does not meet the frozen eligible-donor floor or constitute independent validation.")
Z0 = json.load(open("results/zero_fsk_control.json"))
P.h("4.3c No-forskolin specificity stress test: fails", 2)
P.p("An additional frozen within-accession check asked whether the drug-by-starting-size contrast is more positive at 0.128 uM forskolin than in the same donor and experiment at zero forskolin. Both doses required >=3 objects in every treatment-by-size cell. Of 50 donor-experiment combinations, 24 pairs across 11 donors qualify; six source donors lack an eligible pair. Only 3/11 donor-median paired deltas are positive, exact one-sided sign p=0.967285, so the preregistered specificity H1 FAILS. A no-forskolin size interaction is often larger, weakening a CFTR-stimulation-specific interpretation. This does not prove drug inactivity or identify the cause of the zero-dose pattern: pretreatment, baseline biology and size-dependent tracking remain possible. The same accession and post hoc size cutoffs preclude independent validation (results/zero_fsk_control.json; results/finding_zero_fsk_control.md).")
P.table(["donor", "paired plates", "effect 0 uM", "effect .128 uM", "paired delta"],
        [[d["donor"],d["paired_experiments"],f'{d["median_effect_0"]:+.3f}',f'{d["median_effect_0128"]:+.3f}',f'{d["median_delta"]:+.3f}'] for d in Z0["donors"]],
        "Same-accession specificity check, donor-equal medians: positive delta means greater size interaction at .128 versus zero forskolin. Plate-paired deltas are aggregated separately, so they need not equal the difference of the two displayed median effects. No clinical interpretation.")


EXT = json.load(open("results/external_fis_demo.json"))
P.h("4.3d External single-plate intestinal FIS assay", 2)
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
P.h("Secondary exploratory analyses I: well-level FIS drug-ranking replications", 2)
P.p("Sections 4.3d (well-level ranking) and 4.3e below are off-question secondary analyses, retained for completeness and labeled exploratory.")
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

P.h("4.3f VO3-R2 donor-block magnitude-aware reanalysis", 2)
LMM = json.load(open("results/size_lmm.json")); assert LMM['verdict']=="DISCOVERY" and LMM['n_donors']==12
P.p(f"The original preregistered sign test, 9/12 positive donors and p=0.073, remains a FAIL. A new dated amendment A1 in notes/prereg_size_reanalysis.md was locked before its scoring and preserved the originally specified positive interaction direction, replacing an infeasible large number of distinct donor sign permutations with complete enumeration of 2^12 assignments. The analysis uses {LMM['n_blocks']} eligible donor x experiment x dose blocks from {LMM['n_donors']} donors and {LMM['n_cell_rows']} block-cell means, with the original 722/1400 starting-area cutoffs and at least three objects in every cell. The mixed model is log swelling ~ treatment*large with donor and experiment random intercepts; the treatment-by-large beta is {LMM['lmm']['interaction_beta']:+.6f}, full-versus-reduced likelihood-ratio statistic {LMM['lmm']['lrt_stat']:.4f}, positive-direction one-sided p={LMM['lmm']['lrt_p_one_sided_positive']:.8f}. Both fits converged.")
P.p(f"The separate magnitude-aware exact donor-level sign-flip test averages per-donor median block effects and enumerates {LMM['exact_donor_signflip']['n_assignments']} sign assignments. The observed mean median effect is {LMM['exact_donor_signflip']['observed_mean_median_effect']:+.6f}, positive-direction exact p={LMM['exact_donor_signflip']['p_one_sided_positive']:.8f}. Both new tests meet their pre-scoring thresholds (results/size_lmm.json). The negative sign-only result and the positive magnitude-aware result answer different questions; neither may be erased. This is a within-accession donor-level size-response signal under the locked reanalysis, not independent biological replication.")
P.p("Limitations: size cutoffs were selected on this same accession before VO3-R2; LMM cell means are treated as homoscedastic despite unequal cell counts; selection and track attrition can create a size trend; the prior no-forskolin specificity failure still limits a CFTR-specific story. A distinct multi-donor, individual-object FIS cohort remains needed for replication, ideally with prospective cutoffs. No clinical treatment recommendation follows from the result.")

P.h("4.3g Exploratory measurement-invariance diagnostics after judge review", 2)
MI = json.load(open("results/measurement_invariance.json")); assert MI['n_blocks']==31 and MI['n_donors']==12
P.p(f"A later external critique asked whether the donor-size interaction could reflect unequal cell counts, influential donors, segmentation area error, or selective track retention. We froze notes/prereg_measurement_invariance.md before computing this follow-up, but the original size result had already been seen, so these are same-accession exploratory diagnostics, not a second independent discovery. In the {MI['n_blocks']} original eligible blocks, cell counts span {MI['cell_count_range'][0]} to {MI['cell_count_range'][1]} (median {MI['cell_count_median']:.0f}); this is a real departure from equal-variance cell-mean intuition. We recomputed an equal-donor mean of within-donor median block effects = {MI['equal_donor_mean_median']:+.6f}. A seed-{MI['donor_bootstrap_seed']} donor-resampling percentile interval over {MI['donor_bootstrap_draws']:,} draws is [{MI['donor_bootstrap_ci95'][0]:+.6f}, {MI['donor_bootstrap_ci95'][1]:+.6f}]. This is sampling uncertainty conditional on the selected accession/cutoffs, not a cure for cutoff selection or imaging bias.")
_lo=min(MI['leave_one_donor_out_equal_means'].items(),key=lambda q:q[1]);_hi=max(MI['leave_one_donor_out_equal_means'].items(),key=lambda q:q[1])
P.p(f"Leave-one-donor-out equal means range from {_lo[1]:+.6f} (omitting {_lo[0]}) to {_hi[1]:+.6f} (omitting {_hi[0]}). Within-cell 20% trimmed means, aggregated as donor-median block effects, give {MI['trimmed20_equal_donor_mean_median']:+.6f}; {sum(v>0 for v in MI['trimmed20_per_donor_median'].values())}/{MI['n_donors']} donors are positive, whereas the original sign-only test was 9/12. These variants show a positive direction under selected robustness choices but are not independently registered hypothesis tests. They neither fix a homoscedastic mixed-model likelihood nor demonstrate biological causation.")
P.table(["A0 small/large", "prior blocks eligible", "donors", "positive medians", "equal-donor mean"],
        [[str(q['cutoffs'][0])+"/"+str(q['cutoffs'][1]),q['eligible_prior_blocks'],q['eligible_donors'],q['positive_donors'],f"{q['equal_donor_mean']:+.4f}"] for q in MI['boundary_sensitivity'].values()],
        "Fixed +/-10% cutoff shifts in the same original-block population; eligibility and sign counts may change. These are descriptive perturbations, not substitute discovery gates.")
P.p("The derived analysis CSV alone omits masks, detection confidence and excluded rows. On reviewing the accompanying source DIS database we found paired t=0/t=1 tabular particle rows, including blank endpoints and detection scores. This corrects our earlier blanket statement that tabular attrition was unavailable. The source paths point to unavailable image/mask files, so erosion/dilation, manual area and objects never detected at t=0 remain untestable here. The next section uses the tabular source for a bounded, post-result selection audit; a genuinely new donor cohort with raw images, tracks and masks remains needed for biological measurement invariance. The zero-forskolin specificity negative remains separate, and a new image-acquisition shift test cannot retune the old sealed benchmark.")
P.h("4.3+ Track attrition and selection-bias methodology (exploratory)", 2)
P.p("Sections 4.3h through 4.3p are the methodological spine of this paper: an explicit audit of how track attrition and selection can create or hide a size contrast. They are exploratory methods work, not independent replication of the size finding.")
P.h("4.3h Source-table track retention and reanalysis", 2)
TR = json.load(open("results/tabular_track_retention.json")); assert TR['n_original_eligible_blocks']==31
D={(z['arm'],z['size_bin']):z for z in TR['retention_by_arm_size']}
P.p(f"We inspected the released OrgaSegment DIS_database.csv, which stores {TR['n_source_paired_particles_all_arms_doses']:,} particle keys each with a t=0 and t=1 row, before the {TR['n_derived_particles_all_arms_doses']:,}-row derived table. Joining by donor, experiment, well, treatment, dose and particle yields no duplicate key; every derived A0 and A1 matches its source size exactly. For the 31 original eligible donor-by-experiment-by-dose blocks, among t=0 objects with positive baseline area, the small (<722 px) DMSO arm has {D['DMSO','small']['n_both_positive_area']}/{D['DMSO','small']['n_baseline_detected']} positive-area endpoints and {D['DMSO','small']['n_derived']} included derived objects. The small Trikafta arm has {D['VX445_VX661_VX770','small']['n_both_positive_area']}/{D['VX445_VX661_VX770','small']['n_baseline_detected']} endpoints and {D['VX445_VX661_VX770','small']['n_derived']} derived objects. In the large (>=1400 px) bins those counts are {D['DMSO','large']['n_both_positive_area']}/{D['DMSO','large']['n_baseline_detected']} and {D['DMSO','large']['n_derived']} for DMSO, versus {D['VX445_VX661_VX770','large']['n_both_positive_area']}/{D['VX445_VX661_VX770','large']['n_baseline_detected']} and {D['VX445_VX661_VX770','large']['n_derived']} for treatment. Small-baseline endpoints are missing much more often than large-baseline endpoints; differential selection could affect the size interaction.")
P.p(f"Using every source pair with positive baseline and endpoint areas, before the derived-table filter, all {TR['source_paired_positive_area']['n_blocks']} original blocks and {TR['source_paired_positive_area']['n_donors']} donors remain eligible. Their equal-donor mean of donor-median large-minus-small drug effects is {TR['source_paired_positive_area']['equal_donor_mean_median']:+.6f}, compared with {TR['derived_original']['equal_donor_mean_median']:+.6f} in the derived table; both have {TR['source_paired_positive_area']['positive_donors']}/{TR['source_paired_positive_area']['n_donors']} positive donors. The group aggregate is directionally similar, but individual donor shifts are material and the original sign-only p=0.073 is not repaired. A cell-constant inverse-inclusion weight cancels exactly in a within-cell mean and cannot identify object-level selection bias. Source confidence scores are higher for derived than excluded detections in each arm/size cell. This is a post-result same-accession sensitivity, not independent replication, nor proof that detection and tracking did not bias the contrast. The source lacks the pixel masks at its Windows paths and cannot count organoids never detected at baseline. Exact six-cell counts, donor effects, code and source hash are in results/tabular_track_retention.json and scripts/audit_tabular_track_retention.py.")

P.h("4.3i Post-result modeled inclusion among positive-area source pairs", 2)
CF=json.load(open("results/crossfit_inclusion_sensitivity.json")); assert CF['n_positive_source_pairs_in_original_small_large_cells']==4964
P.p(f"A simulated paper-text critique prompted a narrower selection diagnostic: can variables recorded for the {CF['n_positive_source_pairs_in_original_small_large_cells']:,} source pairs with both positive areas explain which {CF['n_derived_included_in_frame']:,} entered the derived analysis? We froze notes/postresult_inclusion_model_protocol.md after seeing the earlier size and retention results, then fit a five-fold experiment-held-out logistic inclusion model with baseline/end area, both detection scores, dose, donor, arm and well. Its out-of-fold AUROC is {CF['out_of_fold_inclusion_auc']:.4f} and Brier {CF['out_of_fold_inclusion_brier']:.5f} versus {CF['constant_prevalence_brier']:.5f} for constant inclusion prevalence. Crucially, inverse inclusion weights for included observations are proportional to 1/pi, not pi/prevalence as erroneously proposed by the external role-play response. Endpoint area and score make this an observed-selection diagnostic, not a pretreatment causal estimator.")
P.p(f"The weighted equal-donor mean of donor-median treatment-by-size effects is {CF['derived_ipw']['equal_donor_mean_median']:+.6f}, versus {CF['derived_unweighted']['equal_donor_mean_median']:+.6f} unweighted in the derived table and {CF['source_positive_pair']['equal_donor_mean_median']:+.6f} among all source positive-area pairs. The weighted minus derived shift is {CF['delta_weighted_minus_derived']:+.6f}, below the follow-up's 0.10 magnitude flag, but positive donor medians change from {CF['derived_unweighted']['positive_donors']}/12 to {CF['derived_ipw']['positive_donors']}/12, hitting the protocol's sign-count sensitivity flag. Of the source positive pairs, {CF['n_probabilities_clipped']} OOF probabilities were clipped to [0.05,0.995]; included weights have median {CF['included_weight_quantiles']['0.5']:.3f}, maximum {CF['included_weight_quantiles']['1']:.3f}, and effective sample size {CF['included_effective_sample_size']:.0f} of {CF['n_derived_included_in_frame']}. This small aggregate change cannot establish robustness to the major small-organoid missing-endpoint problem, because the outcome is unknown for those rows. A seeded 300-split DMSO pseudo-treatment control yields no absolute equal-donor effect as large as the real weighted contrast, but is a computational null only, not independent treatment evidence. No locked sign or specificity gate changes. Verbatim prompt/answer, corrected method and all per-donor/negative-control draws remain in results/judge_rounds/round07*, scripts/audit_crossfit_inclusion.py and results/crossfit_inclusion_sensitivity.json.")

P.h("4.3j Hypothetical missing-endpoint tipping curve", 2)
MT=json.load(open("results/missing_endpoint_tipping.json")); assert MT['n_original_blocks']==31
P.p(f"The preceding inclusion model cannot recover the large number of small baseline objects whose endpoint area is absent. We therefore froze a separate post-result, explicitly hypothetical tipping calculation on the same 31 original blocks (notes/postresult_missing_endpoint_tipping.md). Among t=0 detected objects in the four original arm/size cells, {MT['observed_reference_by_arm_size']['VX445_VX661_VX770:small']['n_missing']} of {MT['observed_reference_by_arm_size']['VX445_VX661_VX770:small']['n_missing']+MT['observed_reference_by_arm_size']['VX445_VX661_VX770:small']['n_observed']} small treated objects and {MT['observed_reference_by_arm_size']['DMSO:small']['n_missing']} of {MT['observed_reference_by_arm_size']['DMSO:small']['n_missing']+MT['observed_reference_by_arm_size']['DMSO:small']['n_observed']} small DMSO objects have no positive endpoint, versus {MT['observed_reference_by_arm_size']['VX445_VX661_VX770:large']['n_missing']} of {MT['observed_reference_by_arm_size']['VX445_VX661_VX770:large']['n_missing']+MT['observed_reference_by_arm_size']['VX445_VX661_VX770:large']['n_observed']} large treated and {MT['observed_reference_by_arm_size']['DMSO:large']['n_missing']} of {MT['observed_reference_by_arm_size']['DMSO:large']['n_missing']+MT['observed_reference_by_arm_size']['DMSO:large']['n_observed']} large DMSO. Assigning each missing endpoint its pooled observed arm-by-size median log swelling gives a reference equal-donor mean effect {MT['reference_delta0']['equal_donor_mean_median']:+.6f} and {MT['reference_delta0']['positive_donors']}/12 positive donor medians. This reference is an assumption, not an estimate of missing outcomes.")
P.p(f"In the adverse scenario we add a common delta only to hypothetical missing small-treated log responses, keeping the other three missing cells at their pooled observed medians. The equal-donor mean first reaches zero at delta={MT['first_nonpositive_mean']['delta']:.2f} log units (about {np.exp(MT['first_nonpositive_mean']['delta']):.2f}-fold higher imputed swelling for those missing treated-small objects); <=6 of 12 donor medians first remain positive at delta={MT['first_six_or_fewer_positive_donors']['delta']:.2f}, a {np.exp(MT['first_six_or_fewer_positive_donors']['delta']):.2f}-fold multiplier. At mean tipping, the imputed small-treated log response is {MT['observed_reference_by_arm_size']['VX445_VX661_VX770:small']['median']+MT['first_nonpositive_mean']['delta']:.3f}, below the observed small-treated 95th percentile {MT['observed_reference_by_arm_size']['VX445_VX661_VX770:small']['q95']:.3f}; that does not establish plausibility for the missing set, but does show this scenario need not exceed observed values. A negative-delta direction check increases the effect as expected. The whole curve and donor values are in results/missing_endpoint_tipping.json. These results are sensitivity to unidentifiable endpoints, not a new biological effect, imputed evidence, or a revised sign gate. Never-detected baseline objects remain outside even this analysis.")

P.h("4.3k Endpoint observability among baseline-detected organoids", 2)
EO=json.load(open("results/endpoint_observability.json"));assert EO['n_original_blocks']==31
P.p(f"The hypothetical tipping calculation does not estimate the probability that a baseline-detected object has an observed positive endpoint. We froze a further post-result source-table audit (notes/postresult_endpoint_observability_protocol.md) using the same original size bins and 31 blocks: {EO['n_baseline_detected_original_extreme_size_cells']} baseline objects, of which {EO['n_observable_endpoint']} have positive finite A1. A binomial model of endpoint observability includes treatment, large-size status, their interaction, log baseline area, baseline detector score and donor fixed effects, with donor-cluster uncertainty. The interaction coefficient is {EO['endpoint_model']['terms']['treatment:large']['coefficient']:+.4f} (odds ratio {EO['endpoint_model']['terms']['treatment:large']['odds_ratio']:.3f}, 95% cluster interval {EO['endpoint_model']['terms']['treatment:large']['ci95_logit'][0]:+.3f} to {EO['endpoint_model']['terms']['treatment:large']['ci95_logit'][1]:+.3f} on logit scale). Donor-equal medians of block-level treatment minus DMSO observability gaps are {EO['equal_donor_mean_median_gaps']['treatment_gap_small']:+.3f} for small and {EO['equal_donor_mean_median_gaps']['treatment_gap_large']:+.3f} for large, with large-minus-small gap {EO['equal_donor_mean_median_gaps']['gap_interaction_large_minus_small']:+.3f}. A baseline high-detector-score negative-control outcome has interaction odds ratio {EO['baseline_high_score_negative_control']['terms']['treatment:large']['odds_ratio']:.3f}, cluster p={EO['baseline_high_score_negative_control']['terms']['treatment:large']['cluster_p_two_sided']:.3f}; it does not remove unmeasured batch confounding.")
P.p("The treatment-versus-DMSO observability gap is relatively more favorable for small than large baseline objects in this conditional table, after accounting for specified covariates. This is a detection/retention result, not a directionally identified bias in swelling: the unknown A1 of absent endpoints can be high or low, and never-detected baseline objects are outside the frame. Donor and experiment effects were collinear in the restricted original-block design, so the final model uses donor fixed effects without experiment effects; only 12 donor clusters make model p-values descriptive. Per-block four-cell rates, per-donor medians, source hash and controls are in results/endpoint_observability.json and scripts/audit_endpoint_observability.py. Neither this model nor the previous hypothetical curve establishes biological invariance, independent replication or a changed sign gate.")

P.h("4.3l Missing-endpoint percentile envelope", 2)
MB=json.load(open("results/missing_endpoint_bounds.json"));assert MB['n_original_blocks']==31
P.p(f"The prior one-arm tipping scenario probes one missingness pattern; it does not bracket all possible missing outcomes. In a separate post-result partial-identification exercise, we assigned each missing A1 log response a value bounded by the pooled observed 5th and 95th percentiles in its treatment-by-size cell. Choosing the adverse endpoint separately in each of four cells gives an equal-donor mean of donor-median block contrasts as low as {MB['minimize_contrast']['equal_donor_mean_median']:+.3f}, with {MB['minimize_contrast']['positive_donors']}/12 positive donors; reversing every bound gives as high as {MB['maximize_contrast']['equal_donor_mean_median']:+.3f}, with {MB['maximize_contrast']['positive_donors']}/12 positive. Zero is inside this wide conditional envelope. By contrast, the observed-positive-pair mean is {MB['observed_positive_area_pair_mean']:+.3f} and derived-only mean {MB['derived_included_mean']:+.3f}. Even if unobserved responses were restricted to the broad observed central 90% of their own cells, the sign of the size-treatment contrast would not be identified by these source tables.")
P.p("The 5th-95th percentile assumption is hypothetical and chosen after seeing the source distribution. It is not proof that absent objects have outcomes inside those ranges; indeed, segmentation failures may select extremes. The calculations preserve all observed positive A1 outcomes and change only the missing entries in the 31 original blocks. Original masks and objects never detected at baseline remain unavailable. This honest envelope is not a new negative experimental result, a p-value or a license to rewrite the original biological gate. Full per-donor interval endpoints, source hash, protocol and script are in results/missing_endpoint_bounds.json, notes/postresult_missing_endpoint_bounds.md and scripts/bound_missing_endpoint_response.py.")

P.h("4.3m Well-balanced source-table check", 2)
WB=json.load(open("results/well_balance_sensitivity.json")); assert WB['n_original_blocks']==31 and WB['n_donors']==12
flag_word = 'triggered' if WB['predefined_sensitivity_flag'] else 'not triggered'
P.p(f"Source-table organoids from the same well are nested technical replicates. To test whether high-count wells dominate the observed-endpoint contrast, we fixed this post-result comparison before running it (notes/postresult_well_balance_protocol.md): take every positive-A0/positive-A1 source pair in the 31 original blocks and original size bins, compute each well's mean log(A1/A0), then weight observed wells equally within each arm-by-size cell. The original object-weighted mean is recovered at {WB['object_weighted']['equal_donor_mean_median']:+.6f} for the equal-donor mean of median block effects, with {WB['object_weighted']['positive_donors']}/12 positive donors. Equal-well weighting gives {WB['equal_well']['equal_donor_mean_median']:+.6f} with {WB['equal_well']['positive_donors']}/12 positive donors; its difference from object weighting is {WB['delta_equal_well_minus_object']:+.6f} log units. All 31 blocks retain four observed arm/size cells. The predeclared >=0.10 change, donor sign-count change, or missing-cell sensitivity flag was {flag_word}. This only shows little influence from unequal counts among wells with observed endpoints, not robust biological treatment-by-size response.")
P.p(f"We separately retained all {WB['n_baseline_detected']:,} baseline-detected objects, including the {WB['n_baseline_detected']-WB['n_observed_endpoint']:,} without a positive endpoint, to calculate the treatment-by-size interaction in endpoint-observability rates. Object-weighting gives {WB['observability_object_weighted']['equal_donor_mean_median']:+.4f} and equal-well weighting {WB['observability_equal_well']['equal_donor_mean_median']:+.4f}, both with {WB['observability_equal_well']['positive_donors']}/12 positive donor medians. Well-balanced observability does not identify the missing A1 values. Well counts vary from one to three in individual cells, and well grouping does not give independent donor replication. Per-cell counts, source hashes and executable code are in results/well_balance_sensitivity.json and scripts/audit_well_balance.py. This same-accession diagnostic cannot revise the failed original 9/12 sign gate or no-forskolin specificity gate; the prior missing-endpoint envelope still includes zero.")

P.h("4.3n Baseline image-edge proxy and endpoint visibility", 2)
ED=json.load(open("results/edge_observability.json"));assert ED['n_baseline_detected']==6964 and ED['n_blocks']==31
Nedge=sum(z['n'] for z in ED['edge50_counts'] if z['edge50'])
P.p(f"The source's baseline detection boxes permit a narrow test of whether proximity to the 2048-pixel image boundary predicts an absent endpoint. We locked a box-edge clearance below 50 pixels and a fixed 100-pixel sensitivity before calculating the association (notes/postresult_edge_observability_protocol.md). Among {ED['n_baseline_detected']:,} baseline-detected organoids in the original size/arm cells, {Nedge} have a box within 50 pixels of an edge. At 50 pixels, 24 of the 31 original blocks contain both edge and interior objects in all four arm/size cells, spanning 11 donors; the other seven blocks are listed without substitution in results/edge_observability.json. Thus the block-balanced edge comparison is underpowered and incomplete by construction.")
P.p(f"A binomial endpoint-observability model with baseline area, detector score, donor factors, treatment, size, treatment-by-size, edge, and edge-by-size/treatment terms converged. The edge main-effect odds ratio is {ED['edge50_model']['terms']['edge50']['odds_ratio']:.3f} (donor-cluster interval on logit scale {ED['edge50_model']['terms']['edge50']['ci95_logit'][0]:+.3f} to {ED['edge50_model']['terms']['edge50']['ci95_logit'][1]:+.3f}); edge-by-size odds ratio {ED['edge50_model']['terms']['edge50:large']['odds_ratio']:.3f} and edge-by-treatment odds ratio {ED['edge50_model']['terms']['edge50:treatment']['odds_ratio']:.3f}. Their wide intervals cover zero on logit scale. The treatment-by-size endpoint-observability odds ratio stays {ED['edge50_model']['terms']['treatment:large']['odds_ratio']:.3f}, versus {ED['edge100_model']['terms']['treatment:large']['odds_ratio']:.3f} under the fixed 100-pixel proxy. This is not evidence that boundary cropping is absent: box clearance does not show whether pixels or masks were cropped, and only 404/6,964 objects meet the 50-pixel indicator. Baseline-undetected organoids, absent endpoint areas and real image artifacts remain unmeasured. Full four-cell rates, excluded-block identities, hashed source and code are in results/edge_observability.json and scripts/audit_edge_observability.py. No original biological gate changes.")

P.h("4.3o Paired-center displacement and track-selection support", 2)
DS=json.load(open("results/displacement_sensitivity.json"));assert DS['n_original_blocks']==31 and DS['n_observed_pairs_original_extreme_cells']==4964
P.p(f"The source table stores paired x/y centers only when both segmented areas are positive. Before computing this post-result influence check, we fixed displacement bins [0,10), [10,25), [25,40), and [40,50] pixels, then a <40-pixel subset (notes/postresult_displacement_sensitivity.md). Across all {DS['n_positive_pairs_all_source']:,} positive-area pairs under the released source, the maximum displacement is {DS['max_displacement_all_source']:.1f} pixels and none exceed 50. This sharp support boundary suggests a pairing restriction, but source tracking code has not been established and the limit cannot be called its documented algorithm. For the {DS['n_observed_pairs_original_extreme_cells']:,} source-positive pairs in the original 31 size/arm blocks, the all-pair equal-donor mean of median block contrasts is {DS['all_positive_pair_score']['equal_donor_mean_median']:+.6f} ({DS['all_positive_pair_score']['positive_donors']}/12 positive donors). Restricting observed pairs to displacement <40 yields {DS['under40_score']['equal_donor_mean_median']:+.6f} with {DS['under40_score']['positive_donors']}/12 positive donors, but one fixed block loses the minimum three observations in one cell. Comparing exactly the same 30 retained blocks, all pairs give {DS['all_pairs_on_under40_eligible_blocks']['equal_donor_mean_median']:+.6f} and the <40 subset differs by {DS['delta_under40_minus_matched_all']:+.6f}. This is not a newly chosen eligibility gate or an independent replication.")
ctrl={z['size_bin']:z for z in DS['dmso_control']}
P.p(f"Within the DMSO-only observed pairs, the pooled displacement-response Spearman correlations are {ctrl['small']['pooled_spearman']:+.3f} for small and {ctrl['large']['pooled_spearman']:+.3f} for large; equal-donor means of within-well median correlations are {ctrl['small']['equal_donor_median_well_rho']:+.3f} and {ctrl['large']['equal_donor_median_well_rho']:+.3f}. The paired-center boundary is a reproducible property of these selected rows, not proof of correct track identity or absence of selection bias. Without x/y for 2,000 baseline-detected objects missing positive A1, no displacement analysis can recover their swelling responses or test whether long-moving organoids were excluded. The full four-cell displacement histograms, lost-block ID, source hash and script are in results/displacement_sensitivity.json and scripts/audit_displacement_sensitivity.py. This audit leaves the original failed donor sign and zero-forskolin gates untouched and cannot tighten the hypothetical missing-response envelope.")

P.h("4.3p Baseline neighbor spacing and endpoint observability", 2)
NO=json.load(open("results/postresult_neighbor_observability.json")); assert NO['n_eligible_baseline']==6964 and NO['n_observed_endpoint']==4964 and NO['n_eligible_fields']==143
NN=NO['model_terms']['log_nn1']; NI=NO['model_terms']['log_nn1:large']; NF=NO['within_field_close_minus_far']
P.p(f"A post-result, protocol-locked source-table audit measured the distance from each baseline-detected organoid to its nearest OTHER baseline-detected object in the same image, including middle-size objects in the neighbor pool (notes/postresult_neighbor_observability_protocol.md). The analysis frame is the original 31 positive-dose blocks and extreme-size cells: {NO['n_eligible_baseline']:,} baseline objects in {NO['n_eligible_fields']} source image paths, with {NO['n_observed_endpoint']:,} positive finite endpoints. Nearest-distance quartiles are {NO['nn1_quartile_cutoffs_px'][0]:.1f}, {NO['nn1_quartile_cutoffs_px'][1]:.1f}, and {NO['nn1_quartile_cutoffs_px'][2]:.1f} pixels. A descriptive donor-fixed-effects binomial model adjusted for arm, size, their interaction, baseline log area and detector score finds an odds ratio {NN['odds_ratio']:.3f} for a unit increase in log(1 + nearest distance) in small objects (12-donor cluster interval on the coefficient {NN['cluster_ci95_beta'][0]:+.3f} to {NN['cluster_ci95_beta'][1]:+.3f}); the size-by-distance ratio of odds ratios is {NI['odds_ratio']:.3f}, with an interval crossing 1. These are conditional object-level associations, not independent biological replicates.")
P.p(f"Within {NF['n_usable_fields']} image paths, the closer half's mean positive-endpoint rate minus the farther half's rate is {NF['observed_mean_field_difference']:+.3f}; 100 seeded within-field permutations of proximity labels give a descriptive 95% technical-null range [{NF['null_percentile_2_5_97_5'][0]:+.3f}, {NF['null_percentile_2_5_97_5'][1]:+.3f}]. Thus source baseline spacing marks an observability gradient, not a proven causal crowding or tracking mechanism. The source has paired t=0/t=1 keys even when A1 is blank; {NO['n_eligible_baseline']-NO['n_observed_endpoint']:,} endpoint responses remain unobserved, baseline-undetected organoids are outside the frame, and field paths may represent repeated acquisitions. Neither the failed original donor sign and zero-forskolin specificity gates nor the zero-crossing missing-response envelope changes. Full counts, source hash and code are in results/postresult_neighbor_observability.json and scripts/audit_postresult_neighbor_observability.py.")

NS=json.load(open("results/postresult_neighbor_stratified.json")); assert NS["n_frame"]==6964 and NS["n_positive_A1"]==4964
P.p(f"A retrospective composition check held image path, treatment arm and original size bin fixed before dividing each eligible stratum at its own nearest-neighbor median (notes/postresult_neighbor_stratified_sensitivity.md). {NS['n_eligible_strata']} strata cover {NS['n_eligible_objects']:,} of {NS['n_frame']:,} frame objects; {NS['n_excluded_strata']} sparse strata with {NS['n_excluded_objects']} objects are reported, not silently dropped. The closer-minus-farther endpoint-observability difference is {NS['equal_stratum_close_minus_far']:+.3f} with strata weighted equally, {NS['object_weighted_stratum_close_minus_far']:+.3f} when weighted by stratum object count and {NS['donor_equal_close_minus_far']:+.3f} when donor means are weighted equally ({NS['donor_negative_count']}/{NS['donor_n']} donor means negative). The largest cell-level gradient is the treated-small cell at {NS['per_arm_size'][2]['equal_stratum_close_minus_far']:+.3f}; all four cells' counts are recorded. One hundred seeded within-stratum proximity-label permutations yield a technical-null 95% range [{NS['null_95pct'][0]:+.3f}, {NS['null_95pct'][1]:+.3f}] for the equal-stratum mean. The gradient therefore is not solely the pooled arm/size composition in this source frame. It is still a retrospective association, not a treatment-effect correction, tracking-causation test, or an answer for the 2,000 absent positive endpoints.")

P.h("4.3q What the observed tracks identify, and what they do not", 2)
P.p("The selection audit has three distinct denominators. The released paired-particle table contains baseline-detected organoids, including those with a blank or nonpositive endpoint; the derived response table contains a narrower subset with usable A0 and A1; the actual well may also contain organoids never detected at baseline, which neither table enumerates. These populations cannot be treated as interchangeable. The following accounting table uses only the original 31 donor-by-experiment-by-dose blocks and four pre-existing treatment/size cells, not a newly selected donor population. Its observed-endpoint fractions condition on baseline detection; they are not treatment-response rates.")
_cells=[]
for _arm,_size in [("DMSO","small"),("DMSO","large"),("VX445_VX661_VX770","small"),("VX445_VX661_VX770","large")]:
    _z=MB['percentile_reference'][_arm+":"+_size]
    _n=_z['n_missing']+_z['n_observed']
    _cells.append(["Trikafta" if _arm!="DMSO" else "DMSO",_size,_n,_z['n_observed'],_z['n_missing'],f"{_z['n_observed']/_n:.3f}",f"[{_z['q05']:+.3f}, {_z['q95']:+.3f}]"])
P.table(["arm","baseline size","baseline detected","positive A1","missing A1","visible fraction","observed log-response 5-95%"],_cells,
        "Source-frame denominator audit for the four original cells. Percentiles describe selected positive-endpoint tracks and are assumptions, not known support bounds, for absent responses. Source: results/missing_endpoint_bounds.json and tabular_track_retention.json.")
P.p("The difference between small and large organoids is not a minor rounding issue. Conditional positive-endpoint visibility is about 45% in small DMSO and 46% in small treated objects, versus 89% and 83% in the respective large cells. The treatment gap in visibility reverses across these size strata; the donor-equal gap contrast is -0.093. A donor-adjusted observability model that also includes log baseline area and detector score estimates a treatment-by-large odds ratio 0.639. The baseline high-detector-score control has the analogous odds ratio 0.976. These are features of a selected imaging pipeline, not evidence that missing small objects would swell more or less. Donor fixed effects and twelve clusters do not remove plate, imaging, field, or track-assignment confounding.")
P.p("The original biological estimand is the donor-level drug-minus-DMSO contrast in large versus small log swelling. It differs from the estimand of the observability model, which is whether a baseline-detected object receives a positive A1. Conditioning the swelling analysis on a positive A1 may induce selection dependence between initial size and an unmeasured cause of endpoint detection. An inverse-inclusion model fitted among rows with both positive areas can address only selection into the derived table among those already observed; it cannot recover the 2,000 absent positive endpoints. This is why the small shift in cross-fit inclusion-weighted means is not a missing-data correction, and why a baseline score control does not identify unmeasured endpoint values.")
P.p("A separate spacing audit gives a measured example of this selection structure. Within image-path, treatment and size strata, closer baseline neighbors have an equal-stratum positive-endpoint rate 0.082 below farther neighbors; 279 eligible strata retain 6,952 of 6,964 source-frame objects and all twelve donor-equal differences have the same negative direction. Distance is a baseline geometry proxy, not a randomized crowding intervention: object density, image quality, local tracking competition and acquisition repetition can all generate it. A field-stratified association is informative about where data are lost, but does not estimate the missing swelling response or the biological size mechanism.")
P.h("4.3r Conditional missing-response bounds across donors", 2)
P.p("To expose how much the missing responses matter at the donor level rather than only in a pooled interval, The following table lists the fixed-donor contrast under three previously committed analyses. The first uses only source pairs with positive A0 and A1; the latter two impute each missing endpoint at an adverse or favorable observed cell-specific percentile. These endpoints are conditional scenario bounds, not formal identification bounds: the 5th and 95th percentiles of observed tracks need not contain the missing outcomes. The same donor, experiment, dose and size eligibility from the original blocked test is kept throughout.")
_md=TR['source_paired_positive_area']['donor_medians']
_lb=MB['minimize_contrast']['donor_medians'];_ub=MB['maximize_contrast']['donor_medians']
assert set(_md)==set(_lb)==set(_ub) and len(_md)==12
P.table(["donor","observed positive-pair contrast","adverse percentile scenario","favorable percentile scenario"],
        [[k,f"{_md[k]:+.3f}",f"{_lb[k]:+.3f}",f"{_ub[k]:+.3f}"] for k in sorted(_md)],
        "Fixed donor-level median block contrasts in log swelling. Scenario columns do not constitute confidence intervals or imputed outcomes; source: results/tabular_track_retention.json and results/missing_endpoint_bounds.json.")
P.p(f"Across donors, the equal-donor mean of observed positive-pair contrasts is {MB['observed_positive_area_pair_mean']:+.3f}; the adverse and favorable percentile scenarios span {MB['minimize_contrast']['equal_donor_mean_median']:+.3f} to {MB['maximize_contrast']['equal_donor_mean_median']:+.3f}. No donor remains positive in the adverse scenario; all twelve are positive under the favorable one. Because zero lies in this conditional envelope, source-table endpoint observability alone cannot establish the direction of the full baseline-detected-population effect. These scenarios alter only the missing A1 values; they cannot address organoids absent at A0, incorrect track links, size-dependent area error or the failed zero-forskolin specificity test.")
P.p(f"A less adversarial one-arm tipping exercise holds three missing cells at their observed medians and increases only the hypothetical small-treated missing log response. The equal-donor mean reaches zero after a +{MT['first_nonpositive_mean']['delta']:.2f} log-unit shift; this is a scenario threshold, not a fitted estimate of those missing responses. Its direction matters: greater swelling among missing treated-small objects raises the small-treatment effect and can erase the observed large-minus-small difference. The percentile envelope is broader because it permits all four cells to vary together. Neither calculation selects a correction to apply to clinical data.")
P.p("The next discriminating measurement is not another threshold tuned to this accession. A new acquisition should enumerate all baseline objects before filtering, retain raw masks at both time points, identify fields/plates/donors, audit possible one-to-many tracks, and include positive, zero-forskolin and CFTR-inhibitor controls. Lock the size function, donor-level endpoint, missingness model and analysis plan before new outcomes are opened. Until then, the substantive result of the present pipeline is a quantified selection problem and a within-accession size association with failed companion gates, not a verified CFTR-specific size mechanism or a patient predictor.")

P.h("4.4 Does size correction help theratyping? No", 2)
P.p("Locked before running: per-donor standardised Trikafta-vs-DMSO separation with raw, size-adjusted and size-band well readouts. Median "
    "raw 4.73, adjusted 4.35, band 4.39; adjusted better in 8 of 17 donors. Excluding organoids below 1,069 px (a post-hoc threshold): median "
    "4.73 -> 4.75, better in 9 of 17, Wilcoxon p = 0.68. The size trend does not improve calls, and the stricter matched-block H1 fails.")
P.h("Secondary exploratory analyses II: segmentation benchmark and split integrity", 2)
P.p("Sections 4.5 through 4.5j are a secondary arm. Its benchmark-beat claim is withdrawn because repeated acquisition fields cross the splits; the arm is retained as an integrity case study and labeled exploratory.")
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
    f"The filename-distinct {SA['split_counts']['eval']}-image eval split was scored separately, and its mAP@0.5 "
    f"remains {SA['heldout_eval_mAP50']:.3f} versus published {SA['published_mAP50']:.2f}. "
    f"Bootstrap across our eval images gives CI [{SA['eval_image_bootstrap_ci95'][0]:.3f}, "
    f"{SA['eval_image_bootstrap_ci95'][1]:.3f}] for our mean, but leader per-image predictions "
    "and guaranteed identical AP implementation are unavailable. No benchmark break "
    "(results/seg_integrity_audit.json).")

P.h("4.5b First clean train-only and val-only segmentation attempt: negative", 2)
TS = json.load(open("results/trainonly_seg/sealed_eval.json")); VSEL = json.load(open("results/trainonly_seg/val_selection.json"))
P.p(f"A separate protocol was frozen before fitting (results/preregistration_trainonly_seg.md): 184 training images only, "
    f"35 filename-distinct validation images, 12 evaluation files not scored until final selection, but repeated acquisition fields were later found across splits (Section 4.5e). Sixty epochs were fit; "
    f"epoch {VSEL['selected_epoch']} minimized validation cross-entropy. A fixed grid of {VSEL['grid_size']} "
    f"postprocessing choices used the first 20 validation images only, selecting {VSEL['selected_params']} "
    f"at validation AP50 {VSEL['selected_val_mAP50']:.3f}. This validation score is a selection statistic, not an independent benchmark.")
P.p(f"The single sealed final evaluation scored {TS['mean_author']:.5f} AP50 under the authors' released Cellpose-style scorer "
    f"and {TS['mean_ours']:.5f} under our unique-mask-count scorer, below the published 0.76. The only material "
    "matcher difference arises from one released GT mask with a skipped label ID: the original scorer uses the maximum ID "
    "as object count. This small numerical difference is not a model gain. Author baseline per-image predictions were not run in "
    "this 2-CPU, no-TensorFlow/GPU workspace, so we do not assert a paired statistical comparison. Clean selection fixes an "
    "older exact-file train/val overlap, but the later near-duplicate audit finds that it did not secure acquisition-group independence. It is not a benchmark match. Per-image TP/FP/FN and both AP "
    "versions are preserved in results/trainonly_seg/sealed_eval.json; selection frozen before eval at commit 120f944.")
P.table(["image #", "author AP50", "unique-count AP50", "GT masks", "predicted masks"],
        [[i, f"{d['ap_author']:.3f}", f"{d['ap_ours']:.3f}", d['n_gt'], d['n_pred']] for i, d in enumerate(TS['per_image'], 1)],
        "Sealed 12-image evaluation of the clean train-only model; image numbers follow the JSON per_image order, which has full filenames and matching counts. None selected the model or postprocessing.")

P.h("4.5c Locked v2 train-only improvement and narrow published numerical threshold", 2)
V2 = json.load(open("results/trainonly_seg/sealed_eval_v2.json")); V2S = json.load(open("results/trainonly_seg/val_selection_v2.json"))
assert len(V2['per_image']) == 12 and V2S['n_candidates']==108 and V2S['selected_params']==V2['fixed_params']
P.p(f"The negative clean v1 checkpoint above remains the historical result, not a deleted failure. Amendment-free preregistration in notes/prereg_ap50_beat.md locked the v2 gate at mean author-scorer AP50 > 0.76 on the same nominally held-out 12-file eval split. A later scene-overlap audit invalidates its independence (Section 4.5e). Fine-tuning trained solely on 184 train images for 30 extra epochs with stronger augmentation and a lower learning rate; the checkpoint at epoch {V2S['selected_epoch']} minimized cross-entropy on the filename-distinct but scene-overlapping validation split. The unchanged 108-combination grid was scored only on the first 20 validation images. Grid index {V2S['selected_grid_index']} selected {V2S['selected_params']} at validation AP50 {V2S['selected_val_mAP50']:.6f}; selection was frozen at commit {V2['selection_commit']} before one v2 evaluation-file score.")
P.p(f"The released author scorer (OrgaSegment v1.0.1 commit {V2['author_release_tag']}) gives an evaluation-file mean AP50 {V2['mean_author']:.6f}, exceeding the locked numerical comparison 0.760000 by {V2['mean_author'] - .76:.6f}. The unique-mask-count implementation gives {V2['mean_ours']:.6f} and is NOT used for the gate. These are the saved single-run per-image means in results/trainonly_seg/sealed_eval_v2.json, not replicate-run means or a confidence-bound comparison. The released Mask R-CNN predictions were not scored head-to-head here; the published rounded reference 0.76 +/- 0.12 may be an estimate from a differing evaluation protocol. The numerical margin is small relative to image-to-image spread. More importantly, repeated fields of view span train and eval files (Section 4.5e): the >0.76 numerical crossing is preserved as a historical calculation but cannot count as an independent benchmark win. The model cannot be called the best virtual organoid in the world.")
P.table(["image #", "author AP50", "unique-count AP50", "GT masks", "predicted masks"],
        [[i, f"{d['ap_author']:.3f}", f"{d['ap_ours']:.3f}", d['n_gt'], d['n_pred']] for i,d in enumerate(V2['per_image'],1)],
        "V2 single sealed author-scorer evaluation. Full image names, TP, FP, FN, and mask-ID-gap accounting: results/trainonly_seg/sealed_eval_v2.json. This is a frozen-selection evaluation-file comparison, but acquisition-field overlap disqualifies an independent held-out interpretation.")
P.p("Open controls: release the prediction masks and rerun both models under exactly identical code and split; estimate paired uncertainty over images and train seeds; test a separately sourced annotation cohort. No threshold, checkpoint, or postprocessing retuning after the sealed v2 result is permitted for a claim about this same eval split.")

P.h("4.5d What the thin numerical threshold crossing cannot establish", 2)
UA = json.load(open("results/ap50_uncertainty_audit.json")); assert UA['n_images']==12
P.p(f"After the sealed v2 result, we wrote and ran a descriptive uncertainty audit; it is explicitly post-evaluation and cannot change the predeclared >0.76 numerical gate. The 12-image author-scorer mean is {UA['author_scorer_v2_mean']:.6f} with sample image SD {UA['v2_image_sd_sample']:.6f}. A seed-{UA['bootstrap_seed']} image-resampling percentile interval from {UA['bootstrap_draws']:,} draws is [{UA['v2_image_bootstrap_ci95'][0]:.6f}, {UA['v2_image_bootstrap_ci95'][1]:.6f}]. It includes the published rounded 0.76; it does not account for acquisition-field leakage, and {sum(z['v2_mean']<=.76 for z in UA['leave_one_image_out'])}/12 leave-one-image-out v2 means fall at or below 0.76. All image deletion outcomes are preserved in results/ap50_uncertainty_audit.json; none are used to reselect a model. This is evidence of thin numerical margin; the cross-split repeated fields independently disqualify a held-out general-superiority claim.")
P.p(f"The paired v2 minus our own clean v1 image mean is {UA['v2_minus_clean_v1_paired_mean']:+.6f}, with conditional image-bootstrap interval [{UA['v2_minus_v1_image_bootstrap_ci95'][0]:+.6f}, {UA['v2_minus_v1_image_bootstrap_ci95'][1]:+.6f}]. It improves {UA['paired_improved_images']} images, declines {UA['paired_declined_images']} and ties {UA['paired_tied_images']} under a 1e-7 tolerance; the two-sided sign p over nonties is {UA['paired_nonzero_sign_p_two_sided']:.6f}. This within-project comparison is useful debugging evidence, not a head-to-head result against the released Mask R-CNN or a second cohort. Image bootstrap estimates uncertainty conditional on this fixed split and ignores training-seed and domain-shift variation. A new acquisition batch with locked weights/preprocessing and, ideally, paired author checkpoint predictions is the needed confirmatory benchmark.")

P.h("4.5e Post-result audit: repeated acquisition fields cross segmentation splits", 2)
XI = json.load(open("results/seg_split_integrity_v2.json"))
assert XI['split_counts']=={'train':184,'val':35,'eval':12} and XI['cross_split_identical_image_hashes']==0
P.p(f"A new post-result forensic scan found zero identical image or mask SHA256s across the {XI['split_counts']['train']}/{XI['split_counts']['val']}/{XI['split_counts']['eval']} train/val/eval files. But standardized 32-by-32 grayscale thumbnails flagged {len(XI['cross_split_near_duplicate_pairs'])} cross-split pairs above cosine {XI['near_duplicate_threshold_cosine']:.2f}. The top train/eval pair reached 0.996584. Direct visual inspection of the original image pixels confirmed the same organoid layout in this pair and in two train/val pairs despite different filenames, brightness and image detail. A fourth val/eval pair has a partly shifted but similar field; its acquisition identity is not adjudicated. This scan is post hoc and the threshold is an investigation aid, not a validated group identifier. The complete pair list and file hashes are in results/seg_split_integrity_v2.json; screenshots are not needed to reproduce the thumbnail audit from the original public files.")
P.p("This changes the interpretation of Sections 4.5a-d. Train-only means train-filename-only, not acquisition-scene-only. File-untouched means scores were not requested during selection, not that the underlying scenes were unseen. The 0.761108 author-scorer mean and the locked >0.76 numerical arithmetic remain exactly as recorded; the prerequisite of an independently held-out benchmark sample was not satisfied. We retract the benchmark-improvement claim rather than excluding overlapping eval images after seeing outcomes, rewriting the gate, or treating a post hoc subgroup as a new sealed test. The earlier 0.73867 negative remains in the record. We need source acquisition and donor metadata, group-disjoint training/validation, an untouched independently acquired evaluation cohort, and a paired author-baseline scorer before claiming benchmark improvement. Source-derived filenames are not enough to prove independence of remaining images.")

P.h("4.5f A scene-graph audit and its false positives", 2)
SG=json.load(open("results/scene_graph_forensics.json")); G=SG['threshold_sensitivity']
assert SG['image_count']==231 and SG['image_manifest_sha256']==XI['manifest_sha256']
P.p(f"We implemented an image-only graph for all {SG['image_count']} images, with every pair's 32-by-32 standardized thumbnail cosine and an exact-hash check. Connected components are candidate acquisition scenes; their transitive closure retains repeated frames. At threshold 0.95 there are {G[0]['n_components']} components, {G[0]['n_multiframe_components']} multi-image components, and {G[0]['n_cross_split_components']} cross-split components comprising {G[0]['cross_split_component_members']} images. At 0.97 and 0.99, {G[1]['n_cross_split_components']} and {G[2]['n_cross_split_components']} cross-split pairs persist. All cutoffs were selected or inspected AFTER the saved AP50 outcome and so are forensics, not preregistered evaluation choices (scripts/build_scene_graph.py, results/scene_graph_forensics.json).")
P.p(f"A seeded random sample of {SG['negative_control']['n_cross_split_random_pairs']} cross-split pairs has median cosine {SG['negative_control']['cosine_median']:.3f} and 99th percentile {SG['negative_control']['cosine_q99']:.3f}; this descriptive scale check does not yield a calibrated p-value. Visual inspection of the five flagged pairs at 0.95 confirmed three near-identical organoid layouts (one train/eval and two train/val). One val/eval pair may be a shifted scene and remains uncertain. The fifth, a 0003 train/val pair near 0.95065, shows visibly distinct layouts: it is likely a similarity false positive. These adjudications were based on the pixels, not model outcomes; the graph retains even the false-positive candidate to expose threshold limitations. No missing graph edge proves different donors, plates, wells or sessions, because that metadata is unavailable. The old test files have already been seen; excluding image pairs now cannot create a sealed holdout. A new group-disjoint acquisition cohort and published-model paired predictions are still required.")

P.h("4.5g Repeated-field labels confirm shared organoids", 2)
SL=json.load(open("results/scene_label_consistency.json")); assert len(SL['pairs'])==5
P.p(f"The three strongest cross-split image pairs, all detected above 0.99 cosine before opening their masks in this follow-up, also share labeled organoids. We matched raw 512-by-512 instance masks one-to-one at identical pixel coordinates with a maximum-IoU assignment. Respectively, the two train/val pairs contain {SL['pairs'][0]['one_to_one_matches_iou_ge_0_5']} and {SL['pairs'][1]['one_to_one_matches_iou_ge_0_5']} matches at IoU >=0.5 out of {SL['pairs'][0]['n_objects_a']}/{SL['pairs'][0]['n_objects_b']} and {SL['pairs'][1]['n_objects_a']}/{SL['pairs'][1]['n_objects_b']} annotated objects; the train/eval pair contains {SL['pairs'][2]['one_to_one_matches_iou_ge_0_5']} matches out of {SL['pairs'][2]['n_objects_a']}/{SL['pairs'][2]['n_objects_b']} objects. Corresponding foreground Dice values are {SL['pairs'][0]['foreground_dice_same_coordinates']:.3f}, {SL['pairs'][1]['foreground_dice_same_coordinates']:.3f}, and {SL['pairs'][2]['foreground_dice_same_coordinates']:.3f}. This is direct evidence that many identical organoid layouts and corresponding mask annotations cross filenames and splits, not merely similar illumination or texture.")
P.p(f"The visually distinct-layout 0003 train/val candidate at cosine 0.95065 yields {SL['pairs'][3]['one_to_one_matches_iou_ge_0_5']} matched objects and Dice {SL['pairs'][3]['foreground_dice_same_coordinates']:.3f}. A retrospectively selected 512-pixel train/eval distinct-layout comparison also yields {SL['pairs'][4]['one_to_one_matches_iou_ge_0_5']} matches; controls illustrate that thumbnail cosine alone makes false-positive candidate scene groups. The mask checks were performed only after the published-file AP50 was scored; they are annotation consistency evidence, NOT a prospective segmentation test or estimate of benchmark inflation. The same source masks may have been curated together and cannot establish donor independence, either. All original pairs and descriptive controls are in scripts/audit_scene_label_consistency.py and results/scene_label_consistency.json. We have not removed overlapped images or rescored a post hoc uncontaminated subset.")

P.h("4.5h Matched organoid area and boundary variability across repeated fields", 2)
BC=json.load(open("results/paired_scene_boundary_consistency.json")); assert [r['n_matched_ge_0_5'] for r in BC['pairs']]==[61,40,50,0,0]
P.p("The three high-similarity cross-split fields also permit a direct comparison of their separately stored instance boundaries and areas. A fixed post-result descriptive protocol was written for this analysis; it is not an independently held-out test (notes/paired_scene_boundary_protocol.md).")
P.p(f"At identical raw 512-pixel coordinates we matched labels one-to-one by IoU >=0.5. There are {BC['pairs'][0]['n_matched_ge_0_5']}, {BC['pairs'][1]['n_matched_ge_0_5']}, and {BC['pairs'][2]['n_matched_ge_0_5']} matched organoids in the three previously flagged pairs. Their median absolute log area ratios are {BC['pairs'][0]['median_abs_log_area_ratio']:.4f}, {BC['pairs'][1]['median_abs_log_area_ratio']:.4f}, and {BC['pairs'][2]['median_abs_log_area_ratio']:.4f}; the median one-pixel-tolerant boundary F1 values are {BC['pairs'][0]['median_boundary_f1_1px']:.4f}, {BC['pairs'][1]['median_boundary_f1_1px']:.4f}, and {BC['pairs'][2]['median_boundary_f1_1px']:.4f}. An exact-copy metric control returns zero area/centroid drift and boundary F1 1; two visually distinct-layout pair controls have zero matched objects at IoU >=0.5. The code and every per-object measurement are in scripts/audit_paired_scene_boundaries.py and results/paired_scene_boundary_consistency.json. This is a measurable repeated-field annotation/acquisition variability signal, not a human inter-rater experiment: acquisition time, annotator identity, mask-curation history and donor IDs are unknown. No valid AP50 ceiling or CFTR phenotype-robustness result follows from it. Independent image acquisition and blinded repeat annotation would be required for those claims.")

P.h("4.5i Retrospective overlap-exclusion score audit", 2)
SO=json.load(open("results/postresult_scene_overlap_score_audit.json")); assert SO['full_eval_n']==12 and abs(SO['full_eval_v2_mean']-0.7611081451177597)<1e-9
T95,T97,T99=SO['thresholds']; assert [t['thumbnail_cosine_threshold'] for t in (T95,T97,T99)]==[0.95,0.97,0.99]
assert T95['flagged_eval']==T97['flagged_eval'] and len(T95['flagged_eval'])==2 and len(T99['flagged_eval'])==1
P.p(f"A further post-result audit asks the narrowest question still available from already-saved outputs: if the graph-flagged cross-split eval images are removed after all scores were viewed, what happens to the saved author-scorer mean? Exclusions come only from the existing image-only graph, never from outcome values. At thresholds 0.95 and 0.97 the flagged eval images are the train-overlapping 20210727_11562 file (AP50 {T95['flagged_eval'][0]['v2_ap50']:.5f}) and the val-overlapping 20210727_26838 file (AP50 {T95['flagged_eval'][1]['v2_ap50']:.5f}); removing them leaves {T95['n_kept']} images with mean {T95['v2_kept_mean']:.6f}, versus {SO['full_eval_v2_mean']:.6f} on all 12 and {T95['v1_kept_mean']:.6f} for the earlier v1 model on the same retained images. At 0.99 only the train-overlapping file is flagged, leaving {T99['n_kept']} images at mean {T99['v2_kept_mean']:.6f}. Seeded image-bootstrap intervals, [{T95['kept_bootstrap_ci95'][0]:.6f}, {T95['kept_bootstrap_ci95'][1]:.6f}] and [{T99['kept_bootstrap_ci95'][0]:.6f}, {T99['kept_bootstrap_ci95'][1]:.6f}], stay wide and include the rounded published 0.76. Protocol, input hashes and every retained and flagged score are in notes/postresult_scene_overlap_score_audit.md, scripts/audit_postresult_scene_overlap_score.py and results/postresult_scene_overlap_score_audit.json.")
P.p("None of these subset means is a held-out benchmark score: the eval files, model and grid were all fixed and scored before the exclusion rule existed, unflagged files can still share acquisition sessions or donors, and no published author-model predictions exist for a paired comparison. The audit is reported because it removes a tempting misreading: even the most aggressive outcome-aware exclusion does not rescue the comparison, it lowers the retained mean below the rounded 0.76 point. The locked all-image v2 score stays exactly as recorded, the retained-subset numbers are diagnostics, and the requirement for a new group-disjoint acquisition cohort with paired author-baseline scoring is unchanged.")

P.h("4.5j Distinct mouse-intestinal detection transfer: negative", 2)
MX=json.load(open("results/postresult_mouse_yolo_domain.json")); assert len(MX['per_image'])==84 and MX['summary']['tp']==850
P.p(f"We tested the fixed human-organoid segmentation checkpoint on a different image and annotation domain, not as a repaired version of the original human mask benchmark. The Domènech-Moreno mouse-intestinal 4x EVOS release (Zenodo 6768583) contains 756 train and 84 released validation JPEGs with manually drawn YOLO object boxes across four morphology classes. Its 193.9 MB archive matches the source MD5. Before reading target scores we fixed the existing v2 U-Net weights and the original human-validation parameters, converted its instance masks to tight object boxes, collapsed the target classes into one object class and matched boxes once at IoU >=0.5. Across all 84 validation images and {MX['summary']['tp']+MX['summary']['fn']:,} valid labeled boxes, there are {MX['summary']['tp']} true positives, {MX['summary']['fp']} false positives and {MX['summary']['fn']:,} false negatives. Micro precision {MX['summary']['micro_precision']:.3f}, recall {MX['summary']['micro_recall']:.3f}, F1 {MX['summary']['micro_f1']:.3f}; mean per-image F1 {MX['summary']['macro_image_f1']:.3f}. One released zero-width annotation was recorded and excluded from the valid-box denominator. The source hashes, every image count, protocol and implementation are in results/postresult_mouse_yolo_domain.json, notes/postresult_mouse_yolo_domain_protocol.md and scripts/audit_mouse_yolo_domain.py.")
P.p("This is an external-domain negative for the current model's ready-to-use object detection. The target is mouse rather than human, class-specific loose boxes rather than pixel instance masks, and a different microscope/field density, so box F1 cannot be compared with the OrgaSegment paper's human mask AP50 0.76. No published comparator's predictions on the same 84 images are in hand, so this is NOT a leading-tool benchmark. No exact image-hash duplicates occur across target train/validation or original human source, but thumbnail similarity flags possible related target acquisition fields and lacks donor/session metadata. Because no mouse images trained our U-Net, this issue does not explain the low transfer recall, but it prevents an independent mouse split claim. Visual review of example low- and higher-scoring target images showed substantial contrast and object-density differences without identifying the cause of false negatives. The result narrows tool deployment; it neither tests CFTR swelling biology nor alters the failed original donor sign and specificity gates.")

P.h("4.5k Candidate-scene-group-disjoint reassignment and its exploratory test", 2)
GD=json.load(open("results/group_disjoint_seg/exploratory_eval.json"))
GS=json.load(open("results/group_disjoint_seg_split.json"))
assert GS['source_manifest_sha256']==GD['source_manifest_sha256'] and GS['split_counts']=={'train':180,'val':37,'eval':14}
assert GD['selected_epoch']==60 and len(GD['per_image'])==14 and GD['model_sha256']=='cc9f8ddfd07155854c8955ccd3089de777bd53fcf6131d0ef7aadb26142360f4'
P.p(f"We made the image-based leakage diagnosis operational, while preserving the original result as history. Under the 27 September exploratory protocol (notes/prereg_group_disjoint_seg.md, commit b8e5132), an exact image hash or standardized 32-by-32 thumbnail cosine >=0.95 joined images into candidate scene groups, including transitive links. Each group's reassignment followed the highest original priority eval > val > train. Five files moved, including one train file to eval and one val file to eval. The resulting 180/37/14 train/val/eval allocation includes all 231 released files exactly once and puts no graph component across splits (results/group_disjoint_seg_split.json, split implementation commit db30071). At cosine 0.97 and 0.99, the earlier forensic graph found {G[1]['n_cross_split_components']} and {G[2]['n_cross_split_components']} originally cross-split candidate groups, versus {G[0]['n_cross_split_components']} at 0.95; these are descriptive threshold checks, not alternative assignments chosen by AP50. The 0.95 graph is a proxy: one flagged pair has visually different layouts and unflagged files can still share donors, wells or sessions. Without acquisition metadata, this reassignment cannot certify independence.")
P.p(f"We freshly initialized the 3-class base-16 U-Net and trained on only the 180 reassigned train files, with no weights from the old v2 checkpoint. The frozen recipe used 512-pixel grayscale normalization, 256-pixel random crops, batch size two, weighted cross entropy 1/1/3, Adam at 0.002, seeded flips and rotations, and up to 60 epochs. The 37 reassigned validation files selected epoch {GD['selected_epoch']} by lowest weighted cross entropy checked every five epochs (loss {GD['selected_val_loss']:.9f}); their masks did not enter gradient updates. The checkpoint SHA256 is {GD['model_sha256']}. The old postprocessing grid was not rerun: script commit 6b23675 froze the code-default foreground/interior/boundary/minimum-size parameters {GD['params']} before the first reassigned-eval score, as permitted by that protocol. Model checkpoint and full epoch history are retained under results/group_disjoint_seg/.")
P.p(f"On the {len(GD['per_image'])} reassigned eval images, the pinned OrgaSegment v1.0.1 author AP50 implementation gives mean {GD['mean_author']:.6f}; our unique-mask-count implementation gives {GD['mean_ours']:.6f}. A seed-{GD['bootstrap_seed']} image bootstrap with {GD['bootstrap_draws']:,} draws produces a conditional percentile interval [{GD['image_bootstrap_ci95'][0]:.6f}, {GD['image_bootstrap_ci95'][1]:.6f}]. This includes the rounded published 0.76 and omits seed and donor/acquisition variance. The two scorers differ beyond float rounding on two images with missing ground-truth label IDs: the author scorer counts maximum instance label, while ours counts distinct labels. Neither metric is a head-to-head comparison with released Mask R-CNN predictions, which are unavailable here.")
P.table(["eval image #", "author AP50", "TP", "FP", "FN"],
        [[i, f"{d['ap_author']:.3f}", d['tp_author'], d['fp_author'], d['fn_author']] for i,d in enumerate(GD['per_image'],1)],
        "Candidate-group-disjoint exploratory reassignment: all 14 image scores, including the two moved-to-eval files. Full filenames, unique-count scores and mask-ID-gap accounting are in results/group_disjoint_seg/exploratory_eval.json.")
P.p("Failure analysis is heterogeneous rather than a single global threshold effect. The lowest-scoring images had author AP50 0.554 (36 TP, 15 FP, 14 FN) and 0.583 (49 TP, 8 FP, 27 FN); the highest reached 0.958 (68 TP, 1 FP, 2 FN). The two moved-to-eval images scored 0.850 and 0.767; they were original train/val files, not previously scored as eval, but both share image-similarity components with original eval files and were in the known source set. On the 12 original eval images shared with the old v2 report, a descriptive same-file paired mean change is -0.0204 AP50. This is not an isolated effect of leakage removal: weights, training membership, validation set, model selection and extraction parameters all changed. The original 0.761108 result and numerical threshold crossing remain archived; the new 0.750403 estimate is not a corrected benchmark score, an independent replication or proof that the old estimate was inflated by exactly their difference. The eval includes 10 singleton proxy groups and two two-image groups, but these are not verified independent donors or plates; the image-bootstrap interval ignores any hidden nesting. An untouched acquisition cohort with real donor/plate/field IDs and paired published-model predictions is still required for a leading-tool claim (notes/group_disjoint_eval_design_limits.md).")

P.h("4.5k.1 Same-file score reconciliation across the two models", 2)
SD=json.load(open("results/group_disjoint_seg/postresult_samefile_delta.json"))
assert SD['n_shared_original_eval']==12 and SD['n_new_moved_eval']==2 and len(SD['shared_images'])==12
P.p(f"The new 14-file eval mixes two previously unscored-as-eval images with the 12 original eval files. To keep denominators straight, a post-result reconciliation compares only those 12 identical filenames between the old v2 score and the freshly trained group-reassigned model. The older mean on these 12 files is {SD['shared_mean_old']:.6f}; the new mean on the same files is {SD['shared_mean_new']:.6f}, a new-minus-old mean of {SD['mean_new_minus_old']:+.6f}. The median per-file change is {SD['median_new_minus_old']:+.6f}, with range {SD['min_new_minus_old']:+.6f} to {SD['max_new_minus_old']:+.6f}. Three images improve, eight worsen and one ties exactly in saved float values. Old totals across the shared files are {SD['old_tp_fp_fn_sums'][0]} true positives, {SD['old_tp_fp_fn_sums'][1]} false positives and {SD['old_tp_fp_fn_sums'][2]} false negatives; new totals are {SD['new_tp_fp_fn_sums'][0]}, {SD['new_tp_fp_fn_sums'][1]} and {SD['new_tp_fp_fn_sums'][2]}. These totals describe instance matches, not independent experimental units.")
P.table(["shared original eval image #", "old AP50", "new AP50", "change"],
        [[i,f"{r['old_ap50']:.3f}",f"{r['new_ap50']:.3f}",f"{r['new_minus_old_ap50']:+.3f}"] for i,r in enumerate(SD['shared_images'],1)],
        "Same-file author-scorer AP50 changes across the two model fits. Full filenames, proxy components, ground-truth counts and TP/FP/FN are in results/group_disjoint_seg/postresult_samefile_delta.json.")
P.p("The largest per-file decline and the paired average are not measurements of how much leakage inflated the old score. The two fits also changed training membership, validation membership, weights, checkpoint selection and postprocessing. Both use the same masks for these 12 filenames and have unknown donor/plate nesting. The extra two files had been in original train/val and linked to an original eval image by the candidate scene graph; they do not turn this into a new sealed acquisition cohort. The reconciliation protocol was recorded after both outcomes (notes/postresult_samefile_seg_delta_protocol.md). Do not read the table as 12 independent confirmations or retune on its declines.")

P.h("4.5l Post-result candidate-scene-group uncertainty", 2)
GU=json.load(open("results/group_disjoint_seg/postresult_group_cluster_uncertainty.json"))
assert GU['n_groups']==12 and GU['n_images']==14 and GU['n_singleton_groups']==10 and GU['n_two_image_groups']==2
P.p(f"The 14-image bootstrap in Section 4.5k resamples image files as if they were exchangeable. A post-result sensitivity instead resamples the frozen 0.95 thumbnail-graph components. The reassigned eval files form {GU['n_groups']} candidate groups, {GU['n_singleton_groups']} singletons and {GU['n_two_image_groups']} two-image groups. Equal weighting by candidate group gives author AP50 {GU['equal_group_weight_ap50']:.6f}; seed-{GU['seed']} bootstrap of {GU['draws']:,} group-mean resamples gives a descriptive interval [{GU['equal_group_bootstrap_ci95'][0]:.6f}, {GU['equal_group_bootstrap_ci95'][1]:.6f}]. The image-weighted mean is {GU['file_weighted_ap50']:.6f} with image-bootstrap interval [{GU['file_weighted_image_bootstrap_ci95'][0]:.6f}, {GU['file_weighted_image_bootstrap_ci95'][1]:.6f}]. Different point estimates reflect changing the weights of the two duplicated-proxy groups; this is not a controlled variance-only adjustment or an estimate of donor-level uncertainty.")
P.p("This calculation was chosen and run after the AP50 outcomes were known (notes/postresult_group_cluster_uncertainty.md). The image-only graph can merge unrelated fields or miss common donors and plates, and both intervals omit uncertainty over training seeds and annotators. That the two intervals overlap 0.76 has no confirmatory role, because 12 original eval files had already been scored, two moved files were in the known source set, and the published model was not run head-to-head. Full group memberships, scores and deterministic bootstrap implementation are in results/group_disjoint_seg/postresult_group_cluster_uncertainty.json and scripts/audit_group_cluster_uncertainty.py.")

P.h("4.5m Why the two AP50 implementations differ", 2)
GAP=json.load(open("results/group_disjoint_seg/postresult_gt_id_gap_audit.json"))
assert GAP['n_images']==14 and GAP['n_gt_id_gap_images']==2 and GAP['n_other_images_with_difference_gt_1e_6']==0
P.p(f"The pinned released scorer reports mean {GAP['mean_author']:.6f}, while our alternate positive-label-count implementation reports {GAP['mean_unique_count']:.6f}, a difference of {GAP['mean_difference_unique_minus_author']:+.6f}. A post-result ledger check found exactly {GAP['n_gt_id_gap_images']} eval masks with a missing ground-truth label ID; these are the only images whose two AP50 values differ beyond 1e-6. The released function infers its ground-truth object count from the maximum label number, whereas ours counts positive IDs actually present. On the two affected images, the unique-count score is higher by {GAP['gap_images'][0]['difference_unique_minus_author']:+.6f} and {GAP['gap_images'][1]['difference_unique_minus_author']:+.6f}; their author-scorer TP/FP/FN triples are {GAP['gap_images'][0]['author_tp_fp_fn']} and {GAP['gap_images'][1]['author_tp_fp_fn']}. This is a label-index denominator issue, not two biological or model predictions. The author-scored value remains the reported primary exploratory number for compatibility with the released protocol; neither scorer version supplies a paired released-model benchmark. The per-image audit and implementation are results/group_disjoint_seg/postresult_gt_id_gap_audit.json and scripts/audit_seg_gt_id_gap.py.")

P.h("4.5n Object-density error accounting on the exploratory eval", 2)
OD=json.load(open("results/group_disjoint_seg/postresult_object_density.json"))
assert OD['low']['n_images']==OD['high']['n_images']==7 and OD['median_gt_cutoff']==64.5
P.p(f"The per-image AP50 spread motivates a narrower failure audit without another model run. After seeing the outcomes, we split all 14 reassigned eval files at their median ground-truth count ({OD['median_gt_cutoff']:.1f} instances per image): 7 images at or below it and 7 above. The low-count group's author AP50 mean is {OD['low']['mean_author_ap50']:.3f}, with {OD['low']['sum_tp']} TP, {OD['low']['sum_fp']} FP and {OD['low']['sum_fn']} FN; the high-count group's mean is {OD['high']['mean_author_ap50']:.3f}, with {OD['high']['sum_tp']} TP, {OD['high']['sum_fp']} FP and {OD['high']['sum_fn']} FN. Aggregate FN fraction rises from {OD['low']['aggregate_fn_fraction']:.3f} to {OD['high']['aggregate_fn_fraction']:.3f}; mean false positives per image rise from {OD['low']['mean_fp_per_image']:.1f} to {OD['high']['mean_fp_per_image']:.1f}. A continuous GT-count/AP50 Spearman coefficient is {OD['spearman_gt_count_vs_author_ap50']:+.3f}. Both moved-to-eval files land in the low-count group, so original split membership is also confounded with this contrast.")
P.p("This is an outcome-aware descriptive stratification, not a predeclared density subgroup or causal test. Ground-truth count is an annotation property and can covary with magnification, morphology, object overlap, illumination and acquisition batch, none of which was controlled here. The small 7/7 grouping does not establish a general density failure mode, does not identify the cause of false negatives, and cannot guide post-hoc threshold changes on this eval. Every filename, count and TP/FP/FN is in results/group_disjoint_seg/postresult_object_density.json; the rule and implementation are in notes/postresult_seg_object_density_protocol.md and scripts/audit_seg_object_density.py. Independent acquisitions with metadata and blinded masks would be needed to test a density interaction.")

P.h("4.5o Post-result strict-overlap AP75 diagnostic", 2)
T75=json.load(open("results/group_disjoint_seg/postresult_ap75_diagnostic.json"))
assert T75['n_images']==14 and abs(T75['ap50_reproduced_mean']-GD['mean_author'])<1e-7 and T75['checkpoint_sha256']==GD['model_sha256']
P.p(f"The released author metric at IoU 0.5 gives one view of detection, but an object whose overlap barely clears 0.5 can still have a visibly inaccurate boundary. After viewing the AP50 outcomes, we fixed an exploratory stricter-overlap check at IoU 0.75 (notes/postresult_ap75_boundary_protocol.md). The checkpoint, all 14 assigned eval image/mask pairs, three-view preprocessing and extraction thresholds were unchanged; rerunning the author AP50 metric matched all 14 previously saved per-image AP50 and TP/FP/FN values exactly. At IoU 0.75 the same predictions give mean AP {T75['ap75_mean']:.6f}, versus reproduced AP50 {T75['ap50_reproduced_mean']:.6f}: mean per-image AP75-minus-AP50 {T75['mean_ap75_minus_ap50']:+.6f}, median {T75['median_ap75_minus_ap50']:+.6f}. At 0.5 the images total {T75['tp_fp_fn_50_sums'][0]} matched instances, {T75['tp_fp_fn_50_sums'][1]} false positives and {T75['tp_fp_fn_50_sums'][2]} false negatives; at 0.75 the same predictions yield {T75['tp_fp_fn_75_sums'][0]}, {T75['tp_fp_fn_75_sums'][1]} and {T75['tp_fp_fn_75_sums'][2]}. These are alternate matching criteria applied to one prediction set, not new biological observations.")
P.p("The strict-IoU drop is a model-output diagnostic, not proof of one failure mechanism. Boundary errors, instance separation, annotation shifts and acquisition characteristics can all change matching at 0.75; the previously documented repeated-scene annotation/acquisition variability and two missing mask-label IDs complicate simple attribution. We have no paired released author-model AP75 predictions on these same files, so this cannot establish a strict-overlap benchmark win or loss. The 14 images are known-source files and the proxy scene graph is not donor or plate identity. Full per-image AP50/AP75 and TP/FP/FN, checkpoint/source hashes and the executable read-only scorer are in results/group_disjoint_seg/postresult_ap75_diagnostic.json and scripts/audit_group_disjoint_ap75.py. No retuning or original gate change follows.")

P.h("Secondary exploratory analyses III: organoid-to-tissue fidelity across public transcriptomes", 2)
P.p("Sections 4.6 through 4.13 are off-question secondary analyses of organoid-to-tissue fidelity in public expression data, retained and labeled exploratory.")
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
P.p("These series-level ranks were computed before the sample-level audit in Section 4.9. At this stage, the clean text-filtered subset appeared to match cultured fibroblasts often, and the working prediction was that a parenchymal-gene restriction would leave a fibroblast best-match share above 15% in that same text-filtered subset. That prediction is an internal same-series test, not a biological validation. We later checked sample metadata and found many series labeled as organoid-related do not contain organoid samples; Section 4.9 retracts the broad culture-fibroblast interpretation. The tool vorganoid fidelity produces a tissue-rank calculation, not a diagnosis or validated organoid-quality score.")

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
    f"p = {MC['CENT_vs_CENT_PARONLY']['p_exact']:.2g}). These same-series tissue-rank shifts suggest cell-composition sensitivity, but mixed sample labels prevent a general organoid-specific correction claim. Removing the selected culture genes did not materially change this rank result. "
    f"The internal text-filtered-subset prediction (fibroblast best-match share above 15% on parenchymal-only genes) was numerically met before the later sample-label audit: "
    f"{GC['SALL_PARONLY']['fibroblast_best_frac']:.0%} [{GC['SALL_PARONLY']['fibroblast_best_ci95'][0]:.2f}, {GC['SALL_PARONLY']['fibroblast_best_ci95'][1]:.2f}] and "
    f"{GC['CENT_PARONLY']['fibroblast_best_frac']:.0%} [{GC['CENT_PARONLY']['fibroblast_best_ci95'][0]:.2f}, {GC['CENT_PARONLY']['fibroblast_best_ci95'][1]:.2f}]. "
    "The lower bound of the second interval touches the threshold; more importantly, a threshold met on mixed samples cannot validate a general organoid phenomenon. The strict organoid-only reanalysis in Section 4.9 takes precedence over this earlier text-filtered aggregate. "
    "vorganoid fidelity --purity applies the gene restriction to a supplied count matrix without asserting clinical or biological calibration.")

P.h("4.8 What drives the text-filtered fibroblast rank pattern?", 2)
AG = _pd_gp = pd.read_csv("results/attractor_gprofiler.csv"); AG = AG[AG.term_size < 2000]
P.p("For the clean-subset series whose best match is cultured fibroblasts, each gene is scored by how much it pulls the profile toward "
    "fibroblasts and away from the organ of origin, averaged over series:")
P.equation("d_g = (1/|S|) sum_{s in S} z(p_sg) ( z(G_g,fib) - z(G_g,organ(s)) ),   z = rank-based standard score within a profile")
P.p("The top 300 genes were tested with g:Profiler (GO:BP, Reactome, KEGG; g:SCS correction; background = all scored genes).")
sel = pd.concat([AG.head(4), AG[AG.name.str.contains("immune|complement", case=False)].head(3)])
P.table(["source", "term", "adjusted p", "genes"], [[r.source, r.name, f"{r.p_adj:.1e}", int(r.intersection)] for r in sel.itertuples()],
        "Enrichment of attractor-driving genes (results/attractor_gprofiler.csv; terms with more than 2,000 genes omitted).")
P.p("Two components appear: mitotic cell-cycle genes (TOP2A, MKI67, CDK1) that cultures express and tissues do not, and immune and complement "
    "genes that tissues contain and these mixed series profiles lack. Removing either family alone leaves the earlier text-filtered rank pattern in place (Sections 4.6-4.7), but the later strictly organoid-only audit in Section 4.9 overrides a general attractor claim. "
    f"Restricting to {GM['n_protein_coding_genes']:,} HGNC protein-coding genes also retains that text-filtered rank pattern (fibroblast best match "
    f"{GC['SALL_PC']['fibroblast_best_frac']:.0%} [{GC['SALL_PC']['fibroblast_best_ci95'][0]:.2f}, {GC['SALL_PC']['fibroblast_best_ci95'][1]:.2f}] and "
    f"{GC['CENT_PC']['fibroblast_best_frac']:.0%} [{GC['CENT_PC']['fibroblast_best_ci95'][0]:.2f}, {GC['CENT_PC']['fibroblast_best_ci95'][1]:.2f}]) in the text-filtered subset; non-coding genes alone do not explain that aggregate, but sample composition still can. "
    "A future proliferation/purity hypothesis would need organoid and matched primary tissue from the same donors, profiled together. It cannot be inferred as the mechanism from these mixed series.")


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
P.p("A 2021 Cells study also tracked individual iPSC-derived human intestinal organoids every 15 minutes through a four-hour CFTR swelling assay, measuring size against baseline; it tested VX-809/VX-770 rescue in a single named F508del CF line. Its public supplemental archive contains figures, not raw object tracks, and its underlying data are available on request. Thus even individual tracking of intestinal CFTR FIS predates this work. It cannot establish or refute a donor-general drug-by-starting-size effect because it reports only one CF line and no public joinable object-level treatment table (results/ipsc_intestinal_source_suitability.md).")
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

P.h("5. Findings that did not hold", 2)
P.p("The original sign-only donor-block test (9/12, p=0.073) and zero-forskolin specificity test (3/11, p=0.967) failed; the later size mixed-model result is same-accession evidence, not donor-independent validation or a CFTR-specific mechanism. The initial power-law interpretation and size-adjusted donor discrimination claim were not retained. The first segmentation checkpoint fell below the published mean, and the v2 evaluation-file AP50 crossed its rounded comparison only after an image split later shown to share organoid fields; it is not an independent benchmark improvement. Full negative-arm results and protocols remain in the repository results and notes, and in Sections 4.3-4.5.")
P.p("Fidelity claims narrowed after their controls: the broad culture-fibroblast attractor was retracted for strictly organoid-only samples, kidney external profiles contradicted the original organ split, and several mechanism/clinical-target enrichment arms did not meet their registered bars. We retain those numerical tests in the repository rather than expanding each failed arm into the paper's main narrative. Neither the size result nor tissue-rank analysis is a clinical diagnostic.")
P.h("6. Discussion")
P.p("The same-accession eligibility sensitivity weakens from 9/12 (p=0.073) at three measured objects per cell to 8/12 (p=0.194) at five. At ten, the apparently low p comes from only six surviving donors and is uninterpretable under the frozen minimum-eight-donor rule. All three originally nonpositive donors are excluded at minimum ten, so the six-of-six result selects a different donor subset; excluded donors are not imputed as negatives. Selective eligibility is another reason not to claim a donor-general effect (results/finding_block_eligibility_sensitivity.md; results/finding_eligibility_composition.md).")
P.p("A frozen paired no-forskolin specificity test also fails: 3/11 donor-median contrasts are more positive at 0.128 uM than zero forskolin (one-sided p=0.9673). This weakens the proposed CFTR-specific interpretation without proving its absence or explaining the zero-dose response (results/finding_zero_fsk_control.md).")
P.p("In this accession, drug-minus-DMSO swelling tends to rise from small to larger organoids and then plateau. "
    "The original plate-and-dose-matched donor sign test fails at alpha 0.05, while the later locked magnitude-aware donor test and mixed model detect a same-accession positive interaction. This conflict of estimator sensitivity does not establish an independently replicated or CFTR-specific donor-general effect. Uniform surface-flux geometry predicts the reverse direction, but that simple model is not a full "
    "measurement-bias test. Lumen maturity and cell composition are hypotheses, not measured explanations. "
    "Falsifiable prediction: in an independent, matched, time-lapse single-organoid intestinal FIS cohort, with "
    "size thresholds locked before analysis, a donor-level modulator-minus-DMSO large-minus-small effect will "
    "be positive in most donors, with exact donor-magnitude and mixed-model tests agreeing in a new cohort. This needs a second suitable dataset and mechanistic "
    "measurements; the public leads audited here do not supply those controls.")
P.p("Practical note. In this dataset, starting size is associated with some exploratory single-organoid readouts, but the registered matched-block and specificity tests fail, and size correction did not improve per-donor discrimination. Recording the starting-size distribution across passages is a possible research quality-control variable, not a tested intervention or clinical adjustment.")
P.p("Independent novelty and mechanism remain unproved. We did not find the exact intestinal drug-by-starting-size interaction in the bounded literature search, but Calucho et al. 2021 tested starting size versus individual CFTR nasospheroid FSK response and found no association. Berical et al. 2022 tracked individual iPSC-derived airway spheroid FIS and tested modulators across multiple patient lines; the released source workbook does not link patient-treatment outcomes to per-object baseline areas. A 2021 Cells study had already tracked individual iPSC-derived intestinal CFTR-swelling objects in one CF line; its public supplement lacks paired raw object tracks. A nasal 2D-derived-organoid study reported swelling variation in large structures (bioRxiv 2021.07.20.453105). The narrower intestinal question is unconfirmed and this was not a systematic review.")
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
          "High-Throughput Functional Analysis of CFTR and Other Apically Localized Proteins in iPSC-Derived Human Intestinal Organoids. Cells 2021. https://www.mdpi.com/2073-4409/10/12/3419 . Supplement suitability: results/ipsc_intestinal_source_suitability.md",
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
