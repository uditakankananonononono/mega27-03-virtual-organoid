# Candidate finding: size-dependent CFTR swelling (item 3, vorganoid)
Data: OrgaSegment DIS dataset (Lefferts et al. 2024, Zenodo 10610438), 17 CF donors, per-organoid area at t0/t1.
Test: within-well slope of log(fold change) on log(initial area); within-well centering removes donor/plate/well effects; 95% CI by well-cluster bootstrap (500).
Forskolin > 0 only.
| condition | organoids | wells | slope | 95% CI | donors with slope>0 | mean fold |
|---|---|---|---|---|---|---|
| DMSO | 3810 | 94 | 0.011 | [-0.002, 0.028] | 12/17 | 1.43 |
| VX770 | 1496 | 30 | 0.226 | [0.185, 0.269] | 5/5 | 2.14 |
| VX661+VX770 | 4493 | 112 | 0.042 | [0.028, 0.056] | 12/16 | 1.70 |
| VX445+VX661+VX770 | 4984 | 121 | 0.116 | [0.090, 0.142] | 14/17 | 2.44 |
Interpretation: surface-limited secretion (alpha=2/3) predicts slope < 0; observed slope > 0 whenever CFTR is activated and ~0 in DMSO. Area-measurement noise biases slope negative (regression dilution in the denominator), so the positive sign is conservative. Implication: raw fold-change readouts are confounded by organoid size; a size-adjusted readout is needed.
Status: candidate, not yet confirmed. Open checks: segmentation-error model, lumen-presence threshold, replication on FIS time series.

## Follow-up test (locked before running): does a size-adjusted readout improve Trikafta-vs-DMSO separation?
Well-level readouts: raw mean log fold; per-well regression predicted at global median size (adj); size-band mean (band).
Per-donor standardized separation (Trikafta vs DMSO wells), 17 donors: median raw 4.73, adj 4.35, band 4.39; adj > raw in 8/17.
Result: NEGATIVE. Size adjustment does not improve donor-level modulator discrimination; the size effect is real but does not change theratyping calls at this assay's well counts.

## Replication check on FIS set: not possible
FIS_database.csv holds one area value per well per time point (868 rows, 7 time points, one position per well), i.e. summed organoid area, not single organoids. It cannot test per-organoid size dependence. Replication needs another single-organoid dataset (candidate: OrganoID time-lapse data, Matthews et al. 2022).

## Revision (8:43 PM): the size effect is a threshold/saturation, not a power law
Checks: within-well slope restricted to larger organoids turns negative (Trikafta A0>2000 px: -0.117; DMSO: -0.062), and quintile means rise then plateau. Simulated segmentation noise on a size-invariant truth gives slopes of only -0.002 to -0.055 (sd 5-20%), so noise can't create the positive effect.
DMSO-controlled, donor-matched curve (drug minus DMSO mean log fold per size octile; CI = donor bootstrap, 300):
- Trikafta (17 donors): 0.336 [0.262, 0.405] in the smallest octile (<722 px) rising to a plateau ~0.55-0.64 above ~1,400 px; top minus bottom 0.251 [0.184, 0.323].
- VX661+VX770 (16): 0.136 -> ~0.28-0.31 plateau; top minus bottom 0.169 [0.081, 0.259].
- VX770 (5): -0.153 -> 0.31-0.42; top minus bottom 0.462 [0.335, 0.572].
Named candidate: "small-organoid attenuation" - CFTR-modulator-induced swelling is suppressed in the smallest organoids and saturates above a size threshold (~1,000-1,400 px in these images). This replaces the earlier alpha>1 power-law reading, which the large-organoid checks contradict.
Follow-up (post hoc threshold, so exploratory only): excluding organoids <1,069 px does not improve per-donor Trikafta-vs-DMSO separation (median 4.73 -> 4.75, better in 9/17, Wilcoxon p=0.68). Negative for the clinical readout, again.

## Per-donor attenuation (Trikafta minus DMSO, smallest vs largest quartile of starting size)
Attenuation (large-quartile effect minus small-quartile effect) > 0 in 13/17 donors; it does not track overall response size (Spearman rho 0.24, p = 0.35). Exceptions: 2 of 3 X/X donors and one F508del/R117H donor.
