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
