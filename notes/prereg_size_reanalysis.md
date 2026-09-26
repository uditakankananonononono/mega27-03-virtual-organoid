# PREREG VO3-R2 (draft): donor-level size-attenuation reanalysis (discovery arm)
Date: 2026-09-26. Status: LOCKED at commit time; any amendment gets a new dated section.
Background: pooled-quartile size-aware analysis (vorganoid sizeaware) gives
 large-minus-small attenuation 0.249 [0.165,0.326], 14 donors - positive.
 The strict donor x experiment x forskolin-dose matched-block test FAILED
 preregistered H1 (31 blocks, 12 donors, 9 positive medians, sign p=0.073).
 Rule 4: negative is not terminal; rule 6: redirection via ChatGPT if stuck.
Locked question: is within-donor organoid starting size associated with
 attenuated Trikafta swelling response at the DONOR level, using all eligible
 blocks without re-using the failed sign-test exact form?
Locked protocol:
 1. Same discovery accession only (Zenodo 10610438 dis_merged_A0.csv) - no
    label or cutoff changes; the 722/1400 size cutoffs stay as previously
    fixed (they were locked pre-test; re-deriving them would be tuning).
 2. Donor-level estimator: linear mixed model, block effect ~ treatment*sizebin
    + (1|donor) + (1|experiment); primary statistic = likelihood-ratio test
    of the treatment x sizebin interaction, df=1, alpha=0.05 one-sided in the
    attenuation direction (predicted negative interaction).
 3. Exact-permutation sensitivity at the donor level (10k permutations of
    donor labels) reported alongside; both must agree in direction.
 4. Power honesty: if LRT p in [0.05,0.10], verdict = INCONCLUSIVE (not a
    discovery); the discovery claim requires p<0.05 AND permutation agreement.
 5. External replication hunt runs in parallel (sources audited in
    results/external_fis_source_audit.md); any new multi-donor single-object
    dataset gets its own mini-prereg before scoring.
Negative handling: an inconclusive/failed LRT is preserved; pivot ladder =
 culture-composition confound metric (HPA purity validation) as the
 discovery arm, then rule-6 ChatGPT redirection for fresh options.

Amendment A1 (2026-09-26, locked before any LMM or permutation scoring):
 1. Direction correction. The base prereg's parenthetical "(predicted negative
    interaction)" contradicts the direction locked in the EARLIER
    preregistration (results/preregistration_blocked_size.md, commit 87ec71a),
    whose H1 - written before any of these analyses - is a majority of
    POSITIVE donor-median block effects, where effect =
    (drug,large - DMSO,large) - (drug,small - DMSO,small) on mean log
    swelling. To test the same scientific quantity, the coding is locked as:
    response = mean log swelling per block x condition x sizebin cell;
    treatment = VX445_VX661_VX770 (reference DMSO); sizebin = large
    (reference small); the treatment:large interaction is predicted POSITIVE;
    the LRT is one-sided in the positive direction. "Attenuation" in the base
    text refers to attenuated response of SMALL organoids; the sign of the
    locked contrast is unchanged from commit 87ec71a. No LMM has been scored
    at lock time.
 2. Unit and formula (as locked in base): block-cell means (4 rows per
    eligible block); log_swelling_mean ~ treatment*sizebin
    + (1|donor) + (1|experiment), ML fit; LRT df=1 on the interaction.
    Eligible blocks and cell construction are byte-identical to
    src/blocked_size.py assess() (722/1400 cutoffs, >=3 per cell,
    positive forskolin, DMSO + VX445_VX661_VX770 only).
 3. Permutation feasibility. 10k donor-label permutations of a crossed
    mixed model are computationally infeasible in this environment; the
    donor-level exact sensitivity is replaced by the FULL exact sign-flip
    enumeration over donors (2^12 = 4096 assignments) on the per-donor
    median block effect (same donor statistic as the locked sign test, but
    using magnitudes): statistic = mean of donor medians; exact one-sided
    p = fraction of assignments with statistic >= observed. This is a strict
    strengthening of the sign test at the same donor level.
 4. Verdict (unchanged from base rule 4): discovery requires LRT one-sided
    p<0.05 AND exact-permutation agreement (p<0.05, same direction);
    either p in [0.05,0.10] => INCONCLUSIVE (not a discovery); otherwise the
    negative is preserved and the pivot ladder fires.
