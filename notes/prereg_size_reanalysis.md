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
