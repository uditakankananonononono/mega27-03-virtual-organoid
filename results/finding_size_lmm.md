# Finding: donor-level size-dependence of Trikafta FIS response (VO3-R2, amendment A1)
Date: 2026-09-26. Status: DISCOVERY under locked VO3-R2 + A1 criteria; pending judge loop + external replication hunt.

Result (results/size_lmm.json; input csv sha256 a32e5cf8...):
- 31 eligible donor x experiment x dose blocks, 12 donors, 124 block-cell rows
  (construction byte-identical to locked src/blocked_size.py assess(); cutoffs 722/1400, >=3/cell).
- LMM on block-cell mean log swelling, ~ treatment*large + (1|donor) + (1|experiment), ML:
  interaction beta = +0.253 (larger drug response in large organoids), LRT stat 17.92,
  one-sided p = 1.15e-05 (both fits converged).
- Exact donor-level sign-flip (full 2^12 = 4096 enumeration on per-donor median block
  effect, the same donor statistic as the locked sign test but magnitude-aware):
  observed mean of donor medians 0.195, exact one-sided p = 0.00293.
- Both locked criteria agree in the predicted (positive) direction => verdict DISCOVERY.
Interpretation: within donor, larger starting organoids show a LARGER Trikafta swelling
response than small ones (small-organoid response is attenuated). The earlier preregistered
sign test (9/12 positive, p=0.073) failed only on power - it discards magnitudes; the
magnitude-aware exact donor test and the LMM both clear p<0.005.
Caveats (preserved): (1) size cutoffs were selected on this same accession (selection
leakage - carried from commit 87ec71a's warning); (2) cell-mean LMM treats cell means as
homoscedastic though cell n varies 3..99 - verdict does not rest on the LMM alone because
the distribution-free donor-level test agrees; (3) same-accession analysis is not
independent replication - the external multi-donor hunt (results/external_fis_source_audit.md)
remains open with its own mini-prereg requirement; (4) assay-sensitivity result, not a
treatment recommendation.
