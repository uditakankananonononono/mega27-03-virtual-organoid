# Verdict-01 restructure plan (queue items A, B, C) - 27 September 2026

Status: framing amendment after owner-provided verdict 01 (wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMDY0OEIwREIxQjM5MjU3MDQ5MgA=). Original pre-registered gates are unchanged. This is presentation re-framing plus newly pre-registered follow-up work, never a gate redefinition.

## One primary question (locked from verdict item 1 / reframe 1)
"Does starting organoid size modify the CFTR-modulator swelling response of single intestinal organoids after accounting for track attrition?"

## One primary endpoint + one primary analysis (verdict item 6)
- Endpoint: donor-level large-minus-small Trikafta-minus-DMSO log-swelling contrast (the pre-registered donor-block analysis already locked before scoring: 31 donor x experiment x dose blocks, 12 donors).
- Primary analysis: the locked mixed model + 4,096-assignment exact donor-level sign flip, reported WITH its companion pre-registered primary tests that failed (sign-only 9/12 p=0.073; zero-forskolin specificity 3/11 p=0.967). No rescue ordering: failures are stated first.
- Everything else in the paper is labeled exploratory, in text.

## Reframe
Verdict reframe 1 (methodological CBIO: size-aware, missing-data-robust pipeline). The attrition/selection-bias machinery (endpoint observability GLM, neighbor observability + stratified, track retention, well-balance, missing-endpoint bounds/tipping) moves from scattered audits to the methodological spine. This is the verdict's strongest-angle note (top-5 #4).

## Body-page rule reconciliation (owner rule vs verdict item 1/20)
Owner rule: 50+ substantive TEXT body pages excluding headings/appendix/references/figure-only pages. Current ~41. Verdict wants off-question arms moved out of the primary story. Reconciliation: off-question arms (segmentation benchmark, drug-ranking, GEO/organ fidelity) are NOT deleted and NOT pushed to an uncounted appendix; they are consolidated into clearly-marked "Secondary exploratory analyses" body sections after the primary story, each prefixed with an exploratory label. The primary-question material (question, data, attrition methodology, size results, limitations, future experiment) is EXPANDED with existing per-donor, per-block and sensitivity content consolidated from notes into the body. Net body pages must not decrease; check at rebuild. Flagged to parent as a judgment call.

## Section map (old -> new)
- Title/abstract: rewritten (this wake).
- New 1.1: primary question, endpoint, analysis, amendment notice (this wake).
- Sections 2-4.3x (biophysics, data, size results, attrition audits): remain the primary spine; attrition/observability subsections consolidated under one "Track attrition and selection-bias methodology" body section (next wake).
- Segmentation benchmark sections -> "Secondary exploratory analyses I: segmentation integrity" (next wake); benchmark-beat claim already withdrawn at 4c5ae30.
- Drug ranking (Drevinek, VTI/ETI) -> "Secondary exploratory analyses II" (next wake).
- GEO/organ fidelity arms -> "Secondary exploratory analyses III" (next wake).
- Item 19 language pass: "in vitro computational reanalysis only, not a diagnostic or treatment predictor" in abstract + limitations; no clinical framing in headings (this wake in abstract).
- Item C (integrity): group-disjoint segmentation split, pre-registered protocol FIRST (notes/prereg_group_disjoint_seg.md), then implementation; labeled new exploratory evaluation, not a benchmark-gate repair (next wakes).

## 12-slide story
notes/verdict01_isef_story_12slides.md (this wake).
