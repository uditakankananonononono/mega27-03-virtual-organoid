# VO3 provided-verdict 01 - amendment queue, LOCKED 27 September 2026 ~11:14 IST before any execution

Queue discipline: this queue is locked before any verdict-driven work starts. Items execute in the order listed. No gate is redefined after outcomes: the original pre-registered gates (donor-level sign test 9/12 p=.073 FAIL; zero-forskolin specificity 3/11 p=.967 FAIL; sealed segmentation benchmark claim WITHDRAWN for scene overlap) stay in the record unchanged. The verdict is an ISEF-competitiveness critique; executing it changes paper framing, structure and presentation, and adds new pre-registered work - it does not convert failed gates into passes and does not erase preserved negatives. Classification: OPEN-EXEC = doable now computationally/editorially; PARTIAL = the 72pp build 4c5ae30 already covers part of the item, usually as its own stated limitation (no double-claiming the existing coverage as new work); OPEN-EXTERNAL = needs data we do not have; PROCESS = owner-level decision.

## The 20 weaknesses, mapped vs build 4c5ae30

1. No single clear research question - OPEN-EXEC. Restructure around one primary question: "Does starting organoid size modify CFTR-modulator swelling response after accounting for track attrition?" Move segmentation, drug-ranking, GEO, organ-fidelity arms to appendix/backup. (Cheapest first deliverable.)
2. AI research agent / authorship - PROCESS. Title-page line "computational work by an AI research agent" is accurate today. Any student-led recast and AI-disclosure wording is the owner's decision; we prepare options, we do not choose her authorship statement.
3. Novelty weak, prior art acknowledged - OPEN-EXEC. Add a systematic prior-art subsection (Calucho 2021, Berical 2022, 2021 Cells study, OrgaSegment, DIS source papers) and state the narrow gap: donor-matched modulator-by-size interaction with pre-registered cutoffs plus an explicit track-attrition model.
4. No independent replication - PARTIAL. Paper already says every result is within-accession and not independent validation (abstract, size-LMM sections); external FIS missing-tracks bound done (results/external_fis_missingtracks.json). A truly independent individual-organoid dataset is OPEN-EXTERNAL.
5. Key pre-registered tests fail - PARTIAL. Paper already reports sign 9/12 p=.073, zero-fsk 3/11 p=.967, VTI-vs-ETI p=.254 as failures. Verdict asks to reframe as exploratory/negative-first rather than leading with the mixed-model rescue: folded into item 1 restructure; the mixed model becomes secondary, clearly post-lock.
6. Multiple comparisons / forking paths - OPEN-EXEC (with item 1). Freeze one primary endpoint (donor-level large-vs-small log-swelling contrast) and one primary analysis (the pre-registered donor-block mixed model + exact sign flip, as already locked); label every other analysis exploratory in the text.
7. Cutoffs chosen on same accession - PARTIAL. Paper already states 722/1400 px cutoffs were data-derived on this accession and reports octile (near-continuous) results. External-cohort locking is OPEN-EXTERNAL; continuous size modeling extension is OPEN-EXEC (exploratory label).
8. Selection / attrition bias unresolved - PARTIAL (strongest existing coverage). Endpoint-observability GLM (results/endpoint_observability.json), neighbor observability + stratified (postresult_neighbor_observability.json, postresult_neighbor_stratified.json), tabular track retention, well-balance sensitivity, missing-endpoint bounds and tipping analyses all exist. Raw masks/never-detected objects OPEN-EXTERNAL.
9. Measurement error in A0/A1 - PARTIAL. Measurement-invariance exploratory check (scripts/measure_size_invariance.py, results/measurement_invariance.json), displacement sensitivity, paired-scene boundary consistency exist. Formal errors-in-variables modeling OPEN-EXEC (exploratory).
10. Segmentation split leakage - PARTIAL. Scene overlap documented, benchmark-beat WITHDRAWN (notes/seg_split_scene_overlap_correction.md), scene-graph forensics + overlap score audit done. Remaining: acquisition-group-disjoint retrain + rescore - OPEN-EXEC, the biggest integrity item (new pre-registered split; not a repair of the withdrawn claim).
11. Segmentation model not SOTA / thin margin - PARTIAL. Paper already reports sd ~0.13, 12 images, rounded 0.76 reference, AP50 uncertainty audit; identical author scorer used. New independent annotation cohort OPEN-EXTERNAL.
12. No biological mechanism - PARTIAL. Paper already states mechanism untested; falsification conditions for future independent data pre-written (notes/prereg_measurement_invariance.md). Wet-lab mechanism OPEN-EXTERNAL.
13. Zero-forskolin specificity fails - PARTIAL. Failure preserved (results/zero_fsk_control.json). CFTR-inhibitor / non-CF controls OPEN-EXTERNAL.
14. Drug-ranking reproduces published work - PARTIAL. VTI/ETI and Drevinek analyses already framed as reanalysis/validation, not discovery. "Does size correction change patient ranking?" is OPEN-EXEC exploratory.
15. GEO organoid fidelity label contamination - PARTIAL. Fibroblast-attractor retraction recorded; strict organ protocols and tool-fidelity checks exist (results/strict_*.json, tool_fidelity_*.json).
16. Metabolic-parenchyma gap weak - PARTIAL. Kidney external replication failure recorded (strict_organ_*.json); more independent datasets OPEN-EXTERNAL.
17. Statistical model limitations - PARTIAL. Donor-level exact sign flip, well-balance weighting sensitivity, donor-median aggregation exist. Hierarchical Bayesian / weighted-LMM sensitivity OPEN-EXEC (exploratory).
18. Missing endpoints / never-detected objects - PARTIAL. Bounds + tipping analyses exist (missing_endpoint_bounds.json, missing_endpoint_tipping.json); detection-confidence thresholding described. Raw images OPEN-EXTERNAL.
19. Clinical/translational overclaim risk - OPEN-EXEC (with item 1). Title/abstract language pass: state "in vitro computational reanalysis only, not a diagnostic or treatment predictor"; remove clinical-sounding framing from headings.
20. Presentation too dense and negative - OPEN-EXEC. Build the 12-slide story: one question, one method, one key result, one major limitation, one future experiment; negative arms to backup slides. (Cheapest first deliverable pair with item 1.)

## The 3 candidate reframes (verdict offers one pick)

1. Methodological CBIO: size-aware, missing-data-robust pipeline for single-organoid CFTR modulator response analysis. Verdict calls the missing-data/selection-bias model the strongest methodological angle (top-5 #4). ALIGNMENT: items 8/18 coverage already exists and becomes the contribution spine.
2. Validation CBIO: external validation of starting-size effects - blocked on OPEN-EXTERNAL independent cohort.
3. Segmentation-integrity CBIO: acquisition-field leakage inflates benchmarks - our overlap work supports it, but needs a clean independent cohort + paired baselines (OPEN-EXTERNAL for full claim).
DECISION RULE: verdict ranks them; reframe 1 is the only one executable without new external data and carries the verdict's strongest-angle note. Restructure proceeds under reframe 1 with the size question as the primary question. If the owner prefers 2 or 3, the external-data blocker is reported, not worked around.

## Top-5 additions (verdict)

1. One primary question + one primary analysis, pre-registered -> items 1/6 (this queue).
2. Independent cohort -> OPEN-EXTERNAL, report as the key blocker.
3. Fix segmentation split leakage, drop benchmark-win claim unless cleanly validated -> benchmark-win already withdrawn; group-disjoint retrain = item 10.
4. Explicit missing-data and selection-bias model -> items 8/18 consolidation into the methodological contribution.
5. Clear student ownership and AI disclosure + concise ISEF story -> items 2 (PROCESS, owner) and 20.

## Execution order (locked)

A. (items 1,5,6,19) One-question restructure plan + new abstract/section map, marked as post-verdict framing amendment. Original gates untouched.
B. (item 20) 12-slide story document.
C. (item 10) Group-disjoint segmentation split: lock a pre-registered protocol FIRST (new split from the frozen scene graph, identical code, frozen v2 training config; clearly labeled a new exploratory evaluation, not a benchmark-gate repair), then implement.
D. (item 4 of top-5) Consolidated missing-data/selection-bias methods section from existing results.
E. (items 3,9,17,14) Prior-art subsection; exploratory errors-in-variables and Bayesian sensitivity; size-corrected patient-ranking exploratory.
F. OPEN-EXTERNAL items reported as blockers: independent individual-organoid dataset, raw DIS masks/images, CFTR-inhibitor controls, independent annotation cohort, wet-lab mechanism.
G. PROCESS item 2: authorship/AI-disclosure options prepared for owner decision.
