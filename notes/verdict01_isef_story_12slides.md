# Verdict-01 12-slide ISEF story (queue item B) - 27 September 2026

One question, one method, one key result, one major limitation, one future experiment. Negative arms live in backup. Slide text, not design.

1. Title. "Does starting organoid size modify CFTR-modulator response? A size-aware, missing-data-robust reanalysis of single-organoid swelling." Name, program, one-line disclosure footer (AI-assistance disclosure per fair rules - exact wording pending owner decision, queue item 2).
2. The problem. CFTR-modulator response is measured by forskolin-induced swelling of patient intestinal organoids. Organoids in one well differ in size by >10x. If response depends on starting size, well-level readouts depend on well size mix.
3. The question. Does starting size modify the single-organoid response after accounting for track attrition? One sentence. Everything else is backup.
4. The data. Public single-organoid dataset (Zenodo 10610438): 17 CF donors, ~14,800 positive-forskolin organoid observations with baseline and endpoint areas, DMSO and three modulator arms. Strengths and what is missing (no raw masks, no never-detected objects).
5. The integrity problem nobody sees. Tracked organoids are not a random subset: endpoint observability varies with starting size, treatment arm and local crowding (stratified neighbor analysis; donor-clustered observability model). Attrition can create or hide a size effect.
6. The method. A size-aware, attrition-audited pipeline: donor-matched size contrasts; pre-registered cutoffs and gates stated before scoring; explicit observability modeling; donor-level exact tests; sensitivity bounds for missing endpoints.
7. The primary pre-registered tests, honestly. Sign-only donor test: 9/12, p=0.073 - FAIL. Zero-forskolin specificity: 3/11, p=0.967 - FAIL. Matched plate-dose stress test: FAIL. The naive story does not survive its own gates.
8. What survives. Under the locked donor-block analysis (31 blocks, 12 donors), the magnitude-aware treatment-by-large coefficient is +0.253 with exact donor-level sign-flip p=0.00293; exploratory octile contrast 0.336 -> 0.55-0.64 rising with size; top-minus-bottom 0.251 [0.184, 0.323]. Within-accession, not independent replication.
9. The methodological contribution. Quantified, reproducible evidence that track attrition is size- and arm-structured in this assay class - the thing a well-level analysis cannot see. This is the result other CFTR swelling studies can use.
10. The major limitation. One accession, no independent cohort, no raw masks; cutoffs data-derived; mechanism untested. Explicitly: in vitro computational reanalysis only, not a diagnostic or treatment predictor.
11. The future experiment. Pre-registered design: independent individual-organoid cohort with donor metadata and raw masks; lock cutoffs externally; CFTR-inhibitor and zero-forskolin controls; acquisition-group-disjoint segmentation evaluation.
12. Summary. One question, one method, one honest answer: a size interaction is measurable but not yet independently established; the pipeline that proves what attrition can do is the contribution.

Backup slides (not in the 12): segmentation split-leakage forensics and withdrawn benchmark; Drevinek/VTI-ETI drug-ranking replications; GEO organ-fidelity arms; all negative arms and gate ledger.
