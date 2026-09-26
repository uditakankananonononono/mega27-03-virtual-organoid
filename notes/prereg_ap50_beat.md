# PREREG VO3-R1 (draft): benchmark beat on OrgaSegment held-out AP50
Date: 2026-09-26. Status: LOCKED at commit time; any amendment gets a new dated section.
Benchmark: published OrgaSegment mean AP50 = 0.76 (Living-Technologies/
 OrgaSegment v1.0.1 scorer, author release aeaec17). Our sealed clean
 train-only checkpoint: AP50 = 0.738672 (results/trainonly_seg/sealed_eval.json,
 selection commit 120f944, fixed params [0.4,0.5,0.3,40], 512px, 3-view TTA).
Locked gate G1: a NEW model, trained on the SAME train split only, scored by
 the SAME author scorer on the SAME held-out eval split, mean AP50 > 0.76.
Protocol locks:
 - Eval split, scorer, image list, and postprocess selection protocol unchanged.
 - Model selection on the TRAIN-internal validation split only (as before);
   the sealed eval is run ONCE per candidate release, never used for tuning.
 - Permitted levers: pretrained encoder (torchvision weights), stronger
   augmentation, longer schedule, model souping of independent runs,
   input resolution within scorer constraints, architecture swap to a
   stronger open instance-segmentation baseline retrained on the same split.
 - Forbidden: eval-split training, scorer modification, cherry-picking
   per-image params outside the sealed selection protocol.
Fallback pivots if G1 stalls (rule 6, ask ChatGPT for redirection):
 P1: beat the AUTHOR model's own reproduced score on this eval split
   (ap_author column in sealed_eval.json gives per-image author AP; a
   paired per-image win vs author release is a published-comparator beat).
 P2: beat a published same-task baseline on a different locked metric
   (e.g. F1 at the author's operating point) with a fresh prereg.
