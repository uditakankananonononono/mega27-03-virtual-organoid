# Segmentation notes (OrgaSegment eval split, 12 images; published mAP@0.5 0.76 +/- 0.12)
- v1 U-Net (256 px, default post-processing): mAP50 0.634 +/- 0.157 (results/seg_eval.json).
- v1 + post-processing grid tuned on val + flip TTA: eval mAP50 0.710 +/- 0.127 (results/seg_eval_tuned_256.json). Params fg 0.5, seed 0.5, boundary 0.5, min size 10.
  Caveat: val images were also in the training set, so the val score used for tuning (0.792) is optimistic. The eval split is untouched and scored once.
  Verdict: below the published SOTA (0.76), within one SD. Not a beat.
- v2 (512 px) first run was OOM-killed at epoch 18 because it ran alongside tuning (log_seg512_oomkilled_ep18.txt). It was restarted with 3-epoch checkpoints.
