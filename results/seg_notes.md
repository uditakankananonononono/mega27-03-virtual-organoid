# Segmentation notes (OrgaSegment eval split, 12 images; published mAP@0.5 0.76 +/- 0.12)
- v1 U-Net (256 px, default post-processing): mAP50 0.634 +/- 0.157 (results/seg_eval.json).
- v1 + post-processing grid tuned on val + flip TTA: eval mAP50 0.710 +/- 0.127 (results/seg_eval_tuned_256.json). Params fg 0.5, seed 0.5, boundary 0.5, min size 10.
  Caveat: val images were also in the training set, so the val score used for tuning (0.792) is optimistic. The eval split is untouched and scored once.
  Verdict: below the published SOTA (0.76), within one SD. Not a beat.
- v2 (512 px) first run was OOM-killed at epoch 18 because it ran alongside tuning (log_seg512_oomkilled_ep18.txt). It was restarted with 3-epoch checkpoints.

## Split-integrity correction (registered 543bd2d, results/seg_integrity_audit.json)
Network training in `train_seg.py` includes 184 train + 35 val images; `tune_seg.py` chooses post-processing on the first 20 val images, so the quoted validation mAP@0.5 of 0.780 is **not independent** of network fitting and must not be treated as an unbiased validation result. Train/val/eval filenames are pairwise disjoint, and the 12-image eval split was not used in the network-training or post-processing grid. The held-out eval mAP@0.5 is 0.733, below the published OrgaSegment 0.76. Bootstrap over our 12 eval images gives 95% CI [0.663,0.808], difference to the published point mean [-0.097,0.048]. This is not a paired leader comparison: published per-image predictions were not available, and exact metric equivalence is unverified. No benchmark break. Future training must reserve val strictly for model/tuning selection and eval only once; we did not retrain or tune against eval in this audit.
