"""Audit split integrity, val leakage and held-out eval uncertainty.

Pre-registered results/preregistration_seg_audit.md. Never uses eval labels to tune.
"""
import json
from pathlib import Path
import numpy as np


def audit(root,record,n_boot=10000,seed=0):
    sets={k:{p.name for p in (root/k).glob('*_img.jpg')} for k in ('train','val','eval')}
    overlap={f'{a}_and_{b}':len(sets[a]&sets[b]) for a,b in [('train','val'),('train','eval'),('val','eval')]}
    a=np.asarray(record['per_image'],float)
    if len(a)!=len(sets['eval']):raise AssertionError('score count != eval images')
    if not np.isclose(a.mean(),record['eval_mAP50']):raise AssertionError('stored eval mean disagrees')
    b=a[np.random.default_rng(seed).integers(0,len(a),(n_boot,len(a)))].mean(axis=1)
    return {'split_counts':{k:len(v) for k,v in sets.items()},'filename_overlap':overlap,
            'train_script_uses':'train+val (see train_seg.py line 10)',
            'tune_script_uses':'val first 20; post-processing selected on model-training images (see tune_seg.py)',
            'heldout_eval_mAP50':float(a.mean()),'published_mAP50':record['published_mAP50'],
            'difference_to_published':float(a.mean()-record['published_mAP50']),
            'eval_image_bootstrap_ci95':list(map(float,np.quantile(b,[.025,.975]))),
            'eval_image_bootstrap_difference_ci95':list(map(float,np.quantile(b-record['published_mAP50'],[.025,.975]))),
            'eval_image_bootstrap_frac_above_point_baseline':float(np.mean(b>record['published_mAP50'])),
            'verdict':'NO BENCHMARK BREAK; val score is biased by overlap; eval score below published mean',
            'limitations':'Neither matched baseline image predictions nor guaranteed identical AP implementation; bootstrap uses only our 12 images and fixed published point mean.'}

if __name__=='__main__':
    j=json.load(open('results/seg_eval_tuned_512.json'))
    out=audit(Path('data/raw/organoid_basic'),j)
    out['preregistration']='results/preregistration_seg_audit.md (commit 543bd2d)'
    Path('results/seg_integrity_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
