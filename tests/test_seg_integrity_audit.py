from pathlib import Path
import json,numpy as np
from seg_integrity_audit import audit

def test_real_split_and_eval_mean():
    j=json.load(open('results/seg_eval_tuned_512.json'))
    a=audit(Path('data/raw/organoid_basic'),j,n_boot=100,seed=3)
    assert a['split_counts']=={'train':184,'val':35,'eval':12}
    assert a['filename_overlap']=={'train_and_val':0,'train_and_eval':0,'val_and_eval':0}
    assert a['heldout_eval_mAP50']<a['published_mAP50']
