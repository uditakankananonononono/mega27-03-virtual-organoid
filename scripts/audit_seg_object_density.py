"""Post-result descriptive AP50/GT-object-density diagnostic on saved rows."""
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
E=json.loads(Path('results/group_disjoint_seg/exploratory_eval.json').read_text())
S=json.loads(Path('results/group_disjoint_seg_split.json').read_text())
rows=E['per_image'];assert len(rows)==14
moved={x['image'].split('/')[-1] for x in S['moved'] if x['to']=='eval'}
assert len(moved)==2
med=float(np.median([x['n_gt'] for x in rows]))

def summary(z):
    tp=sum(x['tp_author'] for x in z);fp=sum(x['fp_author'] for x in z);fn=sum(x['fn_author'] for x in z)
    assert all(x['tp_author']+x['fn_author'] in (x['n_gt'],x['n_gt']+x['gt_id_gaps']) for x in z)
    return {'n_images':len(z),'n_moved_to_eval':sum(x['image'] in moved for x in z),
            'mean_author_ap50':float(np.mean([x['ap_author'] for x in z])),
            'median_author_ap50':float(np.median([x['ap_author'] for x in z])),
            'sum_tp':tp,'sum_fp':fp,'sum_fn':fn,
            'aggregate_fn_fraction':fn/(tp+fn),
            'per_image_fn_fraction_mean':float(np.mean([x['fn_author']/(x['tp_author']+x['fn_author']) for x in z])),
            'mean_fp_per_image':fp/len(z),
            'images':[{'image':x['image'],'n_gt':x['n_gt'],'author_ap50':x['ap_author'],
                       'tp':x['tp_author'],'fp':x['fp_author'],'fn':x['fn_author'],
                       'moved_to_eval':x['image'] in moved} for x in z]}
lo=[x for x in rows if x['n_gt']<=med];hi=[x for x in rows if x['n_gt']>med]
assert len(lo)+len(hi)==14 and lo and hi
rho=float(spearmanr([x['n_gt'] for x in rows],[x['ap_author'] for x in rows]).statistic)
out={'status':'POST-RESULT descriptive density stratification; not a causal effect or benchmark test',
     'median_gt_cutoff':med,'low_includes_median':True,'low':summary(lo),'high':summary(hi),
     'spearman_gt_count_vs_author_ap50':rho,
     'limits':'14 previously seen eval images; GT counts are annotation density, not validated acquisition conditions; detector and source fields vary; no significance claim, model retuning or benchmark gate upgrade.'}
Path('results/group_disjoint_seg/postresult_object_density.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v if k not in ('low','high') else {q:w for q,w in v.items() if q!='images'} for k,v in out.items()})
