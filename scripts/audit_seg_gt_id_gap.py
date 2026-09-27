"""Explain scorer mean discrepancy using saved per-image labels/counts only."""
import json
from pathlib import Path
import numpy as np
p=Path('results/group_disjoint_seg/exploratory_eval.json')
s=json.loads(p.read_text());r=s['per_image']
assert len(r)==14
rows=[]
for x in r:
 d=x['ap_ours']-x['ap_author']
 rows.append({'image':x['image'],'gt_id_gaps':x['gt_id_gaps'],
              'author_ap50':x['ap_author'],'unique_count_ap50':x['ap_ours'],
              'difference_unique_minus_author':d,
              'author_tp_fp_fn':[x['tp_author'],x['fp_author'],x['fn_author']]})
gap=[x for x in rows if x['gt_id_gaps']>0]
assert len(gap)==2 and all(x['gt_id_gaps']==1 for x in gap)
assert all(abs(x['difference_unique_minus_author'])<1e-6 for x in rows if x['gt_id_gaps']==0)
assert abs(np.mean([x['author_ap50'] for x in rows])-s['mean_author'])<1e-12
assert abs(np.mean([x['unique_count_ap50'] for x in rows])-s['mean_ours'])<1e-12
out={'status':'POST-RESULT metric-accounting audit, no mask or model change',
     'n_images':14,'n_gt_id_gap_images':len(gap),'gap_images':gap,
     'mean_author':s['mean_author'],'mean_unique_count':s['mean_ours'],
     'mean_difference_unique_minus_author':s['mean_ours']-s['mean_author'],
     'n_other_images_with_difference_gt_1e_6':sum(abs(x['difference_unique_minus_author'])>1e-6 for x in rows if not x['gt_id_gaps']),
     'limits':'Author scorer used for primary exploratory AP50; no independent author model predictions or benchmark repair.'}
Path('results/group_disjoint_seg/postresult_gt_id_gap_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k!='gap_images'})
for x in gap:print(x['image'],x['author_tp_fp_fn'],x['difference_unique_minus_author'])
