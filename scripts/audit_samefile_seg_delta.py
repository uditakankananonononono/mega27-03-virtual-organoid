"""Post-result same-file reconciliation of original-v2 and reassigned exploratory eval."""
import json,numpy as np
from pathlib import Path
old=json.load(open('results/trainonly_seg/sealed_eval_v2.json'))
new=json.load(open('results/group_disjoint_seg/exploratory_eval.json'))
split=json.load(open('results/group_disjoint_seg_split.json'))
a={r['image']:r for r in old['per_image']};b={r['image']:r for r in new['per_image']}
assert len(a)==12 and len(b)==14 and set(a)<=set(b) and len(set(b)-set(a))==2
assert split['split_counts']['eval']==14 and split['source_manifest_sha256']==new['source_manifest_sha256']
groups={}
for z in split['moved']:
 if z['to']=='eval':
  for member in z['component']:groups[member.split('/',1)[1]]=z['component']
rows=[]
for image in sorted(a):
 x,y=a[image],b[image]
 assert x['n_gt']==y['n_gt']
 rows.append({'image':image,'old_ap50':x['ap_author'],'new_ap50':y['ap_author'],
    'new_minus_old_ap50':y['ap_author']-x['ap_author'],'n_gt':x['n_gt'],
    'old_tp_fp_fn':[x['tp_author'],x['fp_author'],x['fn_author']],
    'new_tp_fp_fn':[y['tp_author'],y['fp_author'],y['fn_author']],
    'reassigned_eval_proxy_component':groups.get(image,[f'eval/{image}'])})
d=np.array([r['new_minus_old_ap50'] for r in rows])
assert abs(float(np.mean(d))-(-.020398070414861042))<1e-12
moved=sorted(set(b)-set(a))
assert moved==sorted(z['image'].split('/',1)[1] for z in split['moved'] if z['to']=='eval')
result={'status':'post-result same-file descriptive reconciliation; no causal leakage estimate or independent benchmark',
 'protocol':'notes/postresult_samefile_seg_delta_protocol.md',
 'old_result':'results/trainonly_seg/sealed_eval_v2.json',
 'new_result':'results/group_disjoint_seg/exploratory_eval.json',
 'source_manifest_sha256':split['source_manifest_sha256'],
 'n_shared_original_eval':len(rows),'n_new_moved_eval':len(moved),'moved_from_train_val':moved,
 'shared_mean_old':float(np.mean([r['old_ap50'] for r in rows])),
 'shared_mean_new':float(np.mean([r['new_ap50'] for r in rows])),
 'mean_new_minus_old':float(np.mean(d)),'median_new_minus_old':float(np.median(d)),
 'min_new_minus_old':float(np.min(d)),'max_new_minus_old':float(np.max(d)),
 'improve':int(sum(d>0)),'worsen':int(sum(d<0)),'ties':int(sum(d==0)),
 'old_tp_fp_fn_sums':[sum(r['old_tp_fp_fn'][i] for r in rows) for i in range(3)],
 'new_tp_fp_fn_sums':[sum(r['new_tp_fp_fn'][i] for r in rows) for i in range(3)],
 'shared_images':rows,
 'limits':'Outcomes already observed; changed weights, split, selection and postprocessing. Shared files are not independent acquisitions. Per-image scores are correlated and ground truth masks reused.'}
Path('results/group_disjoint_seg/postresult_samefile_delta.json').write_text(json.dumps(result,indent=2)+'\n')
print({k:v for k,v in result.items() if k in ('n_shared_original_eval','n_new_moved_eval','mean_new_minus_old','median_new_minus_old','min_new_minus_old','max_new_minus_old','improve','worsen','ties','old_tp_fp_fn_sums','new_tp_fp_fn_sums')})
