"""Post-evaluation descriptive uncertainty on saved sealed predictions, no eval retuning."""
import json,hashlib
import numpy as np
from scipy.stats import binomtest
v=json.load(open('results/trainonly_seg/sealed_eval_v2.json'))
b=json.load(open('results/trainonly_seg/sealed_eval.json'))
assert v['author_release_tag']=='aeaec17c42fefde61a8b2215a616d419e8126e9f'
old={x['image']:x['ap_author'] for x in b['per_image']};new={x['image']:x['ap_author'] for x in v['per_image']}
assert set(old)==set(new) and len(new)==12
ids=sorted(new);a=np.array([new[i] for i in ids]);o=np.array([old[i] for i in ids]);d=a-o
assert abs(a.mean()-v['mean_author'])<1e-12 and abs(o.mean()-b['mean_author'])<1e-12
rng=np.random.default_rng(20260926);ind=rng.integers(0,12,size=(20000,12));am=a[ind].mean(1);dm=d[ind].mean(1)
nonzero=np.abs(d)>1e-7
out={'status':'POST-EVAL DESCRIPTIVE; cannot upgrade frozen numerical gate','author_scorer_v2_mean':float(a.mean()),'published_rounded_mean':.76,
 'v2_image_sd_sample':float(a.std(ddof=1)),'n_images':12,'bootstrap_seed':20260926,'bootstrap_draws':20000,
 'v2_image_bootstrap_ci95':[float(x) for x in np.quantile(am,[.025,.975])],
 'v2_minus_clean_v1_paired_mean':float(d.mean()),'v2_minus_v1_image_bootstrap_ci95':[float(x) for x in np.quantile(dm,[.025,.975])],
 'paired_improved_images':int((d>1e-7).sum()),'paired_declined_images':int((d< -1e-7).sum()),'paired_tied_images':int((~nonzero).sum()),
 'paired_nonzero_sign_p_two_sided':float(binomtest(int((d>1e-7).sum()),int(nonzero.sum()),.5).pvalue),
 'leave_one_image_out':[{'image':ids[i],'v2_mean':float(np.delete(a,i).mean()),'delta_vs_rounded_0.76':float(np.delete(a,i).mean()-.76)} for i in range(12)],
 'input_sha256':{p:hashlib.sha256(open(p,'rb').read()).hexdigest() for p in ['results/trainonly_seg/sealed_eval_v2.json','results/trainonly_seg/sealed_eval.json']},
 'limit':'image resampling is conditional on 12 published eval images; author-model head-to-head predictions unavailable'}
open('results/ap50_uncertainty_audit.json','w').write(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='leave_one_image_out'},indent=2))
print('LOO below benchmark',sum(z['v2_mean']<=.76 for z in out['leave_one_image_out']))
