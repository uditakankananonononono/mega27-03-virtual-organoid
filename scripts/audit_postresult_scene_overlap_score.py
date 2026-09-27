"""Retrospective segmentation score sensitivity; never a sealed benchmark gate."""
import hashlib, json
from pathlib import Path
import numpy as np
s = Path('results/scene_graph_forensics.json')
v = Path('results/trainonly_seg/sealed_eval_v2.json')
b = Path('results/trainonly_seg/sealed_eval.json')
g = json.loads(s.read_text())
a = json.loads(v.read_text()); old = json.loads(b.read_text())
score = {z['image']: z['ap_author'] for z in a['per_image']}
prior = {z['image']: z['ap_author'] for z in old['per_image']}
assert set(score) == set(prior) and len(score) == 12
assert abs(np.mean(list(score.values()))-a['mean_author']) < 1e-12
rows = []
for h in g['threshold_sensitivity']:
    flagged = set()
    for comp in h['multi_frame_components']:
        members = comp['members']
        if not any(z.startswith('eval/') for z in members): continue
        if any(z.startswith(('train/','val/')) for z in members):
            flagged |= {z.split('/',1)[1] for z in members if z.startswith('eval/')}
    assert flagged <= score.keys()
    kept = sorted(score.keys()-flagged)
    assert kept
    rng=np.random.default_rng(20260927)
    vals=np.array([score[k] for k in kept])
    draws=vals[rng.integers(0,len(kept),size=(20000,len(kept)))].mean(axis=1)
    rows.append({'thumbnail_cosine_threshold':h['threshold'],
       'flagged_eval':[{'image':k,'v2_ap50':score[k], 'v1_ap50':prior[k]} for k in sorted(flagged)],
       'n_kept':len(kept),'v2_kept_mean':float(vals.mean()),
       'v2_flagged_mean':float(np.mean([score[k] for k in flagged])) if flagged else None,
       'v1_kept_mean':float(np.mean([prior[k] for k in kept])),
       'v2_minus_v1_kept':float(np.mean([score[k]-prior[k] for k in kept])),
       'kept_bootstrap_ci95':list(map(float,np.quantile(draws,[.025,.975]))),
       'retained_image_scores':[{'image':k,'v2_ap50':score[k]} for k in kept]})
out={'status':'retrospective threshold-sensitivity audit; no independent benchmark gate',
     'protocol':'notes/postresult_scene_overlap_score_audit.md',
     'source_urls':['https://zenodo.org/records/10278229','https://www.nature.com/articles/s42003-024-05966-4'],
     'inputs_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (s,v,b)},
     'full_eval_n':len(score),'full_eval_v2_mean':float(a['mean_author']),
     'full_eval_v1_mean':float(old['mean_author']), 'published_rounded_mean':.76,
     'thresholds':rows,
     'limits':'All exclusions and scores are post-result. Retained images may share unflagged fields or donors. Published author-model paired predictions and independent acquisition metadata unavailable. Do not infer a clean held-out win or tune the model.'}
Path('results/postresult_scene_overlap_score_audit.json').write_text(json.dumps(out,indent=2)+'\n')
for z in rows:
 print(z['thumbnail_cosine_threshold'], 'flagged', [(x['image'],round(x['v2_ap50'],5)) for x in z['flagged_eval']], 'kept',z['n_kept'], 'v2 mean',round(z['v2_kept_mean'],6), 'CI',z['kept_bootstrap_ci95'],'v1',round(z['v1_kept_mean'],6))
