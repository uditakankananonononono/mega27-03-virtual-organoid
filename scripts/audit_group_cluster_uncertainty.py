"""Post-result 0.95 candidate-scene-group bootstrap; no evaluation or retuning."""
import json
from pathlib import Path
import numpy as np
G=json.loads(Path('results/scene_graph_forensics.json').read_text())
E=json.loads(Path('results/group_disjoint_seg/exploratory_eval.json').read_text())
assert E['status'].startswith('EXPLORATORY')
assert G['image_manifest_sha256']==E['source_manifest_sha256']
score={x['image']:x['ap_author'] for x in E['per_image']}
assert len(score)==14
components=[];seen=set()
for c in G['threshold_sensitivity'][0]['multi_frame_components']:
    members=sorted(set(Path(z).name for z in c['members'])&score.keys())
    if members:components.append(members);seen.update(members)
components.extend([[name] for name in sorted(score.keys()-seen)])
assert len(components)==12 and sorted(map(len,components))==[1]*10+[2]*2
assert sorted(m for c in components for m in c)==sorted(score)
components.sort(key=lambda c:c[0])
group=np.array([np.mean([score[m] for m in c]) for c in components])
rng=np.random.default_rng(20260927)
means=group[rng.integers(0,12,size=(20000,12))].mean(axis=1)
res={'status':'POST-RESULT DESCRIPTIVE: proxy-group resampling does not establish donor independence',
    'source':'results/group_disjoint_seg/exploratory_eval.json',
    'graph_source':'results/scene_graph_forensics.json',
    'threshold':0.95,'groups':[{ 'members':c,'n':len(c),'mean_author_ap50':float(np.mean([score[m] for m in c]))} for c in components],
    'n_groups':12,'n_images':14,'n_singleton_groups':10,'n_two_image_groups':2,
    'file_weighted_ap50':float(np.mean(list(score.values()))),
    'file_weighted_image_bootstrap_ci95':E['image_bootstrap_ci95'],
    'equal_group_weight_ap50':float(group.mean()),
    'equal_group_bootstrap_ci95':list(map(float,np.quantile(means,[.025,.975]))),
    'draws':20000,'seed':20260927,
    'limitations':'Proxy image groups are selected after original AP result, not known donor/plate units; 12 groups and previous test-set exposure; group and file means use different weights; no author-model paired predictions; no benchmark gate upgrade.'}
Path('results/group_disjoint_seg/postresult_group_cluster_uncertainty.json').write_text(json.dumps(res,indent=2)+'\n')
print({k:v for k,v in res.items() if k!='groups'})
