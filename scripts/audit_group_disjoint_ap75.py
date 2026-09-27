"""Frozen post-result strict-IoU diagnostic, never a benchmark gate."""
import hashlib,importlib.util,json,subprocess,sys
from pathlib import Path
import numpy as np,torch,torch.nn.functional as F
sys.path.insert(0,'src')
from vorganoid.seg import UNet,load_pair
from eval_group_disjoint_exploratory import instances,PARAMS,SHA,MANIFEST_SHA,AUTHOR_SHA,SIZE
root=Path('results/group_disjoint_seg');author=Path('/tmp/orgasegment_author')
assert subprocess.check_output(['git','-C',str(author),'rev-parse','HEAD'],text=True).strip()==AUTHOR_SHA
spec=importlib.util.spec_from_file_location('author_metrics',author/'lib'/'metrics.py')
metrics=importlib.util.module_from_spec(spec);spec.loader.exec_module(metrics)
old=json.load(open(root/'exploratory_eval.json'));split=json.load(open('results/group_disjoint_seg_split.json'))
assert old['model_sha256']==SHA and hashlib.sha256((root/'best.pt').read_bytes()).hexdigest()==SHA
assert old['source_manifest_sha256']==split['source_manifest_sha256']==MANIFEST_SHA
assert tuple(old['params'])==PARAMS and len(old['per_image'])==14
ck=torch.load(root/'best.pt',map_location='cpu',weights_only=False);net=UNet();net.load_state_dict(ck['net']);net.eval();torch.set_num_threads(1)
rows=[]
for prev in old['per_image']:
 img=Path('data/derived/group_disjoint_seg/eval')/prev['image'];msk=img.with_name(img.name.replace('_img.jpg','_masks_organoid.png'))
 assert img.is_file() and msk.is_file()
 x,_,lab=load_pair(str(img),str(msk),SIZE);xt=torch.tensor(x)[None,None]
 with torch.no_grad():
  p=(torch.softmax(net(xt),1)+torch.softmax(net(xt.flip(-1)),1).flip(-1)+torch.softmax(net(xt.flip(-2)),1).flip(-2))/3
 pu=F.interpolate(p,size=lab.shape,mode='bilinear')[0].numpy()
 pred=instances(pu,*PARAMS[:3],int(PARAMS[3]*(lab.shape[0]/SIZE)**2));ids=np.unique(pred);pred=np.searchsorted(ids,pred).astype(np.int32)
 ap50,tp50,fp50,fn50=metrics.average_precision(lab.astype(np.int32),pred,.5)
 assert abs(float(ap50[0])-prev['ap_author'])<1e-7 and [int(tp50[0]),int(fp50[0]),int(fn50[0])]==[prev[z] for z in ('tp_author','fp_author','fn_author')],prev['image']
 ap75,tp75,fp75,fn75=metrics.average_precision(lab.astype(np.int32),pred,.75)
 rows.append({'image':prev['image'],'ap50_reproduced':float(ap50[0]),'ap75_author':float(ap75[0]),
   'ap75_minus_ap50':float(ap75[0]-ap50[0]),'n_gt':prev['n_gt'],
   'tp_fp_fn_50':[int(tp50[0]),int(fp50[0]),int(fn50[0])],
   'tp_fp_fn_75':[int(tp75[0]),int(fp75[0]),int(fn75[0])]})
a=np.array([r['ap75_author'] for r in rows]);d=np.array([r['ap75_minus_ap50'] for r in rows]);assert np.all(d<=1e-6)
out={'status':'post-result strict-IoU diagnostic on known source files; not independent benchmark',
 'protocol':'notes/postresult_ap75_boundary_protocol.md','checkpoint_sha256':SHA,
 'source_manifest_sha256':MANIFEST_SHA,'author_release_commit':AUTHOR_SHA,
 'n_images':len(rows),'ap50_reproduced_mean':float(np.mean([r['ap50_reproduced'] for r in rows])),
 'ap75_mean':float(a.mean()),'mean_ap75_minus_ap50':float(d.mean()),
 'median_ap75_minus_ap50':float(np.median(d)),
 'tp_fp_fn_50_sums':[sum(r['tp_fp_fn_50'][i] for r in rows) for i in range(3)],
 'tp_fp_fn_75_sums':[sum(r['tp_fp_fn_75'][i] for r in rows) for i in range(3)],
 'per_image':rows,'limits':'Post-result diagnostic; AP75 has no paired released-model comparison; 14 known source files and proxy acquisition groups; strict-IoU gap cannot isolate model localization from annotation and acquisition effects.'}
(root/'postresult_ap75_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k in ('n_images','ap50_reproduced_mean','ap75_mean','mean_ap75_minus_ap50','median_ap75_minus_ap50','tp_fp_fn_50_sums','tp_fp_fn_75_sums')})
