"""Frozen U-Net mouse-intestinal YOLO class-agnostic detection transfer, no tuning."""
import hashlib,json,zipfile
from pathlib import Path
import numpy as np,torch,torch.nn.functional as F
from PIL import Image
from scipy import ndimage as ndi
from scipy.optimize import linear_sum_assignment
from skimage.segmentation import watershed
root=Path('data/external/mouse_intestinal_yolo');archive=root/'OrganoidDataset.zip'
assert hashlib.md5(archive.read_bytes()).hexdigest()=='98f10253594f37fdee809dd61d618861'
data=root/'OrganoidDataset';files=sorted((data/'val/images').glob('*.jpeg'))
assert len(files)==84 and len(list((data/'train/images').glob('*.jpeg')))==756
assert all((data/'val/labels'/(f.stem+'.txt')).exists() for f in files)
ckfile=Path('results/trainonly_seg/best_v2.pt');cksha=hashlib.sha256(ckfile.read_bytes()).hexdigest()
selectfile=Path('results/trainonly_seg/val_selection_v2.json');sel=json.loads(selectfile.read_text());assert sel['selected_params']==[.5,.5,.3,40] and sel['eval_used'] is False
import sys
sys.path.insert(0,'src')
from vorganoid.seg import UNet
ck=torch.load(ckfile,map_location='cpu',weights_only=True);assert ck['epoch']==29
net=UNet();net.load_state_dict(ck['net']);net.eval();torch.set_num_threads(1)
def bounds(mask):
 ys,xs=np.nonzero(mask);return [int(xs.min()),int(ys.min()),int(xs.max()+1),int(ys.max()+1)]
def parse_yolo(file,W,H):
 out=[]; classes=[]
 for line in file.read_text().splitlines():
  a=line.split();assert len(a)==5
  c,x,y,w,h=map(float,a);assert c in (0,1,2,3) and 0<=x<=1 and 0<=y<=1 and 0<w<=1 and 0<h<=1
  out.append([max(0.,(x-w/2)*W),max(0.,(y-h/2)*H),min(W,(x+w/2)*W),min(H,(y+h/2)*H)]);classes.append(int(c))
 return np.array(out,dtype=float).reshape(-1,4),classes
def pair_score(gt,pred):
 ng,npred=len(gt),len(pred)
 if not ng or not npred:return 0,npred,ng
 ix1=np.maximum(gt[:,None,0],pred[None,:,0]);iy1=np.maximum(gt[:,None,1],pred[None,:,1]);ix2=np.minimum(gt[:,None,2],pred[None,:,2]);iy2=np.minimum(gt[:,None,3],pred[None,:,3])
 inter=np.maximum(0,ix2-ix1)*np.maximum(0,iy2-iy1)
 ag=(gt[:,2]-gt[:,0])*(gt[:,3]-gt[:,1]);ap=(pred[:,2]-pred[:,0])*(pred[:,3]-pred[:,1]);iou=inter/(ag[:,None]+ap[None,:]-inter+1e-12)
 g,p=linear_sum_assignment(-iou);tp=int((iou[g,p]>=.5).sum());return tp,npred-tp,ng-tp
allrows=[]; outpath=Path('results/postresult_mouse_yolo_domain.json')
if outpath.exists():
 prior=json.loads(outpath.read_text());assert prior['archive_md5']=='98f10253594f37fdee809dd61d618861' and prior['checkpoint_sha256']==cksha
 assert [r['image'] for r in prior['per_image']]==[f.name for f in files[:len(prior['per_image'])]]
 allrows=prior['per_image']
for file in files[len(allrows):]:
 source=Image.open(file).convert('L');W,H=source.size
 x=np.asarray(source.resize((512,512),Image.Resampling.BILINEAR),dtype=np.float32)/255.
 x=(x-x.mean())/(x.std()+1e-6);xt=torch.from_numpy(x)[None,None]
 with torch.no_grad():
  prob=(torch.softmax(net(xt),1)+torch.softmax(net(xt.flip(-1)),1).flip(-1)+torch.softmax(net(xt.flip(-2)),1).flip(-2))/3
 pu=F.interpolate(prob,size=(H,W),mode='bilinear')[0].numpy()
 fg=(pu[1]+pu[2])>.5;seeds=(pu[1]>.5)&(pu[2]<.3)
 markers,_=ndi.label(seeds);lab=watershed(-pu[1],markers,mask=fg)
 sizes=np.bincount(lab.ravel());minimum=int(40*W*H/(512*512))
 boxes=np.asarray([bounds(lab==i) for i in range(1,len(sizes)) if sizes[i]>=minimum],dtype=float).reshape(-1,4)
 gt,classes=parse_yolo(data/'val/labels'/(file.stem+'.txt'),W,H)
 tp,fp,fn=pair_score(gt,boxes)
 row={'image':file.name,'image_sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'annotation_sha256':hashlib.sha256((data/'val/labels'/(file.stem+'.txt')).read_bytes()).hexdigest(),
      'size':[W,H],'n_truth':len(gt),'n_pred':len(boxes),'classes':{str(i):classes.count(i) for i in range(4)},
      'tp':tp,'fp':fp,'fn':fn,'precision':tp/(tp+fp) if tp+fp else 0.,'recall':tp/(tp+fn) if tp+fn else 0.,
      'f1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.}
 allrows.append(row)
 result={'status':'post-result external-domain class-agnostic bounding-box detection check, not a mask benchmark win',
  'source_url':'https://zenodo.org/records/6768583','protocol':'notes/postresult_mouse_yolo_domain_protocol.md',
  'archive_md5':'98f10253594f37fdee809dd61d618861','checkpoint_sha256':cksha,
  'selection_sha256':hashlib.sha256(selectfile.read_bytes()).hexdigest(),
  'n_expected_val':84,'per_image':allrows,
  'limits':'Mouse 4x EVOS detection boxes; model trained on different human images/masks. Annotation classes collapsed; no comparable published baseline predictions, acquisition-group independence unverified, not human mask AP50 or biological validation.'}
 outpath.write_text(json.dumps(result,indent=2)+'\n')
 print(len(allrows),file.name,'TP',tp,'FP',fp,'FN',fn,flush=True)
if len(allrows)==84:
 t=sum(r['tp'] for r in allrows);fp=sum(r['fp'] for r in allrows);fn=sum(r['fn'] for r in allrows)
 result['summary']={'tp':t,'fp':fp,'fn':fn,'micro_precision':t/(t+fp),'micro_recall':t/(t+fn),'micro_f1':2*t/(2*t+fp+fn),
    'macro_image_f1':float(np.mean([r['f1'] for r in allrows])),
    'zero_ground_truth_images':sum(r['n_truth']==0 for r in allrows),
    'zero_prediction_images':sum(r['n_pred']==0 for r in allrows)}
 outpath.write_text(json.dumps(result,indent=2)+'\n');print('SUMMARY',result['summary'],flush=True)
