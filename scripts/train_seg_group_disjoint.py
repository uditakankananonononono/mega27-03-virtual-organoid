"""Fresh U-Net on candidate scene group split; resumable per frozen protocol."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np,torch,torch.nn.functional as F
sys.path.insert(0,'src')
from vorganoid.seg import UNet,load_pair
p=argparse.ArgumentParser();p.add_argument('--steps',type=int,default=1);a=p.parse_args();assert 1<=a.steps<=3
S=512;CROP=256;MAX=60
root=Path('results/group_disjoint_seg');root.mkdir(exist_ok=True)
CK=root/'resume.pt';BEST=root/'best.pt';LOG=root/'history.json'
split=Path('data/derived/group_disjoint_seg')
manifest=json.load(open('results/group_disjoint_seg_split.json'))
assert manifest['source_manifest_sha256']=='ba003cd9645067eb455f5a2d0166a53a0a41aa73378f938e848f2ca9c73a9b65'
assert manifest['split_counts']=={'train':180,'val':37,'eval':14}
def pairs(name):
 imgs=sorted((split/name).glob('*_img.jpg'))
 assert len(imgs)==manifest['split_counts'][name]
 return [(i,i.with_name(i.name.replace('_img.jpg','_masks_organoid.png'))) for i in imgs]
train,val,ev=pairs('train'),pairs('val'),pairs('eval')
assert set(i.name for i,m in train).isdisjoint(i.name for i,m in val+ev)
assert set(i.name for i,m in val).isdisjoint(i.name for i,m in ev)
g=json.load(open('results/scene_graph_forensics.json'))
assert g['image_manifest_sha256']==manifest['source_manifest_sha256']
assign={i.name:s for s in ('train','val','eval') for i,m in pairs(s)}
for c in g['threshold_sensitivity'][0]['multi_frame_components']:
 assert len({assign[z.split('/')[-1]] for z in c['members']})==1
tr=[load_pair(str(i),str(m),S)[:2] for i,m in train];va=[load_pair(str(i),str(m),S)[:2] for i,m in val]
X=torch.tensor(np.stack([z[0] for z in tr]))[:,None];Y=torch.tensor(np.stack([z[1].astype(np.uint8) for z in tr]));del tr
VX=torch.tensor(np.stack([z[0] for z in va]))[:,None];VY=torch.tensor(np.stack([z[1].astype(np.uint8) for z in va]));del va
torch.set_num_threads(1)
net=UNet();opt=torch.optim.Adam(net.parameters(),2e-3);w=torch.tensor([1.,1.,3.]);best=float('inf');hist=[]
if CK.exists():
 ck=torch.load(CK,map_location='cpu',weights_only=False)
 assert ck['source_manifest_sha256']==manifest['source_manifest_sha256']
 net.load_state_dict(ck['net']);opt.load_state_dict(ck['opt']);torch.set_rng_state(ck['torch_rng']);np.random.set_state(ck['numpy_rng']);start=ck['epoch']+1;best=ck['best'];hist=ck['history'];print('RESUME',start,flush=True)
else:
 assert not BEST.exists() and not LOG.exists(), 'no silent reset of previous run'
 torch.manual_seed(2026);np.random.seed(2026);net=UNet();opt=torch.optim.Adam(net.parameters(),2e-3);start=0;print('START fresh, graph-disjoint 180/37/14',flush=True)
for ep in range(start,min(start+a.steps,MAX)):
 t0=time.time();net.train();perm=torch.randperm(len(X));tot=0.
 for b in range(0,len(X),2):
  ii=perm[b:b+2];x,y=X[ii],Y[ii]
  oy,ox=np.random.randint(0,S-CROP+1,2);x,y=x[...,oy:oy+CROP,ox:ox+CROP],y[...,oy:oy+CROP,ox:ox+CROP]
  if np.random.rand()<.5:x,y=x.flip(-1),y.flip(-1)
  if np.random.rand()<.5:x,y=x.flip(-2),y.flip(-2)
  k=np.random.randint(4);x,y=torch.rot90(x,k,(-2,-1)),torch.rot90(y,k,(-2,-1))
  opt.zero_grad();loss=F.cross_entropy(net(x),y.long(),weight=w);loss.backward();opt.step();tot+=loss.item()
 vloss=None
 if (ep+1)%5==0:
  net.eval();total=0.
  with torch.no_grad():
   for b in range(0,len(VX),2):total+=F.cross_entropy(net(VX[b:b+2]),VY[b:b+2].long(),weight=w,reduction='mean').item()*len(VX[b:b+2])
  vloss=total/len(VX)
  if vloss<best:best=vloss;torch.save({'net':net.state_dict(),'epoch':ep,'val_loss':vloss,'source_manifest_sha256':manifest['source_manifest_sha256']},BEST)
 hist.append({'epoch':ep+1,'train_loss_sum':tot,'val_loss':vloss,'best_val_loss':best,'seconds':round(time.time()-t0,1)})
 torch.save({'net':net.state_dict(),'opt':opt.state_dict(),'epoch':ep,'torch_rng':torch.get_rng_state(),'numpy_rng':np.random.get_state(),'best':best,'history':hist,'source_manifest_sha256':manifest['source_manifest_sha256']},CK)
 LOG.write_text(json.dumps(hist,indent=2)+'\n')
 print(hist[-1],flush=True)
print('COMPLETE',len(hist),'of',MAX,flush=True)
