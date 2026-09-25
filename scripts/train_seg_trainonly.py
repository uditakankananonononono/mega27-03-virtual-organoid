"""Incremental train-only segmentation: one or two epochs per invocation, val loss only."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np,torch,torch.nn.functional as F
sys.path.insert(0,'src')
from vorganoid.seg import UNet,list_pairs,load_pair

p=argparse.ArgumentParser();p.add_argument('--steps',type=int,default=1);a=p.parse_args()
assert 1<=a.steps<=3
S=512;CROP=256;MAX=60
root=Path('results/trainonly_seg');root.mkdir(exist_ok=True)
CK=root/'resume.pt';BEST=root/'best.pt';LOG=root/'history.json'
torch.set_num_threads(1)
train=list_pairs('train');val=list_pairs('val');ev=list_pairs('eval')
assert (len(train),len(val),len(ev))==(184,35,12)
assert not ({i for i,_ in train}&{i for i,_ in val}) and not ({i for i,_ in train}&{i for i,_ in ev})
tr=[load_pair(i,m,S)[:2] for i,m in train]; va=[load_pair(i,m,S)[:2] for i,m in val]
X=torch.tensor(np.stack([t[0] for t in tr]))[:,None];Y=torch.tensor(np.stack([t[1].astype(np.uint8) for t in tr]));del tr
VX=torch.tensor(np.stack([t[0] for t in va]))[:,None];VY=torch.tensor(np.stack([t[1].astype(np.uint8) for t in va]));del va
net=UNet();opt=torch.optim.Adam(net.parameters(),2e-3);w=torch.tensor([1.,1.,3.]);best=float('inf');history=[]
if CK.exists():
    ck=torch.load(CK,map_location='cpu',weights_only=False);net.load_state_dict(ck['net']);opt.load_state_dict(ck['opt']);torch.set_rng_state(ck['torch_rng']);np.random.set_state(ck['numpy_rng']);start=ck['epoch']+1;best=ck['best'];history=ck['history'];print('RESUME',start,flush=True)
else:
    torch.manual_seed(2026);np.random.seed(2026);start=0;print('START fresh train only 184 val 35 eval untouched 12',flush=True)
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
        net.eval();s=0.
        with torch.no_grad():
            for b in range(0,len(VX),2):s+=F.cross_entropy(net(VX[b:b+2]),VY[b:b+2].long(),weight=w,reduction='mean').item()*len(VX[b:b+2])
        vloss=s/len(VX)
        if vloss<best:best=vloss;torch.save({'net':net.state_dict(),'epoch':ep,'val_loss':vloss},BEST)
    history.append({'epoch':ep+1,'train_loss_sum':tot,'val_loss':vloss,'best_val_loss':best,'seconds':round(time.time()-t0,1)})
    state={'net':net.state_dict(),'opt':opt.state_dict(),'epoch':ep,'torch_rng':torch.get_rng_state(),'numpy_rng':np.random.get_state(),'best':best,'history':history}
    torch.save(state,CK);LOG.write_text(json.dumps(history,indent=2))
    print(history[-1],flush=True)
print('COMPLETE',len(history),'of',MAX,flush=True)
