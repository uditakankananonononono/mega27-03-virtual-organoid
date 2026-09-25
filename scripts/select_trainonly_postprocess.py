"""Select fixed-grid postprocessing on disjoint val only, no eval image access."""
import itertools,json,sys,time,argparse
from pathlib import Path
import numpy as np,torch,torch.nn.functional as F
from scipy import ndimage as ndi
from skimage.segmentation import watershed
sys.path.insert(0,'src')
from vorganoid.seg import UNet,list_pairs,load_pair,average_precision

def instances(p,fg_t,seed_t,bnd_t,min_size):
    fg=(p[1]+p[2])>fg_t
    seeds=(p[1]>seed_t)&(p[2]<bnd_t)
    markers,_=ndi.label(seeds)
    lab=watershed(-p[1],markers,mask=fg)
    sizes=np.bincount(lab.ravel());small=np.where(sizes<min_size)[0]
    lab[np.isin(lab,small[small>0])]=0
    return lab

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--start",type=int,required=True);parser.add_argument("--count",type=int,default=8);a=parser.parse_args()
    torch.set_num_threads(1);S=512
    ck=torch.load('results/trainonly_seg/best.pt',map_location='cpu',weights_only=False)
    assert ck['epoch']==39
    net=UNet();net.load_state_dict(ck['net']);net.eval()
    # The old tune grid was selected on first 20 val images. This new model never trained on val.
    pairs=list_pairs('val')[:20]
    assert len(pairs)==20
    cache=[]
    for im,msk in pairs:
        x,_,lab=load_pair(im,msk,S);xt=torch.tensor(x)[None,None]
        with torch.no_grad():
            p=(torch.softmax(net(xt),1)+torch.softmax(net(xt.flip(-1)),1).flip(-1)+torch.softmax(net(xt.flip(-2)),1).flip(-2))/3
        pu=F.interpolate(p,size=lab.shape,mode='bilinear')[0].numpy()
        cache.append((pu,lab))
        print('cache val',Path(im).name,flush=True)
    grid=list(itertools.product([.4,.5,.6],[.5,.6,.7,.8],[.2,.3,.5],[10,20,40]))
    best=(-1,None);scores=[]
    assert a.start>=0 and 1<=a.count<=8 and a.start+a.count<=len(grid)
    for k in range(a.start,a.start+a.count):
        g=grid[k]
        aps=[average_precision(lab,instances(p,*g[:3],int(g[3]*(lab.shape[0]/S)**2)),.5) for p,lab in cache]
        sc=float(np.mean(aps));scores.append({'params':list(g),'val_mAP50':sc})
        if sc>best[0]:best=(sc,g)
        print('evaluated val grid',k+1,'/',len(grid),'score',sc,'best',best,flush=True)
    out={'dataset':'Zenodo 10278229 train-only fit / val-only postprocessing','selected_epoch':ck['epoch']+1,'selected_val_loss':ck['val_loss'],'val_grid_images':len(cache),'n_candidates':len(grid),'range':[a.start,a.start+a.count],'chunk_best_params':list(best[1]),'chunk_best_val_mAP50':best[0],'grid':scores,'eval_untouched':True}
    Path(f'results/trainonly_seg/val_grid_{a.start:03d}_{a.start+a.count:03d}.json').write_text(json.dumps(out,indent=2))
    print('selected',best,flush=True)
if __name__=='__main__':main()
