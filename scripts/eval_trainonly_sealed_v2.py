"""Single sealed v2 eval, frozen selection and exact author scorer commit checked before scoring."""
"""Single final eval of preregistered train-only model and val-only parameters."""
import importlib.util, json, sys
from pathlib import Path
import numpy as np, torch, torch.nn.functional as F
from PIL import Image
sys.path.insert(0, 'src')
from vorganoid.seg import UNet, list_pairs, load_pair, average_precision
from scipy import ndimage as ndi
from skimage.segmentation import watershed

def instances(p, fg_t, seed_t, bnd_t, min_size):
    fg=(p[1]+p[2])>fg_t
    seeds=(p[1]>seed_t)&(p[2]<bnd_t)
    markers,_=ndi.label(seeds)
    lab=watershed(-p[1],markers,mask=fg)
    sizes=np.bincount(lab.ravel())
    small=np.where(sizes<min_size)[0]
    lab[np.isin(lab,small[small>0])]=0
    return lab


def main():
    import subprocess,hashlib,glob,itertools
    assert len(sys.argv)==2
    assert subprocess.check_output(['git','-C',sys.argv[1],'rev-parse','HEAD'],text=True).strip()=='aeaec17c42fefde61a8b2215a616d419e8126e9f'
    author_path=Path(sys.argv[1])/'lib'/'metrics.py'
    spec=importlib.util.spec_from_file_location('author_metrics',author_path)
    author=importlib.util.module_from_spec(spec);spec.loader.exec_module(author)
    torch.set_num_threads(1)
    selection=json.load(open('results/trainonly_seg/val_selection_v2.json'))
    assert selection['eval_used'] is False and selection['selected_epoch']==30 and selection['selected_grid_index']==41 and selection['n_candidates']==108
    files=glob.glob('results/trainonly_seg/val_grid_v2_*.json'); assert sorted(files)==selection['grid_files']
    assert all(hashlib.sha256(open(f,'rb').read()).hexdigest()==selection['grid_file_sha256'][f] for f in files)
    rows={}
    for f in files:
        d=json.load(open(f)); a,b=d['range']; assert len(d['grid'])==b-a
        for i,row in zip(range(a,b),d['grid']): assert i not in rows; rows[i]=row
    assert sorted(rows)==list(range(108))
    assert rows[selection['selected_grid_index']]['params']==selection['selected_params']
    assert selection['selected_val_mAP50']==max(r['val_mAP50'] for r in rows.values())
    params=selection['selected_params']; assert params==[.5,.5,.3,40]
    ck=torch.load('results/trainonly_seg/best_v2.pt',map_location='cpu',weights_only=False);assert ck['epoch']==29
    S=512;net=UNet();net.load_state_dict(ck['net'])
    net.eval();out=[]
    pairs=list_pairs('eval'); assert len(pairs)==12
    for img,msk in pairs:
        x,_,lab=load_pair(img,msk,S)
        xt=torch.tensor(x)[None,None]
        with torch.no_grad():
            p=(torch.softmax(net(xt),1)+torch.softmax(net(xt.flip(-1)),1).flip(-1)+torch.softmax(net(xt.flip(-2)),1).flip(-2))/3
        pu=F.interpolate(p,size=lab.shape,mode='bilinear')[0].numpy()
        pr=instances(pu,*params[:3],int(params[3]*(lab.shape[0]/S)**2))
        # Author code uses maximum label as instance count. Watershed can leave ID gaps
        # when tiny regions are removed. Relabel IDs contiguously without changing pixels.
        ids=np.unique(pr); pr=np.searchsorted(ids,pr).astype(np.int32)
        ap_author,tp,fp,fn=author.average_precision(lab.astype(np.int32),pr,.5)
        ap_ours=average_precision(lab,pr,.5)
        row={'image':Path(img).name,'ap_author':float(ap_author[0]),'ap_ours':float(ap_ours),'tp_author':int(tp[0]),'fp_author':int(fp[0]),'fn_author':int(fn[0]),'n_gt':int(len(np.unique(lab))-1),'n_pred':int(len(np.unique(pr))-1),'max_gt_id':int(lab.max()),'max_pred_id':int(pr.max()),'gt_id_gaps':int(lab.max()-(len(np.unique(lab))-1))}
        print(row,flush=True);out.append(row)
    result={'author_release':'https://github.com/Living-Technologies/OrgaSegment/tree/v1.0.1',
            'author_release_tag':'aeaec17c42fefde61a8b2215a616d419e8126e9f',
            'our_model':'results/trainonly_seg/best_v2.pt',
            'selection_commit':'7eaa46c17a610ac8173ad1194d992668403a6760',
            'fixed_params':params,
            'per_image':out,'mean_author':float(np.mean([r['ap_author'] for r in out])),
            'mean_ours':float(np.mean([r['ap_ours'] for r in out])),
            'published_comparison':.76,
            'interpretation':'single sealed eval after train-only model and val-only grid; released Mask R-CNN was not run; original evaluator counts max GT ID and one eval mask has a missing instance ID'}
    Path('results/trainonly_seg/sealed_eval_v2.json').write_text(json.dumps(result,indent=2))
    print('means',result['mean_author'],result['mean_ours'],flush=True)
if __name__=='__main__': main()
