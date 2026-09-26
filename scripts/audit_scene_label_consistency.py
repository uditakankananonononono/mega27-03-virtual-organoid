"""Retrospective repeated-field annotation overlap; not a model eval."""
import json
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.optimize import linear_sum_assignment
from vorganoid.seg import iou_matrix
s=json.load(open('results/scene_graph_forensics.json'))
g=next(z for z in s['threshold_sensitivity'] if z['threshold']==.99)
pairs=[(z['members'][0],z['members'][1],'image-only cosine >=0.99') for z in g['multi_frame_components'] if len(z['split_counts'])>1]
assert len(pairs)==3 and all('eval/' in a+' '+b or 'val/' in a+' '+b for a,b,_ in pairs)
pairs.append(('train/0003_3250648_ckmbz0c8v1vke345zguvwuom1_img.jpg','val/0003_6054432_cknj077t115y9345zmqfpq9s0_img.jpg','0.95 graph flag, visual distinct-layout control'))
pairs.append(('train/20210727_11104_ckufl4yfh3f200y77899mcbro_img.jpg','eval/20210727_26838_ckw3ig3xg5ngk0zbp0g7j7dqt_img.jpg','retrospectively selected distinct-layout cross-split control'))
out=[]
for a,b,kind in pairs:
 def get(p):return np.array(Image.open(Path('data/raw/organoid_basic')/p.replace('_img.jpg','_masks_organoid.png')))
 A,B=get(a),get(b);assert A.shape==B.shape and A.ndim==2
 M=iou_matrix(A,B); row,col=linear_sum_assignment(M,maximize=True);z=M[row,col]
 fg1=A>0;fg2=B>0
 out.append({'a':a,'b':b,'pair_type':kind,'n_objects_a':len(np.unique(A))-1,'n_objects_b':len(np.unique(B))-1,
             'one_to_one_matches_iou_ge_0_5':int(sum(z>=.5)),
             'one_to_one_matches_iou_ge_0_75':int(sum(z>=.75)),
             'median_iou_ge_0_5':float(np.median(z[z>=.5])) if any(z>=.5) else None,
             'foreground_dice_same_coordinates':float(2*(fg1&fg2).sum()/(fg1.sum()+fg2.sum()))})
result={'status':'post-result masks-only annotation consistency; no benchmark score or selection',
        'selected_from':'image-only >=0.99 graph for three candidate pairs; two explicitly named descriptive controls',
        'pairs':out,'limitations':'same-coordinate overlap may be reduced by field shift; object matches are not proof of annotation independence or donor identity; these test files have already been seen'}
Path('results/scene_label_consistency.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
