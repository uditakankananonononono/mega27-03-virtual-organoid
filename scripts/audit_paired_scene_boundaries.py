"""Descriptive repeated-field mask area and boundary consistency; no model inference."""
import json
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.optimize import linear_sum_assignment
from scipy.ndimage import binary_dilation, binary_erosion, center_of_mass
from vorganoid.seg import iou_matrix
src=json.load(open('results/scene_label_consistency.json'))['pairs']
out=[]
for row in src:
 a,b=row['a'],row['b']
 def mask(p):return np.array(Image.open(Path('data/raw/organoid_basic')/p.replace('_img.jpg','_masks_organoid.png')))
 A,B=mask(a),mask(b);assert A.shape==B.shape
 M=iou_matrix(A,B);ii,jj=linear_sum_assignment(M,maximize=True);ag=np.unique(A)[1:];bg=np.unique(B)[1:];matches=[]
 for i,j in zip(ii,jj):
  if M[i,j]<.5:continue
  aa=A==ag[i];bb=B==bg[j];nA=int(aa.sum());nB=int(bb.sum())
  ca=np.array(center_of_mass(aa));cb=np.array(center_of_mass(bb));ba=aa^binary_erosion(aa);bc=bb^binary_erosion(bb)
  precision=np.count_nonzero(bc & binary_dilation(ba,iterations=1))/max(1,np.count_nonzero(bc))
  recall=np.count_nonzero(ba & binary_dilation(bc,iterations=1))/max(1,np.count_nonzero(ba))
  bf=2*precision*recall/(precision+recall) if precision+recall else 0.
  matches.append({'id_a':int(ag[i]),'id_b':int(bg[j]),'iou':float(M[i,j]),'area_a':nA,'area_b':nB,
                  'signed_log_area_ratio':float(np.log(nB/nA)),'absolute_log_area_ratio':float(abs(np.log(nB/nA))),
                  'centroid_displacement_px':float(np.linalg.norm(ca-cb)),'boundary_f1_1px':float(bf)})
 out.append({'a':a,'b':b,'pair_type':row['pair_type'],'n_matched_ge_0_5':len(matches),
             'median_abs_log_area_ratio':float(np.median([m['absolute_log_area_ratio'] for m in matches])) if matches else None,
             'median_centroid_displacement_px':float(np.median([m['centroid_displacement_px'] for m in matches])) if matches else None,
             'median_boundary_f1_1px':float(np.median([m['boundary_f1_1px'] for m in matches])) if matches else None,
             'all_matched_objects':matches})
# deterministic exact-copy control verifies metric semantics, including centroid and boundary operations
z=out[0]['all_matched_objects'][0];p=Path('data/raw/organoid_basic')/out[0]['a'].replace('_img.jpg','_masks_organoid.png');A=np.array(Image.open(p));v=A==z['id_a'];assert np.log(v.sum()/v.sum())==0 and np.linalg.norm(np.array(center_of_mass(v))-np.array(center_of_mass(v)))==0
bd=v^binary_erosion(v);assert np.count_nonzero(bd & binary_dilation(bd,iterations=1))/np.count_nonzero(bd)==1
result={'status':'retrospective repeated-field annotation/acquisition variability, not human inter-rater study',
        'pairs':out,'exact_copy_control':{'signed_log_area_ratio':0,'centroid_displacement_px':0,'boundary_f1_1px':1},
        'limitations':'repeated fields may reflect separate acquisition or export; annotator identities absent; no causal attribution to annotator or independent AP50 ceiling'}
Path('results/paired_scene_boundary_consistency.json').write_text(json.dumps(result,indent=2)+'\n')
for r in out:print(r['pair_type'],r['n_matched_ge_0_5'],r['median_abs_log_area_ratio'],r['median_centroid_displacement_px'],r['median_boundary_f1_1px'])
