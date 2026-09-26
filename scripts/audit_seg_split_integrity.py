"""Audit image split provenance by hashes and near-duplicate perceptual similarity.
No model scoring or retuning; flags overlap for review, not proof of independent acquisition.
"""
import hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image
from vorganoid.seg import list_pairs
rows=[]
for sp in ['train','val','eval']:
 for im,mask in list_pairs(sp):
  x=Image.open(im).convert('L').resize((32,32),Image.Resampling.BILINEAR)
  y=np.asarray(x,dtype=float);y=(y-y.mean())/(y.std()+1e-9)
  rows.append({'split':sp,'image':Path(im).name,'mask':Path(mask).name,
               'sha256_image':hashlib.sha256(Path(im).read_bytes()).hexdigest(),
               'sha256_mask':hashlib.sha256(Path(mask).read_bytes()).hexdigest(),
               'thumb_z32':y.ravel()})
assert [sum(r['split']==s for r in rows) for s in ['train','val','eval']]==[184,35,12]
from collections import defaultdict
imagehash=defaultdict(set);maskhash=defaultdict(set)
for r in rows:imagehash[r['sha256_image']].add(r['split']);maskhash[r['sha256_mask']].add(r['split'])
assert not any(len(v)>1 for v in imagehash.values()) and not any(len(v)>1 for v in maskhash.values())
# Cosine near-duplicates of thumbnails; high score is only a lead, not a biological provenance match.
A=np.stack([r['thumb_z32'] for r in rows]);N=A/np.linalg.norm(A,axis=1)[:,None]
near=[]
for i in range(len(rows)):
 for j in range(i+1,len(rows)):
  if rows[i]['split']==rows[j]['split']:continue
  c=float(np.dot(N[i],N[j]))
  if c>=.95:near.append({'a':rows[i]['image'],'split_a':rows[i]['split'],'b':rows[j]['image'],'split_b':rows[j]['split'],'normalized_thumbnail_cosine':c})
near.sort(key=lambda z:-z['normalized_thumbnail_cosine'])
from collections import Counter
out={'status':'split/hash/perceptual scout; no model gate','split_counts':dict(Counter(r['split'] for r in rows)),
 'cross_split_identical_image_hashes':sum(len(v)>1 for v in imagehash.values()),
 'cross_split_identical_mask_hashes':sum(len(v)>1 for v in maskhash.values()),
 'near_duplicate_threshold_cosine':.95,'cross_split_near_duplicate_pairs':near,
 'limit':'no acquisition metadata or donor identity in image filenames; thumbnail similarity cannot prove independence or leak absence',
 'manifest_sha256':hashlib.sha256(''.join(sorted(r['split']+'/'+r['image']+r['sha256_image'] for r in rows)).encode()).hexdigest()}
Path('results/seg_split_integrity_v2.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v if k!='cross_split_near_duplicate_pairs' else len(v) for k,v in out.items()})
print('top similarities',near[:10])
