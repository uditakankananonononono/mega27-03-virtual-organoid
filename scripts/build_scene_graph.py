"""Image-only candidate scene groups; retrospective forensic use, no model scores."""
import glob,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from PIL import Image
ROOT=Path('data/raw/organoid_basic')
rows=[]
for split in ['train','val','eval']:
 for path in sorted(ROOT.joinpath(split).glob('*_img.jpg')):
  px=np.asarray(Image.open(path).convert('L').resize((32,32),Image.Resampling.BILINEAR),dtype=np.float64)
  px-=px.mean(); px/=px.std()+1e-9; px=px.ravel();px/=np.linalg.norm(px)
  rows.append({'split':split,'name':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'vector':px})
assert Counter(r['split'] for r in rows)=={'train':184,'val':35,'eval':12}
X=np.stack([r['vector'] for r in rows]); sims=X@X.T

def graph(threshold):
 parent=list(range(len(rows)))
 def find(x):
  while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
  return x
 def join(i,j):
  a,b=find(i),find(j)
  if a!=b:parent[max(a,b)]=min(a,b)
 edges=[]
 for i in range(len(rows)):
  for j in range(i+1,len(rows)):
   exact=rows[i]['sha256']==rows[j]['sha256']
   if sims[i,j]>=threshold or exact:
    join(i,j)
    edges.append({'a':rows[i]['split']+'/'+rows[i]['name'],'b':rows[j]['split']+'/'+rows[j]['name'],
                  'cosine':round(float(sims[i,j]),9),'same_sha256':exact})
 comps=defaultdict(list)
 for i in range(len(rows)):comps[find(i)].append(i)
 multi=[]
 for ids in comps.values():
  if len(ids)<2:continue
  mult={'members':[rows[i]['split']+'/'+rows[i]['name'] for i in ids],
        'split_counts':dict(Counter(rows[i]['split'] for i in ids))}
  multi.append(mult)
 multi.sort(key=lambda z:z['members'][0])
 return {'threshold':threshold,'n_components':len(comps),'n_multiframe_components':len(multi),
         'n_cross_split_components':sum(len(x['split_counts'])>1 for x in multi),
         'cross_split_component_members':sum(len(x['members']) for x in multi if len(x['split_counts'])>1),
         'edges':edges,'multi_frame_components':multi}
allg=[graph(t) for t in [.95,.97,.99]]
# A seeded negative-control sample establishes scale, not a null p-value.
rng=np.random.default_rng(20260926)
pairs=[(i,j) for i in range(len(rows)) for j in range(i+1,len(rows)) if rows[i]['split']!=rows[j]['split']]
chosen=rng.choice(len(pairs),size=1000,replace=False)
control=np.array([sims[pairs[z][0],pairs[z][1]] for z in chosen])
manifest=''.join(sorted(r['split']+'/'+r['name']+r['sha256'] for r in rows))
out={'status':'retrospective image-only candidate scene grouping, not a repaired held-out evaluation',
     'image_count':len(rows),'split_counts':dict(Counter(r['split'] for r in rows)),
     'image_manifest_sha256':hashlib.sha256(manifest.encode()).hexdigest(),
     'method':'normalized 32x32 grayscale bilinear thumbnail cosine + exact image SHA256; undirected edges and connected components; no mask, label, prediction or AP50 input',
     'negative_control':{'seed':20260926,'n_cross_split_random_pairs':1000,
                         'cosine_median':float(np.median(control)),
                         'cosine_q99':float(np.quantile(control,.99)),
                         'cosine_max':float(np.max(control))},
     'threshold_sensitivity':allg,
     'limits':'Threshold chosen after old AP50 was seen. Components are candidate scenes, not verified donor/well IDs. Missing edges do not establish independence. The old eval cannot be resealed.'}
Path('results/scene_graph_forensics.json').write_text(json.dumps(out,indent=2)+'\n')
print('manifest',out['image_manifest_sha256'],'negative control',out['negative_control'])
for g in allg:print(g['threshold'],g['n_components'],g['n_multiframe_components'],g['n_cross_split_components'],g['cross_split_component_members'])
