"""Frozen candidate-scene group reassignment per notes/prereg_group_disjoint_seg.md."""
import hashlib,json,shutil
from pathlib import Path
from collections import Counter
import numpy as np
from PIL import Image

ROOT=Path('data/raw/organoid_basic'); OUT=Path('data/derived/group_disjoint_seg')
FROZEN='ba003cd9645067eb455f5a2d0166a53a0a41aa73378f938e848f2ca9c73a9b65'
G=json.loads(Path('results/scene_graph_forensics.json').read_text())
assert G['image_manifest_sha256']==FROZEN and G['image_count']==231
rows=[]
for split in ('train','val','eval'):
 for path in sorted((ROOT/split).glob('*_img.jpg')):
  mask=path.with_name(path.name.replace('_img.jpg','_masks_organoid.png'))
  assert mask.is_file()
  rows.append((split,path,mask,hashlib.sha256(path.read_bytes()).hexdigest()))
assert Counter(s for s,*_ in rows)=={'train':184,'val':35,'eval':12}
manifest=''.join(sorted(s+'/'+p.name+h for s,p,m,h in rows))
assert hashlib.sha256(manifest.encode()).hexdigest()==FROZEN
by_name={s+'/'+p.name:(s,p,m,h) for s,p,m,h in rows}
graph=G['threshold_sensitivity'][0];assert graph['threshold']==0.95 and graph['n_components']==215
parent={k:k for k in by_name}
def find(k):
 while parent[k]!=k:parent[k]=parent[parent[k]];k=parent[k]
 return k
for e in graph['edges']:
 a,b=find(e['a']),find(e['b']);parent[max(a,b)]=min(a,b)
components={}
for k in by_name:components.setdefault(find(k),[]).append(k)
assert len(components)==graph['n_components']
assert sorted((sorted(v) for v in components.values() if len(v)>1))==sorted((sorted(v['members']) for v in graph['multi_frame_components']))
priority={'train':0,'val':1,'eval':2};assignment={}
for members in components.values():
 dest=max((by_name[k][0] for k in members),key=lambda x:priority[x])
 for k in members:assignment[k]=dest
assert len(assignment)==231
moved=[]; counts=Counter()
# Do not delete an existing destination: a rerun must not silently replace a trained experiment.
assert not OUT.exists(),f'{OUT} exists; inspect rather than overwrite'
for d in ('train','val','eval'):(OUT/d).mkdir(parents=True)
for key, dest in sorted(assignment.items()):
 src,p,m,h=by_name[key]
 for f in (p,m):
  target=OUT/dest/f.name
  assert not target.exists(),f'collision: {target}'
  target.symlink_to(f.resolve())
 counts[dest]+=1
 if src!=dest:moved.append({'image':key,'from':src,'to':dest,'component':sorted(components[find(key)])})
assert sum(counts.values())==231
for members in components.values():assert len({assignment[k] for k in members})==1
out={'status':'exploratory previously seen eval; no benchmark win','protocol':'notes/prereg_group_disjoint_seg.md','source_manifest_sha256':FROZEN,'threshold':0.95,'component_count':len(components),'split_counts':dict(counts),'n_moved':len(moved),'moved':moved,
     'dest':str(OUT),'limits':'Image-only group proxies, no donor/plate/field IDs. Previously scored eval; not independent.'}
Path('results/group_disjoint_seg_split.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k not in ('moved','limits')});print('moved',moved)
