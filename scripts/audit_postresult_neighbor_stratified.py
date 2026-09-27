"""Retrospective within-image arm-size-stratified neighbor observability check."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
src=Path('data/raw/orgasegment/ResearchData/DIS_database.csv');sha=hashlib.sha256(src.read_bytes()).hexdigest()
assert sha=='0e0e7913947a29c7b32ce87ca73a98614033d41ce63063a1d65053e766b2f315'
raw=pd.read_csv(src); keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
a=raw[raw.t==0][keys+['size','x','y','image']].rename(columns={'size':'A0','x':'x0','y':'y0','image':'image0'})
b=raw[raw.t==1][keys+['size']].rename(columns={'size':'A1'})
assert len(a)==len(b)==47687 and not a.duplicated(keys).any() and not b.duplicated(keys).any()
fixed={(r['donor'],r['experiment'],float(r['forskolin_concentration_µM'])) for r in json.load(open('results/blocked_size.json'))['blocks']};assert len(fixed)==31
z=a.merge(b,on=keys,validate='one_to_one');z=z[z[['donor','experiment','forskolin_concentration_µM']].apply(tuple,axis=1).isin(fixed)&z.condition.isin(['DMSO','VX445_VX661_VX770'])&(z.A0>0)&((z.A0<722)|(z.A0>=1400))].copy()
z['size_bin']=np.where(z.A0<722,'small','large');z['observed']=((z.A1>0)&np.isfinite(z.A1)).astype(int)
assert len(z)==6964 and int(z.observed.sum())==4964 and z.image0.nunique()==143
lookup={}
for image,g in a[(a.A0>0)&np.isfinite(a.A0)&a.x0.notna()&a.y0.notna()].groupby('image0'):
 coords=g[['x0','y0']].to_numpy(float)
 if len(coords)<2:continue
 dist,_=cKDTree(coords).query(coords,k=2)
 for row,dd in zip(g.itertuples(index=False),dist):
  key=tuple(getattr(row,v) for v in keys);assert key not in lookup;lookup[key]=float(dd[1])
z['nn1']=[lookup[tuple(r[v] for v in keys)] for _,r in z.iterrows()]
old=json.load(open('results/postresult_neighbor_observability.json'))
assert np.allclose(np.quantile(z.nn1,[.25,.5,.75]),old['nn1_quartile_cutoffs_px'])
rows=[];small=[]
for (image,arm,size),g in z.groupby(['image0','condition','size_bin']):
 base={'image':image,'donor':g.donor.iloc[0],'arm':arm,'size':size,'n':len(g)}
 if len(g)<4 or g.nn1.nunique()<2:small.append(base);continue
 close=g.nn1<=g.nn1.median()
 if close.nunique()!=2:small.append(base);continue
 rows.append({**base,'close_n':int(close.sum()),'far_n':int((~close).sum()),'close_observed':int(g.loc[close,'observed'].sum()),'far_observed':int(g.loc[~close,'observed'].sum()),'close_minus_far':float(g.loc[close,'observed'].mean()-g.loc[~close,'observed'].mean())})
assert sum(x['n'] for x in rows)+sum(x['n'] for x in small)==6964
assert len({x['image'] for x in rows+small})==143
assert len({x['donor'] for x in rows})<=12
weights=np.array([x['n'] for x in rows]);diffs=np.array([x['close_minus_far'] for x in rows]);donors=sorted({x['donor'] for x in rows});donor_diffs={str(d):float(np.mean([x['close_minus_far'] for x in rows if x['donor']==d])) for d in donors}
groups=[]
for r in rows:
 g=z[(z.image0==r['image'])&(z.condition==r['arm'])&(z.size_bin==r['size'])]
 y=g.observed.to_numpy();close=(g.nn1<=g.nn1.median()).to_numpy();assert len(y)==r['n'];groups.append((y,close))
rng=np.random.default_rng(20260927);null=[]
for _ in range(100):
 v=[]
 for y,close in groups:
  c=rng.permutation(close);v.append(float(y[c].mean()-y[~c].mean()))
 null.append(float(np.mean(v)))
arm_size=[]
for arm in ['DMSO','VX445_VX661_VX770']:
 for size in ['small','large']:
  v=[r for r in rows if r['arm']==arm and r['size']==size]
  arm_size.append({'arm':arm,'size':size,'n_strata':len(v),'n_objects':sum(r['n'] for r in v),'equal_stratum_close_minus_far':float(np.mean([r['close_minus_far'] for r in v])) if v else None})
out={'status':'post-result within-field x arm x size stratified sensitivity, descriptive only','protocol':'notes/postresult_neighbor_stratified_sensitivity.md','source_url':'https://zenodo.org/records/10610438','source_sha256':sha,
 'n_frame':len(z),'n_positive_A1':int(z.observed.sum()),'n_images_frame':z.image0.nunique(),'n_eligible_strata':len(rows),'n_excluded_strata':len(small),'n_eligible_objects':sum(r['n'] for r in rows),'n_excluded_objects':sum(r['n'] for r in small),
 'equal_stratum_close_minus_far':float(diffs.mean()),'object_weighted_stratum_close_minus_far':float(np.average(diffs,weights=weights)),'per_arm_size':arm_size,
 'donor_differences':donor_diffs,'donor_equal_close_minus_far':float(np.mean(list(donor_diffs.values()))),'donor_negative_count':int(sum(v<0 for v in donor_diffs.values())),'donor_n':len(donor_diffs),
 'permutation_seed':20260927,'permutation_n':100,'null_mean':float(np.mean(null)),'null_95pct':[float(v) for v in np.quantile(null,[.025,.975])],
 'eligible_strata':rows,'excluded_strata':small,'limits':'Retrospective technical-source association, not causal crowding or biological replicates; image paths can be repeated fields; missing A1 response remains unknown; no gate repair.'}
Path('results/postresult_neighbor_stratified.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['eligible_strata','excluded_strata','donor_differences']},indent=2))
