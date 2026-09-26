"""Exploratory VO3-R3 donor influence and within-cell robust summary, no new gate."""
import json,hashlib,sys
import numpy as np,pandas as pd
from scipy.stats import trim_mean
sys.path.insert(0,'src')
from blocked_size import assess
source='data/raw/orgasegment/dis_merged_A0.csv';df=pd.read_csv(source)
base=assess(df);blocks=base['blocks'];assert len(blocks)==31 and base['n_eligible_donors']==12
b=pd.DataFrame([{'donor':v['donor'],'experiment':v['experiment'],'dose':float(v['forskolin_concentration_µM']),
                 'effect':v['effect'],**v['cell_counts']} for v in blocks]); assert len(b)==31
median=b.groupby('donor').effect.median().sort_index(); donors=list(median.index)
rng=np.random.default_rng(20260926); boot=rng.integers(0,len(donors),size=(20000,len(donors)))
boot_mean=median.to_numpy()[boot].mean(axis=1)
loo={d:float(median.drop(d).mean()) for d in donors}
# Same previously selected blocks. No search over cutoffs or eligibility.
y=df[(df['forskolin_concentration_µM']>0)&df.condition.isin(['DMSO','VX445_VX661_VX770'])&(df.swelling>0)&np.isfinite(df.swelling)&(df.A0>0)&np.isfinite(df.A0)].copy()
y['log_swelling']=np.log(y.swelling)
keys={(v['donor'],v['experiment'],float(v['forskolin_concentration_µM'])) for v in blocks}
y=y[[tuple(t) in keys for t in y[['donor','experiment','forskolin_concentration_µM']].itertuples(index=False,name=None)]]
def effect(grp,lower,upper):
 z=grp.copy();z['size']=np.select([z.A0<lower,z.A0>=upper],['small','large'],default='middle');z=z[z['size']!='middle']
 g=z.groupby(['condition','size']).log_swelling
 desired=[(c,s) for c in ['DMSO','VX445_VX661_VX770'] for s in ['small','large']]
 if any(t not in g.groups or len(g.get_group(t))<3 for t in desired): return None
 m={t:float(g.get_group(t).mean()) for t in desired};mt={t:float(trim_mean(g.get_group(t),.2)) for t in desired}
 e=lambda q:(q['VX445_VX661_VX770','large']-q['DMSO','large'])-(q['VX445_VX661_VX770','small']-q['DMSO','small'])
 return {'ordinary':e(m),'trimmed20':e(mt),'cell_counts':{':'.join(t):len(g.get_group(t)) for t in desired}}
rows=[]
for key,g in y.groupby(['donor','experiment','forskolin_concentration_µM']):
 if key in keys:
  v=effect(g,722,1400);assert v is not None; rows.append({'donor':key[0],'experiment':key[1],'dose':key[2],**v})
assert len(rows)==31
for row in rows:
 orig=b[(b.donor==row['donor'])&(b.experiment==row['experiment'])&(b.dose==row['dose'])]
 assert len(orig)==1 and abs(orig.effect.iloc[0]-row['ordinary'])<1e-10
s=pd.DataFrame(rows);alt=s.groupby('donor').trimmed20.median().sort_index()
boundary={}
for name,lo,hi in [('lower10',648,1260),('higher10',794,1540)]:
 t=[]
 for key,g in y.groupby(['donor','experiment','forskolin_concentration_µM']):
  if key in keys:
   v=effect(g,lo,hi)
   if v is not None:t.append({'donor':key[0],'ordinary':v['ordinary']})
 dd=pd.DataFrame(t);dm=dd.groupby('donor').ordinary.median()
 boundary[name]={'cutoffs':[lo,hi],'eligible_prior_blocks':len(t),'eligible_donors':len(dm),'equal_donor_mean':float(dm.mean()) if len(dm) else None,'positive_donors':int((dm>0).sum()) if len(dm) else 0}
counts=[n for r in rows for n in r['cell_counts'].values()]
out={'status':'EXPLORATORY same-accession diagnostics, not a confirmatory discovery','input_sha256':hashlib.sha256(open(source,'rb').read()).hexdigest(),
 'original_cutoffs':[722,1400],'n_blocks':len(rows),'n_donors':len(median),'n_cell_rows':len(counts),
 'donor_median_effects':{k:float(v) for k,v in median.items()},'equal_donor_mean_median':float(median.mean()),
 'donor_bootstrap_seed':20260926,'donor_bootstrap_draws':len(boot_mean),'donor_bootstrap_ci95':[float(x) for x in np.quantile(boot_mean,[.025,.975])],
 'leave_one_donor_out_equal_means':loo,'trimmed20_equal_donor_mean_median':float(alt.mean()),
 'trimmed20_per_donor_median':{k:float(v) for k,v in alt.items()},'cell_count_range':[min(counts),max(counts)],
 'cell_count_median':float(np.median(counts)),'boundary_sensitivity':boundary,
 'not_testable_without_raw_tracks':['segmentation erosion/dilation/manual area','excluded-track retention bias']}
open('results/measurement_invariance.json','w').write(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['donor_median_effects','trimmed20_per_donor_median','leave_one_donor_out_equal_means']},indent=2))
