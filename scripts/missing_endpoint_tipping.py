"""Post-result hypothetical A1-missing tipping curve in original donor blocks."""
import hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,'src')
from blocked_size import assess
src=Path('data/raw/orgasegment/ResearchData/DIS_database.csv');d=pd.read_csv('data/raw/orgasegment/dis_merged_A0.csv')
x=pd.read_csv(src);keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
a=x[x.t==0][keys+['size']].rename(columns={'size':'A0'});b=x[x.t==1][keys+['size']].rename(columns={'size':'A1'})
assert len(a)==len(b)==47687
z=a.merge(b,on=keys,validate='one_to_one')
orig=assess(d);assert orig['n_eligible_blocks']==31 and orig['n_eligible_donors']==12
fixed={(q['donor'],q['experiment'],float(q['forskolin_concentration_µM'])) for q in orig['blocks']}
z=z[z[keys[:2]+['forskolin_concentration_µM']].apply(tuple,axis=1).isin(fixed)].copy()
z=z[(z.condition.isin(['DMSO','VX445_VX661_VX770']))&(z.A0>0)&np.isfinite(z.A0)].copy()
z['size_bin']=np.select([z.A0<722,z.A0>=1400],['small','large'],default='middle');z=z[z.size_bin!='middle'].copy()
z['observed_A1']=(z.A1>0)&np.isfinite(z.A1)
z['observed_log_ratio']=np.where(z.observed_A1,np.log(z.A1/z.A0),np.nan)
ref={}
for (arm,size),g in z.groupby(['condition','size_bin']):
 v=g.observed_log_ratio.dropna().to_numpy();assert len(v)>0
 ref[f'{arm}:{size}']={'n_observed':int(g.observed_A1.sum()),'n_missing':int((~g.observed_A1).sum()),
                       'median':float(np.median(v)),'q05':float(np.quantile(v,.05)),'q95':float(np.quantile(v,.95))}
z['base']=z.observed_log_ratio
for key,r in ref.items():
 arm,size=key.split(':');sel=(z.condition==arm)&(z.size_bin==size)&(~z.observed_A1);z.loc[sel,'base']=r['median']
assert z.base.notna().all()
def score(delta):
 z['response']=z.base + delta*((z.condition=='VX445_VX661_VX770')&(z.size_bin=='small')&(~z.observed_A1))
 rows=[]
 for b in orig['blocks']:
  g=z[(z.donor==b['donor'])&(z.experiment==b['experiment'])&(z['forskolin_concentration_µM']==float(b['forskolin_concentration_µM']))]
  means={}
  for arm in ['DMSO','VX445_VX661_VX770']:
   for size in ['small','large']:
    h=g[(g.condition==arm)&(g.size_bin==size)]
    assert len(h)>=3
    means[arm,size]=float(h.response.mean())
  e=(means['VX445_VX661_VX770','large']-means['DMSO','large'])-(means['VX445_VX661_VX770','small']-means['DMSO','small'])
  rows.append({'donor':b['donor'],'experiment':b['experiment'],'dose':float(b['forskolin_concentration_µM']),'effect':e})
 med=pd.DataFrame(rows).groupby('donor').effect.median()
 return {'delta':float(delta),'equal_donor_mean_median':float(med.mean()),'positive_donors':int((med>0).sum()),
         'donor_medians':{k:float(v) for k,v in med.items()}}
curve=[score(round(float(q),2)) for q in np.arange(0,2.001,.01)]
mean_tip=next((q for q in curve if q['equal_donor_mean_median']<=0),None)
sign_tip=next((q for q in curve if q['positive_donors']<=6),None)
negative_control=score(-2)
out={'status':'post-result hypothetical missing-A1 tipping curve, no identification or new gate','protocol':'notes/postresult_missing_endpoint_tipping.md',
 'source_url':'https://zenodo.org/records/10610438','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
 'n_baseline_detected_fixed_small_large_cells':len(z),'n_original_blocks':31,'n_donors':12,
 'observed_reference_by_arm_size':ref,'reference_delta0':curve[0],
 'first_nonpositive_mean':mean_tip,'first_six_or_fewer_positive_donors':sign_tip,
 'negative_delta_direction_check':negative_control,'curve':curve,
 'limits':'Unobserved endpoint values are hypothetical pooled-cell medians plus one-arm delta; not identified by source. Different missingness patterns or never-detected baseline objects can alter inference. Same accession; no original gate repair.'}
Path('results/missing_endpoint_tipping.json').write_text(json.dumps(out,indent=2)+'\n')
print('counts', {k:(v['n_observed'],v['n_missing']) for k,v in ref.items()});print('reference',curve[0]);print('mean tip',mean_tip);print('sign tip',sign_tip);print('direction check',negative_control['equal_donor_mean_median'])
