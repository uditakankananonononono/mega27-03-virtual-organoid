"""Post-result conditional percentile envelopes for absent A1 in original size blocks."""
import json,hashlib,sys
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,'src')
from blocked_size import assess
source=Path('data/raw/orgasegment/ResearchData/DIS_database.csv');raw=pd.read_csv(source);d=pd.read_csv('data/raw/orgasegment/dis_merged_A0.csv')
keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
a=raw[raw.t==0][keys+['size']].rename(columns={'size':'A0'});b=raw[raw.t==1][keys+['size']].rename(columns={'size':'A1'})
z=a.merge(b,on=keys,validate='one_to_one');base=assess(d);assert len(base['blocks'])==31
fixed={(q['donor'],q['experiment'],float(q['forskolin_concentration_µM'])) for q in base['blocks']}
z=z[z[['donor','experiment','forskolin_concentration_µM']].apply(tuple,axis=1).isin(fixed)]
z=z[z.condition.isin(['DMSO','VX445_VX661_VX770'])&(z.A0>0)&np.isfinite(z.A0)].copy()
z['size_bin']=np.select([z.A0<722,z.A0>=1400],['small','large'],default='middle');z=z[z.size_bin!='middle'].copy()
z['observed']=((z.A1>0)&np.isfinite(z.A1));z['log_response']=np.where(z.observed,np.log(z.A1/z.A0),np.nan)
reference={}
for (arm,size),g in z.groupby(['condition','size_bin']):
 v=g.log_response.dropna().values;reference[arm+':'+size]={'n_missing':int((~g.observed).sum()),'n_observed':int(g.observed.sum()),
   'q05':float(np.quantile(v,.05)),'q95':float(np.quantile(v,.95))}
# Coefficient signs for e=treated_large - DMSO_large - treated_small + DMSO_small.
signs={'VX445_VX661_VX770:large':1,'DMSO:large':-1,'VX445_VX661_VX770:small':-1,'DMSO:small':1}
z['lower']=z.log_response;z['upper']=z.log_response
for key,r in reference.items():
 arm,size=key.split(':');ix=(z.condition==arm)&(z.size_bin==size)&(~z.observed)
 z.loc[ix,'lower']=r['q05'] if signs[key]>0 else r['q95']
 z.loc[ix,'upper']=r['q95'] if signs[key]>0 else r['q05']
assert z.lower.notna().all() and z.upper.notna().all()
def aggregate(col):
 rows=[]
 for block in base['blocks']:
  g=z[(z.donor==block['donor'])&(z.experiment==block['experiment'])&(z['forskolin_concentration_µM']==float(block['forskolin_concentration_µM']))]
  means={}
  for arm in ['DMSO','VX445_VX661_VX770']:
   for size in ['small','large']:
    v=g[(g.condition==arm)&(g.size_bin==size)][col];assert len(v)>=3;means[arm+':'+size]=float(v.mean())
  e=sum(signs[key]*means[key] for key in signs)
  rows.append({'donor':block['donor'],'experiment':block['experiment'],'dose':float(block['forskolin_concentration_µM']),'effect':e})
 med=pd.DataFrame(rows).groupby('donor').effect.median()
 return {'equal_donor_mean_median':float(med.mean()),'positive_donors':int((med>0).sum()),'donor_medians':{k:float(v) for k,v in med.items()},'block_effects':rows}
lo=aggregate('lower');hi=aggregate('upper');assert lo['equal_donor_mean_median']<=hi['equal_donor_mean_median']
assert all(lo['donor_medians'][k]<=hi['donor_medians'][k] for k in lo['donor_medians'])
prior=json.load(open('results/tabular_track_retention.json'))
out={'status':'post-result conditional percentile-based missing-endpoint envelope, not identified biological truth',
 'protocol':'notes/postresult_missing_endpoint_bounds.md','source_url':'https://zenodo.org/records/10610438',
 'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'n_baseline_detected':len(z),'n_original_blocks':31,
 'percentile_reference':reference,'minimize_contrast':lo,'maximize_contrast':hi,'zero_contained':bool(lo['equal_donor_mean_median']<=0<=hi['equal_donor_mean_median']),
 'observed_positive_area_pair_mean':prior['source_paired_positive_area']['equal_donor_mean_median'],
 'derived_included_mean':prior['derived_original']['equal_donor_mean_median'],
 'limits':'5th/95th observed quantiles are hypothetical bounds for missing A1, not verified support; all missing outcomes may differ. Baseline-undetected objects outside source; no new p-value or gate.'}
Path('results/missing_endpoint_bounds.json').write_text(json.dumps(out,indent=2)+'\n')
print('frame',len(z),'lower',lo['equal_donor_mean_median'],lo['positive_donors'],'upper',hi['equal_donor_mean_median'],hi['positive_donors'],'zero contained',out['zero_contained']);print(reference)
