"""Fixed donor-by-cell baseline-positive endpoint observability, no gate update."""
import hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,'src')
from blocked_size import assess
src=Path('data/raw/orgasegment/ResearchData/DIS_database.csv');der=Path('data/raw/orgasegment/dis_merged_A0.csv')
assert hashlib.sha256(src.read_bytes()).hexdigest()=='0e0e7913947a29c7b32ce87ca73a98614033d41ce63063a1d65053e766b2f315'
x=pd.read_csv(src);d=pd.read_csv(der); keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
a=x.loc[x.t==0,keys+['size']].rename(columns={'size':'A0'});b=x.loc[x.t==1,keys+['size']].rename(columns={'size':'A1'});assert len(a)==len(b)==47687
z=a.merge(b,on=keys,validate='one_to_one');q=assess(d);assert q['n_eligible_blocks']==31 and q['n_eligible_donors']==12
fixed={(v['donor'],v['experiment'],float(v['forskolin_concentration_µM'])) for v in q['blocks']}
z=z[z[['donor','experiment','forskolin_concentration_µM']].apply(tuple,axis=1).isin(fixed)]
z=z[z.condition.isin(['DMSO','VX445_VX661_VX770'])&np.isfinite(z.A0)&(z.A0>0)].copy()
z['size_bin']=np.select([z.A0<722,z.A0>=1400],['small','large'],default='middle');z=z[z.size_bin!='middle'].copy()
z['positive_A1']=np.isfinite(z.A1)&(z.A1>0)
audit=json.load(open('results/postresult_attrition_model_limit.json'))
assert len(z)==audit['n_baseline_detected']==6964 and int(z.positive_A1.sum())==audit['n_positive_A1']==4964
out={'status':'post-result donor-specific detection/attrition audit; not missing-value recovery or gate','protocol':'notes/postresult_donor_attrition_decomposition_protocol.md','source_url':'https://zenodo.org/records/10610438','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'donors':{},'n_baseline_detected':len(z),'n_positive_A1':int(z.positive_A1.sum()),'n_missing_A1':int((~z.positive_A1).sum()),'limits':'Same selected blocks and previously seen endpoints. Donor differences may confound source batch, image, or genotype; no direct mask audit or missing response identification. No original gate change.'}
for donor,g in z.groupby('donor'):
 cells={}
 for arm in ['DMSO','VX445_VX661_VX770']:
  for size in ['small','large']:
   h=g[(g.condition==arm)&(g.size_bin==size)];assert len(h)>=3
   cells[arm+':'+size]={'baseline_detected':len(h),'positive_A1':int(h.positive_A1.sum()),'missing_A1':int((~h.positive_A1).sum()),'observability_rate':float(h.positive_A1.mean())}
 gs={s:cells['VX445_VX661_VX770:'+s]['observability_rate']-cells['DMSO:'+s]['observability_rate'] for s in ['small','large']}
 out['donors'][donor]={'cells':cells,'n_baseline_detected':len(g),'n_positive_A1':int(g.positive_A1.sum()),'gap_treatment_minus_dmso':gs,'gap_interaction_large_minus_small':gs['large']-gs['small']}
assert len(out['donors'])==12
for key in ['baseline_detected','positive_A1','missing_A1']:
 assert sum(c[key] for don in out['donors'].values() for c in don['cells'].values())==out['n_'+key if key=='baseline_detected' else ('n_positive_A1' if key=='positive_A1' else 'n_missing_A1')]
r=[d['gap_interaction_large_minus_small'] for d in out['donors'].values()]
out['donor_gap_interaction_summary']={'min':float(np.min(r)),'median':float(np.median(r)),'max':float(np.max(r)),'positive':int(np.sum(np.array(r)>0)),'negative':int(np.sum(np.array(r)<0)),'zero':int(np.sum(np.array(r)==0))}
Path('results/postresult_donor_attrition_decomposition.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('fixed source frame',out['n_baseline_detected'],out['n_positive_A1'],out['n_missing_A1'],'donors',len(out['donors']),'summary',out['donor_gap_interaction_summary'])
for k,v in out['donors'].items():print(k,v['n_baseline_detected'],v['n_positive_A1'],round(v['gap_interaction_large_minus_small'],4))
