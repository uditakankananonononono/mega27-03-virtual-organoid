"""Same-accession, post-result well balance audit on fixed original DIS blocks."""
import hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,'src')
from blocked_size import assess
source=Path('data/raw/orgasegment/ResearchData/DIS_database.csv')
derived=Path('data/raw/orgasegment/dis_merged_A0.csv')
x=pd.read_csv(source);d=pd.read_csv(derived)
keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
a=x[x.t==0][keys+['size']].rename(columns={'size':'A0'})
b=x[x.t==1][keys+['size']].rename(columns={'size':'A1'})
assert len(a)==len(b)==47687 and not a.duplicated(keys).any() and not b.duplicated(keys).any()
z=a.merge(b,on=keys,validate='one_to_one');base=assess(d)
assert base['n_eligible_blocks']==31 and base['n_eligible_donors']==12
blocks={(v['donor'],v['experiment'],float(v['forskolin_concentration_µM'])) for v in base['blocks']}
z=z[z[['donor','experiment','forskolin_concentration_µM']].apply(tuple,axis=1).isin(blocks)].copy()
z=z[z.condition.isin(['DMSO','VX445_VX661_VX770'])&(z.A0>0)&np.isfinite(z.A0)].copy()
z['size_bin']=np.select([z.A0<722,z.A0>=1400],['small','large'],default='middle')
z=z[z.size_bin!='middle'].copy()
z['observed']=((z.A1>0)&np.isfinite(z.A1))
z['log_response']=np.where(z.observed,np.log(z.A1/z.A0),np.nan)
group=['donor','experiment','forskolin_concentration_µM','condition','size_bin']
well=group+['well']
w=z.groupby(well,dropna=False).agg(n_baseline=('observed','size'),n_observed=('observed','sum'),observability=('observed','mean'),mean_log_response=('log_response','mean')).reset_index()
cell=w.groupby(group,dropna=False).agg(n_wells=('well','size'),n_baseline=('n_baseline','sum'),n_observed=('n_observed','sum'),equal_well_observability=('observability','mean'),n_wells_observed=('n_observed',lambda v:int((v>0).sum()))).reset_index()
obj=z.groupby(group,dropna=False).agg(object_mean=('log_response','mean'),object_observability=('observed','mean')).reset_index()
cell=cell.merge(obj,on=group,validate='one_to_one')
valid_w=w[w.n_observed>0].groupby(group,dropna=False).mean_log_response.mean().reset_index(name='equal_well_mean')
cell=cell.merge(valid_w,on=group,how='left',validate='one_to_one')
assert len(cell)==124
cv=w.groupby(group,dropna=False).n_observed.agg(lambda v:float(np.std(v,ddof=0)/np.mean(v)) if np.mean(v)>0 else None).reset_index(name='object_count_cv_rekey')
cell=cell.merge(cv,on=group,validate='one_to_one').rename(columns={'object_count_cv_rekey':'object_count_cv'})
signs={'VX445_VX661_VX770:large':1,'DMSO:large':-1,'VX445_VX661_VX770:small':-1,'DMSO:small':1}
rows=[]
for block in base['blocks']:
 g=cell[(cell.donor==block['donor'])&(cell.experiment==block['experiment'])&(cell['forskolin_concentration_µM']==float(block['forskolin_concentration_µM']))]
 assert len(g)==4
 entry={'donor':block['donor'],'experiment':block['experiment'],'dose':float(block['forskolin_concentration_µM']),'cells':{}}
 for r in g.itertuples(index=False):
  k=r.condition+':'+r.size_bin;entry['cells'][k]={col:getattr(r,col) for col in ['n_wells','n_baseline','n_observed','n_wells_observed','object_mean','equal_well_mean','object_observability','equal_well_observability','object_count_cv']}
 assert set(entry['cells'])==set(signs)
 for col,label in [('object_mean','effect_object'),('equal_well_mean','effect_equal_well'),('object_observability','observable_object'),('equal_well_observability','observable_equal_well')]:
  entry[label]=float(sum(signs[k]*entry['cells'][k][col] for k in signs)) if all(np.isfinite(entry['cells'][k][col]) for k in signs) else None
 rows.append(entry)
def aggregate(col):
 v=pd.DataFrame([{k:r[k] for k in ['donor','experiment','dose',col]} for r in rows]);v=v.rename(columns={col:'effect'})
 h=v.groupby('donor').effect.median();return {'n_blocks':int(v.effect.notna().sum()),'n_donors':int(h.notna().sum()),'equal_donor_mean_median':float(h.mean()),'positive_donors':int((h>0).sum()),'donor_medians':{k:float(t) for k,t in h.items()}}
obj=aggregate('effect_object');eq=aggregate('effect_equal_well')
assert obj['n_blocks']==31 and eq['n_blocks']==31
prior=json.load(open('results/tabular_track_retention.json'))
assert abs(obj['equal_donor_mean_median']-prior['source_paired_positive_area']['equal_donor_mean_median'])<1e-10
assert obj['positive_donors']==prior['source_paired_positive_area']['positive_donors']
flag=(abs(eq['equal_donor_mean_median']-obj['equal_donor_mean_median'])>=.10 or eq['positive_donors']!=obj['positive_donors'] or eq['n_blocks']<31)
out={'status':'post-result same-accession well-balance sensitivity; no missing-outcome recovery','protocol':'notes/postresult_well_balance_protocol.md','source_url':'https://zenodo.org/records/10610438','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'derived_sha256':hashlib.sha256(derived.read_bytes()).hexdigest(),'n_original_blocks':31,'n_donors':12,'n_baseline_detected':len(z),'n_observed_endpoint':int(z.observed.sum()),'object_weighted':obj,'equal_well':eq,'delta_equal_well_minus_object':float(eq['equal_donor_mean_median']-obj['equal_donor_mean_median']),'observability_object_weighted':aggregate('observable_object'),'observability_equal_well':aggregate('observable_equal_well'),'predefined_sensitivity_flag':bool(flag),'blocks':rows,'limits':'Wells are nested technical replicates; endpoint-negative objects absent from log-response analysis, never-detected baseline objects outside frame. Well weighting cannot repair size-selection, segmenter error or original sign/specificity failures.'}
Path('results/well_balance_sensitivity.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print({k:out[k] for k in ['n_baseline_detected','n_observed_endpoint','object_weighted','equal_well','delta_equal_well_minus_object','observability_object_weighted','observability_equal_well','predefined_sensitivity_flag']})
