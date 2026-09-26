"""Post-result fixed displacement-bin sensitivity on observed source pairs only."""
import hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0,'src')
from blocked_size import assess
source=Path('data/raw/orgasegment/ResearchData/DIS_database.csv');derived=Path('data/raw/orgasegment/dis_merged_A0.csv')
x=pd.read_csv(source);d=pd.read_csv(derived);keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
a=x[x.t==0][keys+['size','x','y']].rename(columns={'size':'A0','x':'x0','y':'y0'})
b=x[x.t==1][keys+['size','x','y']].rename(columns={'size':'A1','x':'x1','y':'y1'})
assert len(a)==len(b)==47687 and not a.duplicated(keys).any() and not b.duplicated(keys).any()
allpairs=a.merge(b,on=keys,validate='one_to_one')
valid=allpairs[(allpairs.A0>0)&(allpairs.A1>0)&np.isfinite(allpairs.A0)&np.isfinite(allpairs.A1)].copy()
valid['displacement']=np.hypot(valid.x1-valid.x0,valid.y1-valid.y0)
assert np.isfinite(valid.displacement).all() and len(valid)==29648
max_all=float(valid.displacement.max());n_above50=int((valid.displacement>50).sum())
base=assess(d);assert base['n_eligible_blocks']==31 and base['n_eligible_donors']==12
fixed={(q['donor'],q['experiment'],float(q['forskolin_concentration_µM'])) for q in base['blocks']}
z=valid[valid[['donor','experiment','forskolin_concentration_µM']].apply(tuple,axis=1).isin(fixed)].copy()
z=z[z.condition.isin(['DMSO','VX445_VX661_VX770'])].copy()
z['size_bin']=np.select([z.A0<722,z.A0>=1400],['small','large'],default='middle');z=z[z.size_bin!='middle'].copy()
z['log_response']=np.log(z.A1/z.A0)
z['bin']=pd.cut(z.displacement,[0,10,25,40,50,np.inf],right=False,include_lowest=True,labels=['0-10','10-25','25-40','40-50','>50'])
# The exact ==50 boundary belongs to the last bin; report rather than discard it.
assert z['bin'].notna().all()
counts=[]
for (arm,size,bin_),g in z.groupby(['condition','size_bin','bin'],observed=True):
 counts.append({'arm':arm,'size_bin':size,'displacement_bin':str(bin_),'n':len(g),'median_log_response':float(g.log_response.median()),'mean_log_response':float(g.log_response.mean())})
summary=[]
for (arm,size),g in z.groupby(['condition','size_bin']):
 summary.append({'arm':arm,'size_bin':size,'n':len(g),'median_displacement':float(g.displacement.median()),'n_ge40':int((g.displacement>=40).sum()),'fraction_ge40':float((g.displacement>=40).mean())})
blocks=[]
for bl in base['blocks']:
 g=z[(z.donor==bl['donor'])&(z.experiment==bl['experiment'])&(z['forskolin_concentration_µM']==float(bl['forskolin_concentration_µM']))]
 row={'donor':bl['donor'],'experiment':bl['experiment'],'dose':float(bl['forskolin_concentration_µM']),'cells':{}}
 for arm in ['DMSO','VX445_VX661_VX770']:
  for size in ['small','large']:
   cell=g[(g.condition==arm)&(g.size_bin==size)];low=cell[cell.displacement<40]
   row['cells'][arm+':'+size]={'n_all':len(cell),'n_under40':len(low),'mean_all':float(cell.log_response.mean()) if len(cell) else None,'mean_under40':float(low.log_response.mean()) if len(low) else None}
 assert all(c['n_all']>=3 for c in row['cells'].values())
 row['eligible_under40']=all(c['n_under40']>=3 for c in row['cells'].values())
 signs={'VX445_VX661_VX770:large':1,'DMSO:large':-1,'VX445_VX661_VX770:small':-1,'DMSO:small':1}
 row['effect_all']=sum(signs[k]*row['cells'][k]['mean_all'] for k in signs)
 row['effect_under40']=sum(signs[k]*row['cells'][k]['mean_under40'] for k in signs) if row['eligible_under40'] else None
 blocks.append(row)
def aggregate(effect):
 v=pd.DataFrame([{'donor':r['donor'],'effect':r[effect]} for r in blocks]);med=v.groupby('donor').effect.median()
 return {'n_blocks':int(v.effect.notna().sum()),'n_donors':int(med.notna().sum()),'equal_donor_mean_median':float(med.mean()),'positive_donors':int((med>0).sum()),'donor_medians':{k:float(t) for k,t in med.items()}}
allscore=aggregate('effect_all');low40=aggregate('effect_under40')
matched=pd.DataFrame([{'donor':r['donor'],'effect':r['effect_all']} for r in blocks if r['eligible_under40']]);matched_med=matched.groupby('donor').effect.median()
matched_all={'n_blocks':len(matched),'n_donors':len(matched_med),'equal_donor_mean_median':float(matched_med.mean()),'positive_donors':int((matched_med>0).sum()),'donor_medians':{k:float(t) for k,t in matched_med.items()}}

prior=json.load(open('results/tabular_track_retention.json'))
assert abs(allscore['equal_donor_mean_median']-prior['source_paired_positive_area']['equal_donor_mean_median'])<1e-10
control=[]
for size in ['small','large']:
 g=z[(z.condition=='DMSO')&(z.size_bin==size)]
 s=spearmanr(g.displacement,g.log_response)
 wells=[]
 for k,h in g.groupby(['donor','experiment','well','forskolin_concentration_µM']):
  if len(h)<5:continue
  r=spearmanr(h.displacement,h.log_response)
  if np.isfinite(r.statistic):wells.append({'donor':k[0],'experiment':k[1],'well':k[2],'dose':float(k[3]),'n':len(h),'rho':float(r.statistic)})
 donor=pd.DataFrame(wells).groupby('donor').rho.median()
 control.append({'size_bin':size,'pooled_n':len(g),'pooled_spearman':float(s.statistic),'n_wells_with_ge5':len(wells),'n_donors':len(donor),'median_well_rho_by_donor':{k:float(v) for k,v in donor.items()},'equal_donor_median_well_rho':float(donor.mean())})
out={'status':'post-result observed positive-pair displacement sensitivity; not missing-outcome recovery','protocol':'notes/postresult_displacement_sensitivity.md','source_url':'https://zenodo.org/records/10610438','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'derived_sha256':hashlib.sha256(derived.read_bytes()).hexdigest(),'n_positive_pairs_all_source':len(valid),'max_displacement_all_source':max_all,'n_above50_all_source':n_above50,'n_original_blocks':31,'n_observed_pairs_original_extreme_cells':len(z),'displacement_bins':counts,'arm_size_summary':summary,'all_positive_pair_score':allscore,'under40_score':low40,'all_pairs_on_under40_eligible_blocks':matched_all,'delta_under40_minus_matched_all':float(low40['equal_donor_mean_median']-matched_all['equal_donor_mean_median']),'delta_under40_minus_all_unmatched':float(low40['equal_donor_mean_median']-allscore['equal_donor_mean_median']),'excluded_blocks_under40':[{'donor':r['donor'],'experiment':r['experiment'],'dose':r['dose']} for r in blocks if not r['eligible_under40']],'dmso_control':control,'blocks':blocks,'limits':'Observed source pairs only; no x/y for missing endpoints. A <=50 radius may be tracking-algorithm support but source code not established. Restriction to <40 is a sensitivity, not a new eligibility gate; within-radius mistakes and never-detected baseline objects remain. No original biological gate repair.'}
Path('results/displacement_sensitivity.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print('max',max_all,'above50',n_above50,'n',len(z),'summaries',summary)
print('all',allscore['equal_donor_mean_median'],allscore['positive_donors'],'under40',low40['equal_donor_mean_median'],low40['positive_donors'],'matched_all',matched_all['equal_donor_mean_median'],matched_all['positive_donors'],'excluded',out['excluded_blocks_under40'])
print('controls',[(v['size_bin'],v['pooled_spearman'],v['equal_donor_median_well_rho']) for v in control])
