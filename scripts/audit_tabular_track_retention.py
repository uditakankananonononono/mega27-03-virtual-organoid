"""Post-result source-table retention and detection-score sensitivity for fixed DIS blocks."""
import hashlib,json
import numpy as np,pandas as pd
from pathlib import Path
from scipy.stats import mannwhitneyu
import sys
sys.path.insert(0,'src')
from blocked_size import assess
source=Path('data/raw/orgasegment/ResearchData/DIS_database.csv');derived=Path('data/raw/orgasegment/dis_merged_A0.csv')
x=pd.read_csv(source);d=pd.read_csv(derived)
keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
assert set(x.t.dropna().unique())=={0.,1.}
a=x[x.t==0][keys+['size','score']].rename(columns={'size':'A0_source','score':'score0'})
b=x[x.t==1][keys+['size','score']].rename(columns={'size':'A1_source','score':'score1'})
assert len(a)==len(b)==47687 and not a.duplicated(keys).any() and not b.duplicated(keys).any()
z=a.merge(b,on=keys,validate='one_to_one',indicator=False)
ids=d[keys+['A0','A1']].copy();ids['included_derived']=True
z=z.merge(ids,on=keys,how='left',validate='one_to_one');z['included_derived']=z.included_derived.fillna(False).astype(bool)
assert (z.loc[z.included_derived,'A0_source']==z.loc[z.included_derived,'A0']).all()
assert (z.loc[z.included_derived,'A1_source']==z.loc[z.included_derived,'A1']).all()
orig=assess(d);assert orig['n_eligible_blocks']==31 and orig['n_eligible_donors']==12
blockkeys={(v['donor'],v['experiment'],float(v['forskolin_concentration_µM'])) for v in orig['blocks']}
z=z[z[['donor','experiment','forskolin_concentration_µM']].apply(tuple,axis=1).isin(blockkeys)].copy()
z=z[z.condition.isin(['DMSO','VX445_VX661_VX770'])];z=z[np.isfinite(z.A0_source)&(z.A0_source>0)].copy()
z['size_bin']=np.select([z.A0_source<722,z.A0_source>=1400],['small','large'],default='middle')
z['both_finite']=np.isfinite(z.A1_source)&(z.A1_source>0)
summary=[]
for (arm,size),g in z.groupby(['condition','size_bin']):
 both=g[g.both_finite]
 summary.append({'arm':arm,'size_bin':size,'n_baseline_detected':len(g),'n_both_positive_area':len(both),
  'n_derived':int(g.included_derived.sum()),'fraction_both_positive_area':float(g.both_finite.mean()),
  'fraction_derived_given_both':float(both.included_derived.mean()),
  'baseline_score_median_derived':float(g.loc[g.included_derived,'score0'].median()),
  'baseline_score_median_excluded':float(g.loc[~g.included_derived,'score0'].median()),
  'endpoint_score_median_derived':float(g.loc[g.included_derived,'score1'].median()),
  'endpoint_score_median_excluded':float(g.loc[~g.included_derived,'score1'].median())})
q=z[z.both_finite & z.size_bin.isin(['small','large'])].copy();q['log_swelling']=np.log(q.A1_source/q.A0_source)
group=['donor','experiment','forskolin_concentration_µM','condition','size_bin']
counts=q.groupby(group).agg(n_source=('log_swelling','size'),n_included=('included_derived','sum')).reset_index()
q=q.merge(counts,on=group,validate='many_to_one');q['weight']=(q.n_source/q.n_included).clip(upper=10).where(q.included_derived,0)
# Weighted diagnostic is on derived subset; constant cell weights cancel at cell-mean level, recorded honestly.
def effects(data,weighted=False):
 out=[]
 for key,g in data.groupby(group[:3]):
  cells={}
  for arm in ['DMSO','VX445_VX661_VX770']:
   for size in ['small','large']:
    h=g[(g.condition==arm)&(g.size_bin==size)]
    if len(h)<3:break
    cells[arm,size]=float(np.average(h.log_swelling,weights=h.weight)) if weighted else float(h.log_swelling.mean())
  if len(cells)!=4:continue
  e=(cells['VX445_VX661_VX770','large']-cells['DMSO','large'])-(cells['VX445_VX661_VX770','small']-cells['DMSO','small'])
  out.append({'donor':key[0],'experiment':key[1],'dose':key[2],'effect':e})
 v=pd.DataFrame(out); med=v.groupby('donor').effect.median()
 return {'n_blocks':len(v),'n_donors':len(med),'equal_donor_mean_median':float(med.mean()),'positive_donors':int((med>0).sum()),
         'donor_medians':{k:float(h) for k,h in med.items()},'block_effects':out}
raw=effects(q);derived_effect=effects(q[q.included_derived]);weighted=effects(q[q.included_derived],True)
assert derived_effect['n_blocks']==31 and abs(derived_effect['equal_donor_mean_median']-np.mean([r['median_effect'] for r in orig['per_donor'] if r['n_blocks']]))<1e-10
assert weighted['block_effects']==derived_effect['block_effects'] or all(abs(v['effect']-w['effect'])<1e-12 for v,w in zip(weighted['block_effects'],derived_effect['block_effects']))
out={'status':'post-result same-accession tabular track-retention sensitivity, not independent validation',
 'source_url':'https://zenodo.org/records/10610438','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
 'derived_sha256':hashlib.sha256(derived.read_bytes()).hexdigest(),'protocol':'notes/postresult_track_retention_protocol.md',
 'n_source_paired_particles_all_arms_doses':len(a),'n_derived_particles_all_arms_doses':len(d),
 'n_original_eligible_blocks':len(blockkeys),'retention_by_arm_size':summary,
 'source_paired_positive_area':raw,'derived_original':derived_effect,'cell_constant_inverse_inclusion_weight':weighted,
 'limits':'Source only contains detected/tracked particle rows, not all image objects. Blank endpoints and selection after detection may vary by swelling. Original image and mask paths are unavailable. Within-cell constant inverse inclusion weights exactly cancel in a within-cell mean, so they do not identify retention bias. No donor-independent cohort or preregistered biological confirmation.'}
Path('results/tabular_track_retention.json').write_text(json.dumps(out,indent=2)+'\n')
print('source:',raw['n_blocks'],raw['n_donors'],raw['equal_donor_mean_median'],raw['positive_donors'])
print('derived:',derived_effect['n_blocks'],derived_effect['n_donors'],derived_effect['equal_donor_mean_median'],derived_effect['positive_donors'])
print('summary',json.dumps(summary,indent=2))
