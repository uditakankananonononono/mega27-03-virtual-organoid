"""Post-result association of baseline neighbor spacing and positive endpoints."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
import statsmodels.api as sm
import statsmodels.formula.api as smf
src=Path('data/raw/orgasegment/ResearchData/DIS_database.csv');sha=hashlib.sha256(src.read_bytes()).hexdigest()
assert sha=='0e0e7913947a29c7b32ce87ca73a98614033d41ce63063a1d65053e766b2f315'
raw=pd.read_csv(src)
keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
a=raw[raw.t==0][keys+['size','score','x','y','image']].rename(columns={'size':'A0','score':'score0','x':'x0','y':'y0','image':'image0'})
b=raw[raw.t==1][keys+['size']].rename(columns={'size':'A1'})
assert len(a)==len(b)==47687 and not a.duplicated(keys).any() and not b.duplicated(keys).any()
z=a.merge(b,on=keys,validate='one_to_one');assert len(z)==47687
blocks=json.load(open('results/blocked_size.json'))['blocks'];assert len(blocks)==31
fixed={(q['donor'],q['experiment'],float(q['forskolin_concentration_µM'])) for q in blocks}
zero=z[z['forskolin_concentration_µM']==0]
zero_matched=z[(z['forskolin_concentration_µM']==0)&z[['donor','experiment']].apply(tuple,axis=1).isin({(t[0],t[1]) for t in fixed})]
z=z[z[['donor','experiment','forskolin_concentration_µM']].apply(tuple,axis=1).isin(fixed)].copy()
z=z[z.condition.isin(['DMSO','VX445_VX661_VX770'])&(z.A0>0)&np.isfinite(z.A0)].copy()
z['size_bin']=np.select([z.A0<722,z.A0>=1400],['small','large'],default='middle');z=z[z.size_bin!='middle'].copy()
assert len(z)==6964 and z[['x0','y0','score0']].notna().all().all()
z['observable']=((z.A1>0)&np.isfinite(z.A1)).astype(int)
assert int(z.observable.sum())==4964
z['nn1']=np.nan;z['nn3']=np.nan;z['field_n']=0
# Neighbors from the whole positive-A0 source baseline field, including middle size and other arms as applicable.
fields=a[(a.A0>0)&np.isfinite(a.A0)&a.x0.notna()&a.y0.notna()].groupby('image0')
lookup={}
for image,g in fields:
 coords=g[['x0','y0']].to_numpy(float);k=min(len(g),4)
 dist,_=cKDTree(coords).query(coords,k=k)
 if dist.ndim==1:dist=dist[:,None]
 for row,dd in zip(g.itertuples(index=False),dist):
  key=tuple(getattr(row, v) for v in keys)
  assert key not in lookup
  lookup[key]=(float(dd[1]) if k>=2 else np.nan,float(dd[3]) if k>=4 else np.nan,len(g))
for i,row in z.iterrows():
 d=lookup[tuple(row[v] for v in keys)];z.at[i,'nn1']=d[0];z.at[i,'nn3']=d[1];z.at[i,'field_n']=d[2]
assert z.nn1.notna().all();z['log_nn1']=np.log1p(z.nn1);z['large']=(z.size_bin=='large').astype(int);z['treatment']=(z.condition=='VX445_VX661_VX770').astype(int)
z['log_A0']=np.log1p(z.A0)
# This model is post-result descriptive: donor clusters few, field dependence and cell-level selection.
form='observable ~ treatment*large + log_A0 + score0 + log_nn1*large + C(donor)'
fit=smf.glm(formula=form,data=z,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':z.donor})
terms=['treatment','large','treatment:large','log_A0','score0','log_nn1','log_nn1:large']
def summary(m):return {v:{'beta':float(m.params[v]),'odds_ratio':float(np.exp(m.params[v])),
 'cluster_se':float(m.bse[v]),'cluster_p_two_sided':float(m.pvalues[v]),
 'cluster_ci95_beta':[float(t) for t in m.conf_int().loc[v]]} for v in terms}
cut=np.quantile(z.nn1,[.25,.5,.75]);z['neighbor_quartile']=np.searchsorted(cut,z.nn1,side='right')+1
rates=[]
for (arm,size,q),g in z.groupby(['condition','size_bin','neighbor_quartile']):rates.append({'arm':arm,'size':size,'quartile':int(q),'n':len(g),'observed':int(g.observable.sum()),'rate':float(g.observable.mean()),'median_nn1':float(g.nn1.median())})
# Within-field close-minus-far rates, bounded to fields with both halves; 100 shuffles of exposure within field.
z['close']=(z.nn1<=z.groupby('image0').nn1.transform('median')).astype(int)
usable=[g for _,g in z.groupby('image0') if len(g)>=4 and g.close.nunique()==2]
actual=np.mean([g[g.close==1].observable.mean()-g[g.close==0].observable.mean() for g in usable])
rng=np.random.default_rng(20260927); null=[]
for i in range(100):
 diffs=[]
 for g in usable:
  shuffled=rng.permutation(g.close.to_numpy());y=g.observable.to_numpy();diffs.append(float(y[shuffled==1].mean()-y[shuffled==0].mean()))
 null.append(float(np.mean(diffs)))
# Return no causal inference or significance test, just empirical scale vs technical null.
out={'status':'post-result source-baseline neighbor-spacing vs endpoint observability, same accession and donor set',
 'protocol':'notes/postresult_neighbor_observability_protocol.md','source_url':'https://zenodo.org/records/10610438','source_sha256':sha,
 'n_particles_source':len(a),'n_source_images':a.image0.nunique(),'n_original_blocks':len(blocks),'n_eligible_baseline':len(z),'n_observed_endpoint':int(z.observable.sum()),
 'n_donors':z.donor.nunique(),'n_eligible_fields':z.image0.nunique(),'n_missing_nn1':int(z.nn1.isna().sum()),'n_missing_nn3':int(z.nn3.isna().sum()),
 'n_exact_coincident_neighbor':int((z.nn1==0).sum()),'nn1_quartile_cutoffs_px':[float(q) for q in cut],
 'quartile_rates':rates,'model_terms':summary(fit),'model_converged':bool(fit.converged),
 'within_field_close_minus_far':{'n_usable_fields':len(usable),'observed_mean_field_difference':float(actual),
  'within_field_permutation_seed':20260927,'n_permutations':100,'null_mean':float(np.mean(null)),
  'null_percentile_2_5_97_5':[float(q) for q in np.quantile(null,[.025,.975])]},
 'zero_fsk_source_rows':len(zero),'zero_fsk_rows_same_donor_experiment_as_original_blocks':len(zero_matched),
 'limits':'Observability means positive finite source A1, not t1 row presence. Donor fixed-effects descriptive association with 12 clusters and within-field dependence. Zero FSK no clean mechanical control; no missing A1 identification, causality, original sign-gate repair or biology.'}
Path('results/postresult_neighbor_observability.json').write_text(json.dumps(out,indent=2)+'\n')
print('n',len(z),'observed',z.observable.sum(),'fields',z.image0.nunique(),'NN quartiles',cut)
print('terms',json.dumps({k:v for k,v in out['model_terms'].items() if 'nn' in k},indent=2))
print('within-field',out['within_field_close_minus_far'])
