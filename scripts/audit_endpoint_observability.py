"""Post-result t1 observability conditional on detected baseline objects."""
import hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
sys.path.insert(0,'src')
from blocked_size import assess
src=Path('data/raw/orgasegment/ResearchData/DIS_database.csv');d=pd.read_csv('data/raw/orgasegment/dis_merged_A0.csv')
x=pd.read_csv(src);keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
a=x[x.t==0][keys+['size','score']].rename(columns={'size':'A0','score':'score0'});b=x[x.t==1][keys+['size']].rename(columns={'size':'A1'})
assert len(a)==len(b)==47687;z=a.merge(b,on=keys,validate='one_to_one')
base=assess(d);assert base['n_eligible_blocks']==31
fixed={(q['donor'],q['experiment'],float(q['forskolin_concentration_µM'])) for q in base['blocks']}
z=z[z[['donor','experiment','forskolin_concentration_µM']].apply(tuple,axis=1).isin(fixed)].copy()
z=z[z.condition.isin(['DMSO','VX445_VX661_VX770'])&(z.A0>0)&np.isfinite(z.A0)].copy()
z['size_bin']=np.select([z.A0<722,z.A0>=1400],['small','large'],default='middle');z=z[z.size_bin!='middle'].copy()
z['large']=(z.size_bin=='large').astype(int);z['treatment']=(z.condition=='VX445_VX661_VX770').astype(int)
z['observable']=((z.A1>0)&np.isfinite(z.A1)).astype(int);z['log_A0']=np.log1p(z.A0)
z['high_score0']=(z.score0>=.99).astype(int)
z['score0_imputed']=z.score0.fillna(z.score0.median())
form='observable ~ treatment * large + log_A0 + score0_imputed + C(donor)'
fit=smf.glm(formula=form,data=z,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':z.donor})
control=smf.glm('high_score0 ~ treatment * large + log_A0 + C(donor)',data=z,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':z.donor})
count=z.groupby(['donor','experiment','forskolin_concentration_µM','condition','size_bin']).agg(n=('observable','size'),seen=('observable','sum'),n_high0=('high_score0','sum')).reset_index()
assert not count.n.eq(0).any()
blocks=[]
for b in base['blocks']:
 g=count[(count.donor==b['donor'])&(count.experiment==b['experiment'])&(count['forskolin_concentration_µM']==float(b['forskolin_concentration_µM']))]
 cells={}
 for arm in ['DMSO','VX445_VX661_VX770']:
  for size in ['small','large']:
   h=g[(g.condition==arm)&(g.size_bin==size)];assert len(h)==1 and int(h.n.iloc[0])>=3
   cells[arm+':'+size]={'n':int(h.n.iloc[0]),'seen':int(h.seen.iloc[0]),'rate':float(h.seen.iloc[0]/h.n.iloc[0]),
                         'baseline_high_score_rate':float(h.n_high0.iloc[0]/h.n.iloc[0])}
 diff_small=cells['VX445_VX661_VX770:small']['rate']-cells['DMSO:small']['rate']
 diff_large=cells['VX445_VX661_VX770:large']['rate']-cells['DMSO:large']['rate']
 blocks.append({'donor':b['donor'],'experiment':b['experiment'],'dose':float(b['forskolin_concentration_µM']),
                'cells':cells,'treatment_gap_small':diff_small,'treatment_gap_large':diff_large,
                'gap_interaction_large_minus_small':diff_large-diff_small})
frame=pd.DataFrame([{k:v for k,v in b.items() if k!='cells'} for b in blocks]);dmed=frame.groupby('donor').agg({c:'median' for c in ['treatment_gap_small','treatment_gap_large','gap_interaction_large_minus_small']})
def result(m):
 return {'terms':{term:{'coefficient':float(m.params[term]),'odds_ratio':float(np.exp(m.params[term])),
                 'cluster_se':float(m.bse[term]),'cluster_p_two_sided':float(m.pvalues[term]),
                 'ci95_logit':[float(v) for v in m.conf_int().loc[term]]} for term in ['treatment','large','treatment:large','log_A0']},
         'converged':bool(m.converged),'nobs':int(m.nobs),'deviance':float(m.deviance)}
out={'status':'post-result conditional endpoint-observability and t0 score control, same accession','protocol':'notes/postresult_endpoint_observability_protocol.md',
 'source_url':'https://zenodo.org/records/10610438','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
 'n_baseline_detected_original_extreme_size_cells':len(z),'n_observable_endpoint':int(z.observable.sum()),
 'n_original_blocks':len(blocks),'n_donors':len(dmed),'score0_missing':int(z.score0.isna().sum()),
 'endpoint_model':result(fit),'baseline_high_score_negative_control':result(control),
 'equal_donor_mean_median_gaps':{c:float(dmed[c].mean()) for c in dmed.columns},
 'donor_medians':{k:{c:float(v) for c,v in row.items()} for k,row in dmed.to_dict(orient='index').items()},
 'blocks':blocks,'limits':'Conditioned on source-detected t0 objects; no missing outcome identification or unseen baseline objects. Donor fixed effects and donor-cluster standard errors, not experiment fixed effects due to collinearity in original block subset; only 12 clusters and batch confounding. P-values descriptive and no original gate repair.'}
Path('results/endpoint_observability.json').write_text(json.dumps(out,indent=2)+'\n')
print('n',len(z),'seen',z.observable.sum(),'score0 missing',z.score0.isna().sum())
print('model',json.dumps(out['endpoint_model']['terms'],indent=2));print('control',json.dumps(out['baseline_high_score_negative_control']['terms'],indent=2));print('donor gaps',out['equal_donor_mean_median_gaps'])
