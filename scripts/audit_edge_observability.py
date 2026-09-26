"""Post-result image-boundary proxy check for positive endpoint observation."""
import hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
sys.path.insert(0,'src')
from blocked_size import assess
source=Path('data/raw/orgasegment/ResearchData/DIS_database.csv');derived=Path('data/raw/orgasegment/dis_merged_A0.csv')
x=pd.read_csv(source);d=pd.read_csv(derived);keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
a=x[x.t==0][keys+['size','score','x1','y1','x2','y2']].rename(columns={'size':'A0','score':'score0'})
b=x[x.t==1][keys+['size']].rename(columns={'size':'A1'})
assert len(a)==len(b)==47687 and not a.duplicated(keys).any() and not b.duplicated(keys).any()
z=a.merge(b,on=keys,validate='one_to_one');base=assess(d);assert base['n_eligible_blocks']==31
fixed={(q['donor'],q['experiment'],float(q['forskolin_concentration_µM'])) for q in base['blocks']}
z=z[z[['donor','experiment','forskolin_concentration_µM']].apply(tuple,axis=1).isin(fixed)]
z=z[z.condition.isin(['DMSO','VX445_VX661_VX770'])&(z.A0>0)&np.isfinite(z.A0)].copy()
z['size_bin']=np.select([z.A0<722,z.A0>=1400],['small','large'],default='middle');z=z[z.size_bin!='middle'].copy()
z['edge_clearance']=np.minimum.reduce([z.x1,z.y1,2048-z.x2,2048-z.y2]);assert np.isfinite(z.edge_clearance).all() and (z.edge_clearance>=0).all()
z['edge50']=(z.edge_clearance<50).astype(int);z['edge100']=(z.edge_clearance<100).astype(int)
z['observable']=((z.A1>0)&np.isfinite(z.A1)).astype(int);z['large']=(z.size_bin=='large').astype(int)
z['treatment']=(z.condition=='VX445_VX661_VX770').astype(int);z['log_A0']=np.log1p(z.A0)
assert z.score0.notna().all()
def model(edge):
 f=f'observable ~ treatment * large + {edge} + {edge}:large + {edge}:treatment + log_A0 + score0 + C(donor)'
 m=smf.glm(f,data=z,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':z.donor})
 terms=[edge,f'{edge}:large',f'{edge}:treatment','treatment:large']
 return {'nobs':int(m.nobs),'converged':bool(m.converged),'deviance':float(m.deviance),'terms':{t:{'logit_coef':float(m.params[t]),'odds_ratio':float(np.exp(m.params[t])),'cluster_se':float(m.bse[t]),'cluster_p_two_sided':float(m.pvalues[t]),'ci95_logit':[float(v) for v in m.conf_int().loc[t]]} for t in terms}}
counts=[]
for (arm,size,edge),v in z.groupby(['condition','size_bin','edge50']):counts.append({'arm':arm,'size_bin':size,'edge50':int(edge),'n':len(v),'seen':int(v.observable.sum()),'rate':float(v.observable.mean())})
block=[]
for b in base['blocks']:
 g=z[(z.donor==b['donor'])&(z.experiment==b['experiment'])&(z['forskolin_concentration_µM']==float(b['forskolin_concentration_µM']))]
 cells={}
 for arm in ['DMSO','VX445_VX661_VX770']:
  for size in ['small','large']:
   for edge in [0,1]:
    h=g[(g.condition==arm)&(g.size_bin==size)&(g.edge50==edge)]
    cells[f'{arm}:{size}:{edge}']={'n':len(h),'seen':int(h.observable.sum()),'rate':float(h.observable.mean()) if len(h) else None}
 eligible=all(v['n']>0 for v in cells.values())
 block.append({'donor':b['donor'],'experiment':b['experiment'],'dose':float(b['forskolin_concentration_µM']),'eligible_all_edge_cells':eligible,'cells':cells,
               'edge_minus_interior_four_cell':{f'{arm}:{size}':cells[f'{arm}:{size}:1']['rate']-cells[f'{arm}:{size}:0']['rate'] for arm in ['DMSO','VX445_VX661_VX770'] for size in ['small','large']} if eligible else None})
eligible=[v for v in block if v['eligible_all_edge_cells']]
agg={}
for arm in ['DMSO','VX445_VX661_VX770']:
 for size in ['small','large']:
  k=f'{arm}:{size}'
  if eligible:
   dmed=pd.DataFrame([{'donor':v['donor'],'diff':v['edge_minus_interior_four_cell'][k]} for v in eligible]).groupby('donor')['diff'].median()
   agg[k]={'equal_donor_mean_of_median_block_edge_minus_interior':float(dmed.mean()),'n_donors':len(dmed),'donor_medians':{s:float(t) for s,t in dmed.items()}}
  else:agg[k]=None
out={'status':'post-result image-boundary proxy audit; no missing outcome identification','protocol':'notes/postresult_edge_observability_protocol.md','source_url':'https://zenodo.org/records/10610438','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'derived_sha256':hashlib.sha256(derived.read_bytes()).hexdigest(),'n_baseline_detected':len(z),'n_observable':int(z.observable.sum()),'n_blocks':len(block),'edge50_model':model('edge50'),'edge100_model':model('edge100'),'edge50_counts':counts,'n_blocks_with_all_eight_edge_cells':len(eligible),'excluded_blocks':[{'donor':v['donor'],'experiment':v['experiment'],'dose':v['dose']} for v in block if not v['eligible_all_edge_cells']],'equal_donor_edge50_rate_differences':agg,'blocks':block,'limits':'Baseline box edge clearance is a crop-risk proxy, not actual truncation. Conditional associations do not identify absent A1, image masks or never-detected baseline objects; 12 donor clusters and fixed prior blocks. No original biological gate repair.'}
Path('results/edge_observability.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print({k:out[k] for k in ['n_baseline_detected','n_observable','n_blocks','n_blocks_with_all_eight_edge_cells','edge50_counts','equal_donor_edge50_rate_differences']})
print('model50',out['edge50_model']['terms']);print('model100',out['edge100_model']['terms'])
