"""Locked synthetic-mask method benchmark; no real missing outcome identification."""
import hashlib,json,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy.special import expit

SOURCE=Path('data/raw/orgasegment/ResearchData/DIS_database.csv')
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()=='0e0e7913947a29c7b32ce87ca73a98614033d41ce63063a1d65053e766b2f315'
raw=pd.read_csv(SOURCE)
keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
a=raw.loc[raw.t==0,keys+['size','score']].rename(columns={'size':'A0','score':'score0'})
b=raw.loc[raw.t==1,keys+['size']].rename(columns={'size':'A1'})
z=a.merge(b,on=keys,validate='one_to_one')
blocks=json.load(open('results/blocked_size.json'))['blocks'];assert len(blocks)==31
fixed={(r['donor'],r['experiment'],float(r['forskolin_concentration_µM'])) for r in blocks}
z=z[z[['donor','experiment','forskolin_concentration_µM']].apply(tuple,axis=1).isin(fixed)]
z=z[z.condition.isin(['DMSO','VX445_VX661_VX770'])&(z.A0>0)&(z.A1>0)&np.isfinite(z.A0)&np.isfinite(z.A1)].copy()
z['size_bin']=np.select([z.A0<722,z.A0>=1400],['small','large'],default='middle');z=z[z.size_bin!='middle'].copy().reset_index(drop=True)
z['log_response']=np.log(z.A1/z.A0)
assert len(z)>4000
sign={'VX445_VX661_VX770:large':1,'DMSO:large':-1,'VX445_VX661_VX770:small':-1,'DMSO:small':1}
z['cell']=z.condition+':'+z.size_bin
index={(r['donor'],r['experiment'],float(r['forskolin_concentration_µM']),cell):v.index.to_numpy()
       for (d,e,dose,cell),v in z.groupby(['donor','experiment','forskolin_concentration_µM','cell'])
       for r in [{'donor':d,'experiment':e,'forskolin_concentration_µM':dose}]}
for block in blocks:
 for cell in sign: assert len(index[(block['donor'],block['experiment'],float(block['forskolin_concentration_µM']),cell)])>=3
values=z.log_response.to_numpy();bycell={c:z.index[z.cell==c].to_numpy() for c in sign}
loga=np.log1p(z.A0.to_numpy());sc=z.score0.to_numpy()
def standardize(v):
 return (v-v.mean())/v.std(ddof=0)
base_logit=-1.5+.45*(z.size_bin.to_numpy()=='small')+.35*(z.condition.to_numpy()=='VX445_VX661_VX770')+.25*standardize(loga)-.25*standardize(sc)
mar=np.clip(expit(base_logit),.05,.60)
mechanisms={'MCAR':np.full(len(z),.20),'MAR':mar,'MNAR':np.clip(expit(base_logit+.8*standardize(values)),.05,.60)}

def aggregate(cellmean):
 per_donor={}
 for block in blocks:
  d,e,dose=block['donor'],block['experiment'],float(block['forskolin_concentration_µM'])
  effect=sum(s*cellmean[(d,e,dose,c)] for c,s in sign.items())
  per_donor.setdefault(d,[]).append(effect)
 return float(np.mean([np.median(v) for v in per_donor.values()])),int(sum(np.median(v)>0 for v in per_donor.values()))

def means_from(visible, weights=None, lo=None, hi=None):
 vals={}
 for block in blocks:
  d,e,dose=block['donor'],block['experiment'],float(block['forskolin_concentration_µM'])
  for cell in sign:
   ids=index[(d,e,dose,cell)];have=ids[visible[ids]];missing=ids[~visible[ids]]
   if lo is None:
    w=None if weights is None else weights[have]
    vals[(d,e,dose,cell)]=float(np.average(values[have],weights=w))
   else:
    bound=lo[cell] if sign[cell]>0 else hi[cell]
    vals[(d,e,dose,cell)]=float((values[have].sum()+len(missing)*bound)/len(ids))
 return aggregate(vals)

truth=means_from(np.ones(len(z),dtype=bool))[0]
summary={}
for mech,p_miss in mechanisms.items():
 records=[];hidden_outside=[];missing_rates={c:[] for c in sign};failed=0
 for seed in range(20260928,20261028):
  rng=np.random.default_rng(seed);visible=rng.random(len(z))>=p_miss
  for block in blocks:
   d,e,dose=block['donor'],block['experiment'],float(block['forskolin_concentration_µM'])
   for cell in sign:
    ids=index[(d,e,dose,cell)]; tries=0
    while visible[ids].sum()<3 and tries<100:
     visible[ids]=rng.random(len(ids))>=p_miss[ids];tries+=1
  if any(visible[ids].sum()<3 for ids in index.values()):failed+=1;continue
  cc=means_from(visible)[0]
  proxy=np.full(len(z),.20) if mech=='MCAR' else mar
  weights=1/np.maximum(1-proxy,.05)
  ipw=means_from(visible,weights=np.minimum(weights,20))[0]
  selected=values[visible]
  g05,g95=np.quantile(selected,[.05,.95]);local_lo={};local_hi={};outside=0;hidden=0
  for cell,ids in bycell.items():
   local_lo[cell],local_hi[cell]=np.quantile(values[ids[visible[ids]]],[.05,.95])
   hidden_ids=ids[~visible[ids]];hidden+=len(hidden_ids)
   outside+=int(((values[hidden_ids]<local_lo[cell])|(values[hidden_ids]>local_hi[cell])).sum())
   missing_rates[cell].append(float(len(hidden_ids)/len(ids)))
  lower=means_from(visible,lo=local_lo,hi=local_hi)[0]
  upper=means_from(visible,lo=local_hi,hi=local_lo)[0]
  glo={c:g05 for c in sign};ghi={c:g95 for c in sign}
  glower=means_from(visible,lo=glo,hi=ghi)[0]
  gupper=means_from(visible,lo=ghi,hi=glo)[0]
  assert lower<=upper and glower<=gupper
  records.append({'seed':seed,'complete_case':cc,'ipw_known_baseline_propensity':ipw,
                  'cell_interval':[lower,upper],'global_interval':[glower,gupper],
                  'hidden_outside_cell_percentile':outside,'n_hidden':hidden})
 summary[mech]={'eligible_masks':len(records),'ineligible_masks':failed,'missing_fraction_by_cell':{c:float(np.mean(r)) for c,r in missing_rates.items()},
  'complete_case_bias':float(np.mean([r['complete_case']-truth for r in records])),
  'complete_case_mae':float(np.mean([abs(r['complete_case']-truth) for r in records])),
  'ipw_bias':float(np.mean([r['ipw_known_baseline_propensity']-truth for r in records])),
  'ipw_mae':float(np.mean([abs(r['ipw_known_baseline_propensity']-truth) for r in records])),
  'hidden_outside_cell_percentile_fraction':float(sum(r['hidden_outside_cell_percentile'] for r in records)/sum(r['n_hidden'] for r in records)),
  'intervals':{key:{'containment':float(np.mean([r[key][0]<=truth<=r[key][1] for r in records])),
                     'median_width':float(np.median([r[key][1]-r[key][0] for r in records])),
                     'mean_width':float(np.mean([r[key][1]-r[key][0] for r in records])),
                     'zero_contained_fraction':float(np.mean([r[key][0]<=0<=r[key][1] for r in records]))}
               for key in ('cell_interval','global_interval')},
  'replicates':records}
 summary[mech]['median_cell_to_global_width_ratio']=summary[mech]['intervals']['cell_interval']['median_width']/summary[mech]['intervals']['global_interval']['median_width']
pass_gate=all(x['intervals']['cell_interval']['containment']>=.90 and x['median_cell_to_global_width_ratio']<=.90 for x in summary.values())
out={'status':'same-accession synthetic observed-track masking; not real-missing endpoint identification','protocol':'notes/prereg_masked_observability_method_benchmark_20260928.md','source_url':'https://zenodo.org/records/10610438','source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
 'n_full_observed_tracks':len(z),'n_fixed_blocks':len(blocks),'n_fixed_donors':len({b['donor'] for b in blocks}),'full_known_equal_donor_mean_median':truth,
 'mechanisms':summary,'locked_method_gate':'PASS' if pass_gate else 'FAIL',
 'limitations':'Same-source post-result simulated masks over preselected positive tracks; true missing endpoints and never-detected objects remain unknown. 5/95 percentiles are scenarios not support bounds; IPW knows simulator propensity and MNAR is misspecified. Original biological/segmentation gates unchanged.'}
Path('results/masked_observability_method_benchmark.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'truth':truth,'n':len(z),'gate':out['locked_method_gate'],'summary':{k:{kk:vv for kk,vv in x.items() if kk!='replicates'} for k,x in summary.items()}},indent=2))
