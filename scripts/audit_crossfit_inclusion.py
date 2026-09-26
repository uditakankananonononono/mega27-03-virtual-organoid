"""Post-result cross-fitted positive-pair derived-inclusion sensitivity, not a new gate."""
import hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss
sys.path.insert(0,'src')
from blocked_size import assess
source=Path('data/raw/orgasegment/ResearchData/DIS_database.csv');derived=Path('data/raw/orgasegment/dis_merged_A0.csv')
x=pd.read_csv(source);d=pd.read_csv(derived);keys=['donor','experiment','well','condition','forskolin_concentration_µM','particle']
a=x[x.t==0][keys+['size','score']].rename(columns={'size':'A0','score':'score0'});b=x[x.t==1][keys+['size','score']].rename(columns={'size':'A1','score':'score1'})
assert len(a)==len(b)==47687
z=a.merge(b,on=keys,validate='one_to_one');r=d[keys].copy();r['included']=1
z=z.merge(r,on=keys,how='left',validate='one_to_one');z['included']=z.included.fillna(0).astype(int)
blocks=assess(d)['blocks'];fixed={(v['donor'],v['experiment'],float(v['forskolin_concentration_µM'])) for v in blocks}
z=z[z[['donor','experiment','forskolin_concentration_µM']].apply(tuple,axis=1).isin(fixed)].copy()
z=z[z.condition.isin(['DMSO','VX445_VX661_VX770']) & (z.A0>0) & (z.A1>0) & np.isfinite(z.A0) & np.isfinite(z.A1)].copy()
z['size_bin']=np.select([z.A0<722,z.A0>=1400],['small','large'],default='middle');z=z[z.size_bin!='middle'].copy()
z['log_A0']=np.log1p(z.A0);z['log_A1']=np.log1p(z.A1);z['log_swelling']=np.log(z.A1/z.A0)
num=['log_A0','log_A1','score0','score1','forskolin_concentration_µM'];cat=['donor','condition','well']
cv=GroupKFold(n_splits=5);prob=np.zeros(len(z));fold=[]
for train,test in cv.split(z,z.included,z.experiment):
 trans=ColumnTransformer([('numeric',make_pipeline(SimpleImputer(strategy='median'),StandardScaler()),num),
                          ('category',OneHotEncoder(handle_unknown='ignore'),cat)])
 pipe=make_pipeline(trans,LogisticRegression(C=1,max_iter=1000,class_weight=None))
 pipe.fit(z.iloc[train][num+cat],z.iloc[train].included)
 prob[test]=pipe.predict_proba(z.iloc[test][num+cat])[:,1]
 fold.append({'n_train':len(train),'n_test':len(test),'n_test_included':int(z.iloc[test].included.sum()),
              'test_experiments':sorted(z.iloc[test].experiment.unique().tolist())})
clipped=np.clip(prob,.05,.995);z['weight']=1/clipped;z['pred']=prob
assert len(z[z.included==1])>0
# Preserve ordered original blocks and all 31 outcomes.
def calculate(df,weighted=False):
 rows=[]
 for b in blocks:
  g=df[(df.donor==b['donor'])&(df.experiment==b['experiment'])&(df['forskolin_concentration_µM']==float(b['forskolin_concentration_µM']))]
  means={}
  for arm in ['DMSO','VX445_VX661_VX770']:
   for size in ['small','large']:
    h=g[(g.condition==arm)&(g.size_bin==size)]
    if len(h)<3:raise RuntimeError('original block cell became ineligible')
    means[arm,size]=float(np.average(h.log_swelling,weights=h.weight)) if weighted else float(h.log_swelling.mean())
  e=(means['VX445_VX661_VX770','large']-means['DMSO','large'])-(means['VX445_VX661_VX770','small']-means['DMSO','small'])
  rows.append({'donor':b['donor'],'experiment':b['experiment'],'dose':float(b['forskolin_concentration_µM']),'effect':e})
 med=pd.DataFrame(rows).groupby('donor').effect.median()
 return {'n_blocks':len(rows),'n_donors':len(med),'positive_donors':int((med>0).sum()),
         'equal_donor_mean_median':float(med.mean()),'donor_medians':{k:float(v) for k,v in med.items()},'block_effects':rows}
full=calculate(z);unweighted=calculate(z[z.included==1]);weighted=calculate(z[z.included==1],True)
assert len(blocks)==31 and len(unweighted['donor_medians'])==12
prior=json.load(open('results/tabular_track_retention.json'))
assert abs(unweighted['equal_donor_mean_median']-prior['derived_original']['equal_donor_mean_median'])<1e-10
assert abs(full['equal_donor_mean_median']-prior['source_paired_positive_area']['equal_donor_mean_median'])<1e-10
rng=np.random.default_rng(20260926);neg=[];inc=z[z.included==1]
for rep in range(300):
 # DMSO-only split into pseudo arms within donor/experiment/dose/size; same split in each group without outcome inspection.
 dm=inc[inc.condition=='DMSO'].copy();dm['pseudo']='DMSO'
 group=['donor','experiment','forskolin_concentration_µM','size_bin']
 for _,ids in dm.groupby(group).groups.items():
  ids=np.array(list(ids));sample=rng.choice(ids,size=len(ids)//2,replace=False);dm.loc[sample,'pseudo']='FAKE'
 rows=[]
 for b in blocks:
  g=dm[(dm.donor==b['donor'])&(dm.experiment==b['experiment'])&(dm['forskolin_concentration_µM']==float(b['forskolin_concentration_µM']))]
  vals={}
  for arm in ['DMSO','FAKE']:
   for size in ['small','large']:
    h=g[(g.pseudo==arm)&(g.size_bin==size)]
    if len(h)==0:raise RuntimeError('pseudo cell empty')
    vals[arm,size]=float(np.average(h.log_swelling,weights=h.weight))
  rows.append((b['donor'],(vals['FAKE','large']-vals['DMSO','large'])-(vals['FAKE','small']-vals['DMSO','small'])))
 med=pd.DataFrame(rows,columns=['donor','effect']).groupby('donor').effect.median()
 neg.append({'mean_donor_median':float(med.mean()),'positive_donors':int((med>0).sum())})
ws=inc.weight.to_numpy();ess=float(ws.sum()**2/(ws**2).sum())
out={'status':'post-result same-accession observed-positive-pair inclusion sensitivity','protocol':'notes/postresult_inclusion_model_protocol.md',
 'source_url':'https://zenodo.org/records/10610438','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
 'derived_sha256':hashlib.sha256(derived.read_bytes()).hexdigest(),'n_positive_source_pairs_in_original_small_large_cells':len(z),
 'n_derived_included_in_frame':int(z.included.sum()),'folds':fold,'out_of_fold_inclusion_auc':float(roc_auc_score(z.included,prob)),
 'out_of_fold_inclusion_ap':float(average_precision_score(z.included,prob)),
 'out_of_fold_inclusion_brier':float(brier_score_loss(z.included,prob)),
 'constant_prevalence_brier':float(brier_score_loss(z.included,np.full(len(z),z.included.mean()))),
 'raw_probability_range':[float(prob.min()),float(prob.max())],
 'n_probabilities_clipped':int((clipped!=prob).sum()),'inclusion_rate':float(z.included.mean()),
 'included_weight_quantiles':{str(q):float(v) for q,v in zip([0,.5,.9,.99,1],np.quantile(ws,[0,.5,.9,.99,1]))},
 'included_effective_sample_size':ess,'source_positive_pair':full,'derived_unweighted':unweighted,'derived_ipw':weighted,
 'delta_weighted_minus_derived':weighted['equal_donor_mean_median']-unweighted['equal_donor_mean_median'],
 'negative_control_seed':20260926,'negative_control_repetitions':len(neg),'negative_control_pseudo_effects':neg,
 'negative_control_abs_ge_real_fraction':float(np.mean([abs(r['mean_donor_median'])>=abs(weighted['equal_donor_mean_median']) for r in neg])),
 'limits':'Endpoint variables used in selection model; not pretreatment causal weighting. Positive-area source pairs only. No baseline-undetected objects, pixel masks or missing A1 outcomes. Post-result observational sensitivity cannot repair sign-only/specificity gates or independent validation.'}
Path('results/crossfit_inclusion_sensitivity.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['n_positive_source_pairs_in_original_small_large_cells','n_derived_included_in_frame','out_of_fold_inclusion_auc','out_of_fold_inclusion_brier','constant_prevalence_brier','raw_probability_range','n_probabilities_clipped','included_weight_quantiles','included_effective_sample_size','delta_weighted_minus_derived','negative_control_abs_ge_real_fraction']},indent=2))
for name in ['source_positive_pair','derived_unweighted','derived_ipw']:
 a=out[name];print(name,a['equal_donor_mean_median'],a['positive_donors'])
