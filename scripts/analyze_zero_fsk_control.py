"""Frozen 0 vs .128 uM forskolin within-donor/experiment specificity contrast."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest

source=Path('data/raw/orgasegment/dis_merged_A0.csv')
df=pd.read_csv(source)
assert df.donor.nunique()==17
D='VX445_VX661_VX770'; C='DMSO'
z=df[df.condition.isin([D,C]) & df['forskolin_concentration_µM'].isin([0.0,.128])].copy()
z=z[np.isfinite(z.A0)&np.isfinite(z.swelling)&(z.A0>0)&(z.swelling>0)].copy()
z['bin']=np.select([z.A0<722,z.A0>=1400],['small','large'],default='middle')
z=z[z.bin!='middle'].copy();z['logsw']=np.log(z.swelling)
cells=z.groupby(['donor','experiment','forskolin_concentration_µM','condition','bin']).logsw.agg(['mean','size'])
records=[];all_exps=z[['donor','experiment']].drop_duplicates()
for donor,experiment in all_exps.itertuples(index=False,name=None):
    effects={};min_counts={}
    for dose in (0.,.128):
        cc={(arm,bin):cells.loc[(donor,experiment,dose,arm,bin)]
            for arm in (D,C) for bin in ('small','large')
            if (donor,experiment,dose,arm,bin) in cells.index}
        if len(cc)!=4 or any(c['size']<3 for c in cc.values()):continue
        effects[str(dose)]=float((cc[(D,'large')]['mean']-cc[(C,'large')]['mean'])-
                                  (cc[(D,'small')]['mean']-cc[(C,'small')]['mean']))
        min_counts[str(dose)]=int(min(c['size'] for c in cc.values()))
    if len(effects)!=2:continue
    records.append({'donor':str(donor),'experiment':str(experiment),'effect_0':effects['0.0'],
                    'effect_0128':effects['0.128'],'delta':effects['0.128']-effects['0.0'],
                    'min_cell_0':min_counts['0.0'],'min_cell_0128':min_counts['0.128']})
matched=pd.DataFrame(records)
donors=[]
for donor,group in matched.groupby('donor') if len(matched) else []:
    donors.append({'donor':donor,'paired_experiments':len(group),
                   'median_effect_0':float(group.effect_0.median()),
                   'median_effect_0128':float(group.effect_0128.median()),
                   'median_delta':float(group.delta.median())})
n=len(donors);pos=sum(d['median_delta']>0 for d in donors)
p=float(binomtest(pos,n,.5,alternative='greater').pvalue) if n else None
out={'source':'https://zenodo.org/records/10610438',
     'preregistration':'results/preregistration_zero_fsk_control.md',
     'donors_source':17,'donor_experiments_source':len(all_exps),
     'paired_experiments':len(records),'eligible_donors':n,'excluded_donors':17-n,
     'positive_paired_delta_donors':pos,'exact_one_sided_sign_p':p,
     'verdict':'UNINTERPRETABLE' if n<8 else 'PASS' if pos>n/2 and p<.05 else 'FAIL',
     'warning':'same accession, post hoc source-derived size cutoffs; not independent validation',
     'donors':donors,'experiments':records}
Path('results/zero_fsk_control.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('donor_experiments_source','paired_experiments','eligible_donors','excluded_donors','positive_paired_delta_donors','exact_one_sided_sign_p','verdict')},indent=2))
print('donors:',[(x['donor'],round(x['median_effect_0'],3),round(x['median_effect_0128'],3),round(x['median_delta'],3)) for x in donors])
