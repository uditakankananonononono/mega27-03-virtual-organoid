"""Internal experiment split of discovery accession (never independent cohort)."""
import json,re,pandas as pd,numpy as np
from pathlib import Path
from scipy.stats import binomtest
PATH='data/raw/orgasegment/dis_merged_A0.csv'

def date_key(s):
    m=re.search(r'20\d{6}',s)
    if not m:raise ValueError('Experiment date not found: '+s)
    return m.group(0)

def assess(x):
    y=x[(x['forskolin_concentration_µM']>0)& x.condition.isin(['DMSO','VX445_VX661_VX770']) & (x.swelling>0)].copy()
    y['date']=y.experiment.map(date_key)
    # Ties at the earliest date are all "early"; every later date is independent of early plates.
    earliest=y.groupby('donor').date.transform('min');y['portion']=np.where(y.date==earliest,'earliest','later')
    y['size']=np.select([y.A0<722,y.A0>=1400],['small','large'],default='middle')
    y=y[y['size']!='middle'].copy();y['log_swelling']=np.log(y.swelling)
    rows=[]
    for (portion,donor),z in y.groupby(['portion','donor']):
        q=z.groupby(['condition','size']).log_swelling.agg(['mean','size'])
        cells=[(c,s) for c in ['DMSO','VX445_VX661_VX770'] for s in ['small','large']]
        if any(i not in q.index or q.loc[i,'size']<5 for i in cells):
            rows.append({'portion':portion,'donor':donor,'valid':False,'cell_counts':{f'{c}:{s}':int(q.loc[(c,s),'size']) if (c,s) in q.index else 0 for c,s in cells}});continue
        effect=(q.loc[('VX445_VX661_VX770','large'),'mean']-q.loc[('DMSO','large'),'mean'])-(q.loc[('VX445_VX661_VX770','small'),'mean']-q.loc[('DMSO','small'),'mean'])
        rows.append({'portion':portion,'donor':donor,'valid':True,'attenuation':float(effect),'cell_counts':{f'{c}:{s}':int(q.loc[(c,s),'size']) for c,s in cells}})
    out={'preregistration':'results/preregistration_withincohort_replication.md (commit 8d80877)',
         'source':'OrgaSegment DIS, Zenodo 10610438, same accession as original discovery',
         'size_cutoffs_px':{'small_lt':722,'large_gte':1400},'rows':rows,'portions':{}}
    for portion in ['earliest','later']:
        v=[r['attenuation'] for r in rows if r['portion']==portion and r['valid']]
        pos=sum(a>0 for a in v);n=len(v);p=binomtest(pos,n,.5,alternative='greater').pvalue if n else None
        out['portions'][portion]={'n_donors':n,'n_positive':pos,'median_attenuation':float(np.median(v)) if n else None,'sign_p_one_sided':p}
    v=out['portions']['later'];out['H1']={'verdict':'UNINTERPRETABLE' if v['n_donors']<8 else 'PASS' if v['median_attenuation']>0 and v['sign_p_one_sided']<.05 else 'FAIL',
                                            'reason':'fewer than eight eligible later-experiment donors' if v['n_donors']<8 else None}
    return out

if __name__=='__main__':
    out=assess(pd.read_csv(PATH));Path('results/withincohort_replication.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['portions'],out['H1'])
