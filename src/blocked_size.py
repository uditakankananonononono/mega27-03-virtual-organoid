"""Same-accession plate and dose matched robustness test."""
import json,pandas as pd,numpy as np
from scipy.stats import binomtest
from pathlib import Path


def assess(x):
    y=x[(x['forskolin_concentration_µM']>0)&x.condition.isin(['DMSO','VX445_VX661_VX770'])&(x.swelling>0)].copy()
    y['size']=np.select([y.A0<722,y.A0>=1400],['small','large'],default='middle')
    y=y[y['size']!='middle'].copy();y['log_swelling']=np.log(y.swelling)
    groups=['donor','experiment','forskolin_concentration_µM']
    blocks=[]
    for key,z in y.groupby(groups):
        q=z.groupby(['condition','size']).log_swelling.agg(['mean','size']); cells=[(c,s) for c in ['DMSO','VX445_VX661_VX770'] for s in ['small','large']]
        if any(c not in q.index or q.loc[c,'size']<3 for c in cells):continue
        effect=(q.loc[('VX445_VX661_VX770','large'),'mean']-q.loc[('DMSO','large'),'mean'])-(q.loc[('VX445_VX661_VX770','small'),'mean']-q.loc[('DMSO','small'),'mean'])
        blocks.append(dict(zip(groups,[str(v) for v in key]),effect=float(effect),cell_counts={f'{c}:{s}':int(q.loc[(c,s),'size']) for c,s in cells}))
    by={d:[] for d in sorted(y.donor.unique())}
    for b in blocks:by[b['donor']].append(b['effect'])
    rows=[{'donor':d,'n_blocks':len(v),'median_effect':float(np.median(v)) if v else None} for d,v in by.items()]
    eligible=[r for r in rows if r['n_blocks']];n=len(eligible);pos=sum(r['median_effect']>0 for r in eligible)
    p=binomtest(pos,n,.5,alternative='greater').pvalue if n else None
    out={'source':'OrgaSegment DIS, Zenodo 10610438, same discovery accession',
         'preregistration':'results/preregistration_blocked_size.md (commit 87ec71a)',
         'blocks':blocks,'per_donor':rows,'n_eligible_blocks':len(blocks),
         'n_eligible_donors':n,'n_positive_donors':pos,'n_excluded_donors':len(rows)-n,
         'median_block_effect':float(np.median([b['effect'] for b in blocks])) if blocks else None,
         'sign_p_one_sided':float(p) if p is not None else None}
    out['H1']={'verdict':'UNINTERPRETABLE' if n<8 else 'PASS' if pos>n/2 and p<.05 else 'FAIL'}
    return out

if __name__=='__main__':
    out=assess(pd.read_csv('data/raw/orgasegment/dis_merged_A0.csv'))
    Path('results/blocked_size.json').write_text(json.dumps(out,indent=2)+'\n')
    print({k:v for k,v in out.items() if k not in ('blocks','per_donor')});print(out['per_donor'])
