"""One-plate FIS demo well-level and attrition sensitivity, prereg bf80ddc."""
import json,glob,numpy as np,pandas as pd
from pathlib import Path
from external_fis_demo import tracks,COLUMNS
LOW=2042.626642;HIGH=6700.811789

def assess(tables,boots=2000):
    z,audit=tracks(tables)
    z['size']=np.select([z.A0<=LOW,z.A0>=HIGH],['small','large'],default='middle')
    wells=[]
    for (dose,condition,well),v in z.groupby(['dose','condition','well']):
        counts={size:int(sum(v['size']==size)) for size in ['small','large']}
        if min(counts.values())<5:effect=None
        else:effect=float(v.loc[v['size']=='large','log_swelling'].mean()-v.loc[v['size']=='small','log_swelling'].mean())
        wells.append({'dose_uM':float(dose),'condition':condition,'well':well,'size_counts':counts,'large_minus_small':effect})
    eligible=[]
    for dose in sorted(set(r['dose_uM'] for r in wells)):
        w=[r for r in wells if r['dose_uM']==dose and r['large_minus_small'] is not None]
        by={c:[r['large_minus_small'] for r in w if r['condition']==c] for c in ['fsk','fsk_770_809']}
        if min(map(len,by.values()))<2:continue
        eligible.append({'dose_uM':dose,'effect':float(np.mean(by['fsk_770_809'])-np.mean(by['fsk'])),'well_effects':by})
    rng=np.random.default_rng(27)
    draws=[]
    if eligible:
        for _ in range(boots):
            e=[np.mean(rng.choice(row['well_effects']['fsk_770_809'],size=2,replace=True))-np.mean(rng.choice(row['well_effects']['fsk'],size=2,replace=True)) for row in eligible]
            draws.append(float(np.median(e)))
    ci=np.quantile(draws,[.025,.975]).tolist() if draws else None
    R1={'n_eligible_doses':len(eligible),'n_positive':sum(x['effect']>0 for x in eligible),
        'median_effect':float(np.median([x['effect'] for x in eligible])) if eligible else None,
        'well_bootstrap_ci95':ci,'dose_effects':eligible,'well_rows':wells,
        'verdict':'UNINTERPRETABLE' if len(eligible)<4 else 'SUPPORTED' if ci[0]>0 and np.median([x['effect'] for x in eligible])>0 else 'FAIL'}
    attr={c:{s:{'eligible_baseline':0,'retained':0} for s in ['small','large']} for c in ['fsk','fsk_770_809']}
    missing={'baseline_rows_missing_track_id':0,'duplicate_baseline_track_rows':0,'duplicate_endpoint_track_rows':0}
    for path,t in tables:
        d=t[COLUMNS].copy();d['track']=pd.to_numeric(d.TrackObjects_Label_4,errors='coerce')
        d['area']=pd.to_numeric(d.Math_area_micronsq,errors='coerce')
        d['time']=pd.to_numeric(d.Metadata_timeNum,errors='coerce')
        d['condition']=d.Metadata_compound.astype(str)
        b=d[(d.time==0)&np.isfinite(d.area)&(d.area>0)];missing['baseline_rows_missing_track_id']+=int(b.track.isna().sum())
        b=b[b.track.notna()].copy();duplicates=b.duplicated('track',keep=False);missing['duplicate_baseline_track_rows']+=int(duplicates.sum());b=b[~duplicates]
        e=d[(d.time==6)&d.track.notna()&np.isfinite(d.area)&(d.area>0)].copy()
        duplicates=e.duplicated('track',keep=False);missing['duplicate_endpoint_track_rows']+=int(duplicates.sum());e=e[~duplicates]
        ids=set(e.track);c=str(d.condition.iloc[0]);assert c in attr
        for size,sel in [('small',b.area<=LOW),('large',b.area>=HIGH)]:
            v=b[sel];attr[c][size]['eligible_baseline']+=len(v);attr[c][size]['retained']+=int(v.track.isin(ids).sum())
    for c in attr:
        for s in attr[c]:
            a=attr[c][s];a['rate']=a['retained']/a['eligible_baseline'] if a['eligible_baseline'] else None
    within={c:(attr[c]['small']['rate']-attr[c]['large']['rate']) if all(attr[c][size]['rate'] is not None for size in ['small','large']) else None for c in attr}
    cross=(within['fsk_770_809']-within['fsk']) if all(v is not None for v in within.values()) else None
    flag=None if cross is None else any(abs(v)>.10 for v in within.values()) or abs(cross)>.10
    R2={'by_condition_size':attr,'small_minus_large_retention':within,'cross_condition_difference':cross,
        'material_attrition_flag':flag,'missing_and_duplicates':missing}
    return {'source':'https://github.com/hmbotelho/FIS_image_analysis','preregistration':'results/preregistration_external_fis_robustness.md (commit bf80ddc)',
            'size_area_micronsq':{'small_lte':LOW,'large_gte':HIGH},'R1':R1,'R2':R2,'source_audit':audit}

if __name__=='__main__':
    prior=json.load(open('results/external_fis_demo.json'))
    tables=[(p,pd.read_csv('data/external_fis_demo/'+__import__('re').search(r'W\d{4}',p).group(0)+'.csv',usecols=COLUMNS)) for p in prior['paths']]
    out=assess(tables)
    Path('results/external_fis_robustness.json').write_text(json.dumps(out,indent=2)+'\n')
    print({k:v for k,v in out['R1'].items() if k not in ['dose_effects','well_rows']})
    print(out['R2'])
