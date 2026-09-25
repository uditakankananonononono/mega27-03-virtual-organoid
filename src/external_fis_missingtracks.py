"""Preregistered missing-track audit and 5th/95th extreme envelope."""
import json,re,numpy as np,pandas as pd
from pathlib import Path
from external_fis_demo import COLUMNS,tracks
LOW=2042.626642;HIGH=6700.811789

def assess(tables):
    z,audit=tracks(tables)
    z['size']=np.select([z.A0<=LOW,z.A0>=HIGH],['small','large'],default='middle')
    z['key']=list(zip(z.well,z.track))
    base=[]
    for path,t in tables:
        d=t[COLUMNS].copy();d['area']=pd.to_numeric(d.Math_area_micronsq,errors='coerce');d['track']=pd.to_numeric(d.TrackObjects_Label_4,errors='coerce')
        d['time']=pd.to_numeric(d.Metadata_timeNum,errors='coerce');d['dose']=pd.to_numeric(d.Metadata_concentration,errors='coerce')
        d['well']=d.Metadata_wellNum.astype(str).str.zfill(4);d['condition']=d.Metadata_compound.astype(str)
        b=d[(d.time==0)&np.isfinite(d.area)&(d.area>0)].copy();b['size']=np.select([b.area<=LOW,b.area>=HIGH],['small','large'],default='middle')
        b['duplicate']=b.track.notna()&b.duplicated('track',keep=False)
        b['usable_label']=b.track.notna()&~b['duplicate']
        base.append(b[['condition','dose','well','track','area','size','usable_label']])
    base=pd.concat(base,ignore_index=True);base['key']=list(zip(base.well,base.track))
    observed=set(z.key);base['paired']=base.usable_label&base.key.isin(observed)
    m1={}
    for (c,size),g in base.groupby(['condition','size']):
        lost=g[g.usable_label&~g.paired];ret=g[g.paired]
        m1.setdefault(c,{})[size]={'all_baseline_objects':len(g),'usable_label':int(g.usable_label.sum()),
            'missing_or_duplicate_label':int((~g.usable_label).sum()),'paired':int(g.paired.sum()),
            'paired_fraction_all_baseline':float(g.paired.mean()),
            'median_area_lost_labeled':float(lost.area.median()) if len(lost) else None,
            'median_area_retained':float(ret.area.median()) if len(ret) else None}
    doses=[]
    for dose,g in base.groupby('dose'):
        cells={}
        for c in ['fsk','fsk_770_809']:
            for size in ['small','large']:
                b=g[(g.condition==c)&(g['size']==size)];matched=z[(z.dose==dose)&(z.condition==c)&(z['size']==size)]
                key=f'{c}:{size}';cells[key]={'n_baseline':len(b),'n_observed':len(matched),'n_missing':len(b)-len(matched)}
                if len(matched)<5:continue
                obs=matched.log_swelling.to_numpy();lo,hi=np.quantile(obs,[.05,.95]);nmiss=len(b)-len(matched)
                cells[key].update({'mean_observed':float(np.mean(obs)),'p05':float(lo),'p95':float(hi),
                    'mean_adverse_low':float((obs.sum()+nmiss*lo)/len(b)),
                    'mean_adverse_high':float((obs.sum()+nmiss*hi)/len(b))})
        eligible=all('mean_observed' in q for q in cells.values())
        adverse=best=None
        if eligible:
            def v(c,size,bound):return cells[f'{c}:{size}']['mean_adverse_'+bound]
            adverse=(v('fsk_770_809','large','low')-v('fsk','large','high'))-(v('fsk_770_809','small','high')-v('fsk','small','low'))
            best=(v('fsk_770_809','large','high')-v('fsk','large','low'))-(v('fsk_770_809','small','low')-v('fsk','small','high'))
        doses.append({'dose_uM':float(dose),'eligible':eligible,'adverse_effect':float(adverse) if adverse is not None else None,
                      'best_effect':float(best) if best is not None else None,'cells':cells})
    v=[x['adverse_effect'] for x in doses if x['eligible']]
    return {'source':'https://github.com/hmbotelho/FIS_image_analysis','preregistration':'results/preregistration_external_fis_missingtracks.md (commit 6a0aed7)',
            'M1':m1,'M2':{'doses':doses,'n_eligible':len(v),'n_adverse_positive':int(sum(x>0 for x in v)),
                        'adverse_median':float(np.median(v)) if v else None,
                        'best_median':float(np.median([x['best_effect'] for x in doses if x['eligible']])) if v else None,
                        'verdict':'UNINTERPRETABLE' if not v else 'ROBUST_WITHIN_ENVELOPE' if np.median(v)>0 else 'VULNERABLE_TO_MISSINGNESS'}}

if __name__=='__main__':
    prior=json.load(open('results/external_fis_demo.json'))
    tables=[(p,pd.read_csv('data/external_fis_demo/'+re.search(r'W\d{4}',p).group(0)+'.csv',usecols=COLUMNS)) for p in prior['paths']]
    out=assess(tables);Path('results/external_fis_missingtracks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['M1']);print({k:v for k,v in out['M2'].items() if k!='doses'});print([(d['dose_uM'],d['adverse_effect'],d['best_effect']) for d in out['M2']['doses']])
