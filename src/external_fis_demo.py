"""External single-plate FIS tracked-object analysis, preregistration_fis_external.md."""
import json,re,requests,numpy as np,pandas as pd
from pathlib import Path
from scipy.stats import binomtest
API='https://api.github.com/repos/hmbotelho/FIS_image_analysis/git/trees/HEAD?recursive=1'
RAW='https://raw.githubusercontent.com/hmbotelho/FIS_image_analysis/master/'
PREFIX='demo_dataset/05-images_analysis/demoplate_01--cellprofiler/'
COLUMNS=['Metadata_compound','Metadata_concentration','Metadata_timeNum','Metadata_wellNum','TrackObjects_Label_4','Math_area_micronsq']

def selected(paths):
    files=[p for p in paths if p.startswith(PREFIX) and p.endswith('/objects.csv') and ('--fsk--' in p or '--fsk_770_809--' in p)]
    return sorted(files)

def tracks(tables):
    chunks=[]; audit=[]
    for path,t in tables:
        d=t[COLUMNS].copy();d['path']=path
        d['track']=pd.to_numeric(d['TrackObjects_Label_4'],errors='coerce')
        d['time']=pd.to_numeric(d['Metadata_timeNum'],errors='coerce')
        d['area']=pd.to_numeric(d['Math_area_micronsq'],errors='coerce')
        d['dose']=pd.to_numeric(d['Metadata_concentration'],errors='coerce')
        d['well']=d['Metadata_wellNum'].astype(str).str.zfill(4)
        d['condition']=d['Metadata_compound'].astype(str)
        valid=d.track.notna()&d.time.isin([0,6])&np.isfinite(d.area)&(d.area>0)&d.dose.notna()
        d=d[valid].copy()
        dup=d.duplicated(['well','track','time'],keep=False)
        audit.append({'path':path,'rows':len(t),'candidate_rows':len(d),'duplicate_track_time_rows':int(dup.sum())})
        chunks.append(d[~dup])
    a=pd.concat(chunks,ignore_index=True)
    keys=['well','track','condition','dose']
    z=a[a.time.isin([0,6])].pivot_table(index=keys,columns='time',values='area',aggfunc='first')
    z=z.dropna(subset=[0,6]).reset_index();z=z.rename(columns={0:'A0',6:'A60'})
    z['log_swelling']=np.log(z.A60/z.A0)
    return z,audit

def assess(z,audit,paths):
    assert set(z.condition)<=set(['fsk','fsk_770_809'])
    low,high=z.A0.quantile([.25,.75]);z=z.copy()
    z['size']=np.select([z.A0<=low,z.A0>=high],['small','large'],default='middle')
    doses=[]
    for dose,g in z.groupby('dose'):
        counts={}
        for condition in ['fsk','fsk_770_809']:
            for size in ['small','large']:
                t=g[(g.condition==condition)&(g['size']==size)]
                counts[f'{condition}:{size}']={'tracks':int(len(t)),'wells':int(t.well.nunique()),'mean_log_swelling':float(t.log_swelling.mean()) if len(t) else None}
        ok=all(v['tracks']>=5 for v in counts.values()) and all(g[g.condition==c].well.nunique()>=2 for c in ['fsk','fsk_770_809'])
        effect=None
        if ok:
            m=lambda c,s:counts[f'{c}:{s}']['mean_log_swelling']
            effect=(m('fsk_770_809','large')-m('fsk','large'))-(m('fsk_770_809','small')-m('fsk','small'))
        doses.append({'dose_uM':float(dose),'eligible':bool(ok),'effect_log':float(effect) if effect is not None else None,'cells':counts,
                      'n_wells':{c:int(g[g.condition==c].well.nunique()) for c in ['fsk','fsk_770_809']}})
    effects=[d['effect_log'] for d in doses if d['eligible']];n=len(effects);pos=sum(v>0 for v in effects)
    p=binomtest(pos,n,.5,alternative='greater').pvalue if n else None
    return {'source':'https://github.com/hmbotelho/FIS_image_analysis','preregistration':'results/preregistration_fis_external.md (commit bdc3fcb)',
            'paths':paths,'source_audit':audit,'n_tracked_pairs':int(len(z)),'quartile_cutoffs_area_micronsq':{'small_lte':float(low),'large_gte':float(high)},
            'doses':doses,'n_eligible_doses':n,'n_positive_doses':pos,'median_dose_effect':float(np.median(effects)) if n else None,
            'sign_p_one_sided':float(p) if p is not None else None,'H1':{'verdict':'UNINTERPRETABLE' if n<4 else 'PASS' if np.median(effects)>0 and p<.05 else 'FAIL'}}

if __name__=='__main__':
    d=Path('data/external_fis_demo');d.mkdir(parents=True,exist_ok=True)
    resp=requests.get(API,timeout=20);resp.raise_for_status();tree=resp.json();assert not tree.get('truncated')
    paths=selected(x['path'] for x in tree['tree']);print('source csvs',len(paths),flush=True)
    tables=[]
    for p in paths:
        dest=d/(re.search(r'W\d{4}',p).group(0)+'.csv')
        if not dest.exists():
            r=requests.get(RAW+p,timeout=12);r.raise_for_status();dest.write_bytes(r.content)
        tables.append((p,pd.read_csv(dest,usecols=COLUMNS)))
    z,audit=tracks(tables);out=assess(z,audit,paths)
    Path('results/external_fis_demo.json').write_text(json.dumps(out,indent=2)+'\n')
    print('pairs',len(z),'H1',out['H1'],'positive',out['n_positive_doses'],'of',out['n_eligible_doses'],'p',out['sign_p_one_sided'])
    print([(a['dose_uM'],a['eligible'],a['effect_log']) for a in out['doses']])
