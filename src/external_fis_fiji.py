"""Same-plate Fiji tracking pipeline sensitivity, preregistration_fiji_crosspipeline.md."""
import json,re,requests,pandas as pd,numpy as np
from pathlib import Path
from external_fis_demo import tracks,assess
from external_fis_robustness import assess as well_assess
API='https://api.github.com/repos/hmbotelho/FIS_image_analysis/git/trees/HEAD?recursive=1'
RAW='https://raw.githubusercontent.com/hmbotelho/FIS_image_analysis/master/'
PREFIX='demo_dataset/05-images_analysis/demoplate_01--ij/'
ORIG=['Metadata_compound','Metadata_concentration','Metadata_timeNum','Metadata_wellNum','TrackObjects_Label','Math_area_micronsq']
COLUMNS=['Metadata_compound','Metadata_concentration','Metadata_timeNum','Metadata_wellNum','TrackObjects_Label_4','Math_area_micronsq']

def selected(paths):return sorted(p for p in paths if p.startswith(PREFIX) and p.endswith('/objects.csv') and ('--fsk--' in p or '--fsk_770_809--' in p))

def analyze(tables):
    normalized=[]
    for p,t in tables:normalized.append((p,t.rename(columns={'TrackObjects_Label':'TrackObjects_Label_4'})))
    z,audit=tracks(normalized);main=assess(z,audit,[p for p,_ in tables])
    well=well_assess(normalized,boots=2000)
    elig=[x for x in main['doses'] if x['eligible']]
    main['H1_fiji']={'verdict':'UNINTERPRETABLE' if len(elig)<4 else 'PASS' if np.median([x['effect_log'] for x in elig])>0 and sum(x['effect_log']>0 for x in elig)>=6 else 'FAIL',
                     'n_eligible':len(elig),'n_positive':sum(x['effect_log']>0 for x in elig),
                     'median_effect':float(np.median([x['effect_log'] for x in elig])) if elig else None}
    main['well_and_retention_descriptive']={'n_eligible_doses':well['R1']['n_eligible_doses'],
        'n_positive_doses':well['R1']['n_positive'],'median_equal_well_effect':well['R1']['median_effect'],
        'retention_by_condition_size':well['R2']['by_condition_size'],'missing_and_duplicates':well['R2']['missing_and_duplicates']}
    main['preregistration']='results/preregistration_fiji_crosspipeline.md (commit 4aaf86f)'
    main['interpretation']='same donor/plate/images; alternative analysis pipeline, not biological replication'
    return main

if __name__=='__main__':
    d=Path('data/external_fis_fiji');d.mkdir(parents=True,exist_ok=True)
    r=requests.get(API,timeout=15);r.raise_for_status();tree=r.json();assert not tree.get('truncated')
    paths=selected(x['path'] for x in tree['tree']);print('selected',len(paths),flush=True)
    tables=[]
    for p in paths:
        dest=d/(re.search(r'W\d{4}',p).group(0)+'.csv')
        if not dest.exists():
            r=requests.get(RAW+p,timeout=12);r.raise_for_status();dest.write_bytes(r.content)
        tables.append((p,pd.read_csv(dest,usecols=ORIG)))
    out=analyze(tables);Path('results/external_fis_fiji.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['H1_fiji'],out['n_tracked_pairs'],out['quartile_cutoffs_area_micronsq'])
    print([(x['dose_uM'],x['effect_log']) for x in out['doses']])
    print(out['well_and_retention_descriptive'])
