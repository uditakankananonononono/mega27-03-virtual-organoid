"""HPA v25.1 adult tissue mass-spectrometry validation, preregistration_hpa_protein.md."""
import json
from zipfile import ZipFile
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import fisher_exact,mannwhitneyu
PATH='data/hpa_ms/ms_tissue.tsv.zip'
URL='https://www.proteinatlas.org/download/tsv/ms_tissue.tsv.zip'

def assess(deficits, ms):
    p=ms[ms.Tissue.isin(['liver','cerebral cortex'])].pivot(index='Gene',columns='Tissue',values='Intensity')
    out={'source':URL,'preregistration':'results/preregistration_hpa_protein.md (commit 4f05337)',
         'units':'HPA MS protein intensity, undefined absolute scale','per_organ':{}}
    if not {'liver','cerebral cortex'}<=set(p):
        out['H1']={'verdict':'UNINTERPRETABLE','reason':'liver or cerebral cortex absent'}
        out['H2']=out['H1'];return out
    values={}
    for o,g in deficits.groupby('organ'):
        idx=list(g.ens.drop_duplicates()); a=p.reindex(idx)
        liver=a['liver']; brain=a['cerebral cortex']; both=(liver>0)&(brain>0)
        contrast=np.log2(liver[both]/brain[both])
        out['per_organ'][o]={'n_input':len(idx),'n_mapped':int(a.notna().any(axis=1).sum()),
                             'n_liver_detected':int((liver>0).sum()),'n_both_detected':int(both.sum()),
                             'median_log2_liver_over_brain_when_both':float(contrast.median()) if len(contrast) else None,
                             'liver_detected_genes':g.loc[(liver>0).to_numpy(),'symbol'].dropna().tolist()}
        values[o]=contrast
    L=out['per_organ']['Liver'];B=out['per_organ']['Brain - Cortex']
    if min(L['n_mapped'],B['n_mapped'])<50:
        out['H1']=out['H2']={'verdict':'UNINTERPRETABLE','reason':'<50 measured proteins in liver or brain list'}
    else:
        tab=[[L['n_liver_detected'],L['n_input']-L['n_liver_detected']],
             [B['n_liver_detected'],B['n_input']-B['n_liver_detected']]]
        OR,pv=fisher_exact(tab,alternative='greater')
        out['H1']={'table':tab,'odds_ratio':OR,'p_one_sided':pv,'verdict':'PASS' if pv<.05 else 'FAIL'}
        if min(L['n_both_detected'],B['n_both_detected'])<50:
            out['H2']={'verdict':'UNINTERPRETABLE','reason':'<50 pairwise measured proteins in one list'}
        else:
            u,pv=mannwhitneyu(values['Liver'],values['Brain - Cortex'],alternative='greater')
            out['H2']={'U':float(u),'p_one_sided':float(pv),'verdict':'PASS' if pv<.05 else 'FAIL'}
    return out

if __name__=='__main__':
    with ZipFile(PATH) as z:ms=pd.read_csv(z.open('ms_tissue.tsv'),sep='\t')
    d=pd.read_csv('results/strict_organ_deficits.csv');out=assess(d,ms)
    Path('results/deficit_hpa_ms.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['H1'],out['H2']);print({k:{t:v[t] for t in ('n_mapped','n_liver_detected','n_both_detected','median_log2_liver_over_brain_when_both')} for k,v in out['per_organ'].items()})
