"""Direct HPO metabolism/homeostasis term on top organoid deficits, preregistration_hpo.md."""
import json
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact
SOURCE='https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt'
TERM='HP:0001939'


def assess(deficits, hpo):
    annotated=set(hpo.gene_symbol.dropna())
    metabolic=set(hpo.loc[hpo.hpo_id==TERM,'gene_symbol'].dropna())
    rows={}
    for organ,g in deficits.groupby('organ'):
        genes=set(g.symbol.dropna())
        rows[organ]={'n_input':len(genes),'n_any_hpo':len(genes&annotated),
                     'n_direct_metabolism':len(genes&metabolic),
                     'direct_metabolism_genes':sorted(genes&metabolic)}
    L,B=rows['Liver'],rows['Brain - Cortex']
    out={'source':SOURCE,'release':'v2026-09-01 (release resolved by source redirect)',
         'preregistration':'results/preregistration_hpo.md (commit 16771b5)',
         'term':TERM,'per_organ':rows,'limitations':'direct term only; no ancestors; absence from annotation is not evidence of absence'}
    if min(L['n_any_hpo'],B['n_any_hpo'])<50:
        out['H1']={'verdict':'UNINTERPRETABLE','reason':'<50/100 genes covered by HPO in liver or brain'}
    else:
        tab=[[L['n_direct_metabolism'],L['n_input']-L['n_direct_metabolism']],
             [B['n_direct_metabolism'],B['n_input']-B['n_direct_metabolism']]]
        OR,p=fisher_exact(tab,alternative='greater')
        out['H1']={'table':tab,'odds_ratio':OR,'p_one_sided':p,'verdict':'PASS' if p<.05 else 'FAIL'}
    return out

if __name__=='__main__':
    d=pd.read_csv('results/strict_organ_deficits.csv')
    h=pd.read_csv('data/hpo/genes_to_phenotype.txt',sep='\t',usecols=['gene_symbol','hpo_id'])
    out=assess(d,h)
    Path('results/deficit_hpo.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
