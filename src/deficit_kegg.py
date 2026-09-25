"""KEGG REST direct pathway-gene annotation, preregistration_kegg.md."""
import json,requests,pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact
BASE='https://rest.kegg.jp'
TERMS=('fatty acid','bile acid','xenobiotic','drug metabolism')

def paths(text):
    return [{'id':l.split('\t')[0], 'name':l.split('\t')[1],
             'url':BASE+'/link/hsa/'+l.split('\t')[0]}
            for l in text.splitlines() if '\t' in l and any(t in l.split('\t')[1].lower() for t in TERMS)]

def assess(deficits,pathways,ncbi):
    E=set(g for p in pathways for g in p['gene_ids'])
    subset=ncbi[ncbi.GeneID.astype(str).isin(E)]
    symbols=set(subset.loc[~subset.Symbol.duplicated(keep=False),'Symbol'])
    known=set(ncbi.Symbol)
    out={'source':BASE,'preregistration':'results/preregistration_kegg.md (commit 32f57b9)',
         'n_selected_pathways':len(pathways),'selected_pathways':pathways,
         'n_union_entrez':len(E),'n_union_symbols':len(symbols),'per_organ':{}}
    for o,g in deficits.groupby('organ'):
        s=set(g.symbol.dropna());out['per_organ'][o]={'n_input':len(s),'n_mapped_ncbi':len(s&known),
             'n_union':len(s&symbols),'symbols':sorted(s&symbols)}
    L,B=out['per_organ']['Liver'],out['per_organ']['Brain - Cortex']
    if len(pathways)<2 or min(L['n_mapped_ncbi'],B['n_mapped_ncbi'])<80:
        out['H1']={'verdict':'UNINTERPRETABLE','reason':'pathways/mapping below preregistered floor'}
    else:
        tab=[[L['n_union'],L['n_input']-L['n_union']],[B['n_union'],B['n_input']-B['n_union']]]
        odds,p=fisher_exact(tab,alternative='greater');out['H1']={'table':tab,'odds_ratio':odds,'p_one_sided':p,'verdict':'PASS' if p<.05 else 'FAIL'}
    return out

if __name__=='__main__':
    list_url=BASE+'/list/pathway/hsa';r=requests.get(list_url,timeout=20);r.raise_for_status()
    P=paths(r.text)
    for p in P:
        x=requests.get(p['url'],timeout=20);x.raise_for_status()
        p['gene_ids']=sorted({z.split('\t')[1].replace('hsa:','') for z in x.text.splitlines() if '\t' in z})
        print(p['id'],p['name'],len(p['gene_ids']),flush=True)
    d=pd.read_csv('results/strict_organ_deficits.csv')
    mp=pd.read_csv('data/wikipathways/Homo_sapiens.gene_info.gz',sep='\t',usecols=['GeneID','Symbol'],dtype=str)
    out=assess(d,P,mp);out['list_url']=list_url
    Path('results/deficit_kegg.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['H1']);print({k:v['n_union'] for k,v in out['per_organ'].items()})
