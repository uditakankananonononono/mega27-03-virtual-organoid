"""Fixed WikiPathways metabolic union test, preregistration_wikipathways.md."""
import json
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact
ROOT=Path('data/wikipathways')
TERMS=('bile acid','fatty acid','drug metabolism','xenobiotic metabolism')
URL='https://data.wikipathways.org/20260910/gmt/wikipathways-20260910-gmt-Homo_sapiens.gmt'
MAPURL='https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'

def select(lines):
    return [x for x in lines if any(t in x[0].lower() for t in TERMS)]

def assess(deficits,pathways,mapper):
    chosen=select(pathways)
    geneids=set(v for x in chosen for v in x[2:])
    # A symbol must have a unique NCBI GeneID in human gene_info.
    m=mapper[mapper.GeneID.astype(str).isin(geneids)].copy()
    unique=m.loc[~m.Symbol.duplicated(keep=False)]
    names=set(unique.Symbol)
    all_known=set(mapper.Symbol)
    result={'source':URL,'gene_map_source':MAPURL,
      'preregistration':'results/preregistration_wikipathways.md (commit 06f309c)',
      'selected_pathways':[{'name':x[0],'url':x[1],'n_entrez':len(set(x[2:]))} for x in chosen],
      'n_selected_pathways':len(chosen),'n_union_gene_ids':len(geneids),
      'n_union_symbols':len(names),'per_organ':{}}
    for o,g in deficits.groupby('organ'):
        s=set(g.symbol.dropna())
        result['per_organ'][o]={'n_input':len(s),'n_mapped_ncbi':len(s&all_known),
                               'n_union':len(s&names),'matching_symbols':sorted(s&names)}
    l,b=result['per_organ']['Liver'],result['per_organ']['Brain - Cortex']
    if len(chosen)<2 or min(l['n_mapped_ncbi'],b['n_mapped_ncbi'])<80:
        result['H1']={'verdict':'UNINTERPRETABLE','reason':'insufficient pathways or gene mapping'}
    else:
        tab=[[l['n_union'],l['n_input']-l['n_union']],[b['n_union'],b['n_input']-b['n_union']]]
        odds,p=fisher_exact(tab,alternative='greater')
        result['H1']={'table':tab,'odds_ratio':odds,'p_one_sided':p,'verdict':'PASS' if p<.05 else 'FAIL'}
    return result

if __name__=='__main__':
    ps=[l.rstrip('\n').split('\t') for l in open(ROOT/'human_20260910.gmt')]
    mp=pd.read_csv(ROOT/'Homo_sapiens.gene_info.gz',sep='\t',usecols=['GeneID','Symbol'],dtype=str)
    d=pd.read_csv('results/strict_organ_deficits.csv')
    out=assess(d,ps,mp)
    Path('results/deficit_wikipathways.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out['H1'],indent=2));print(out['per_organ'])
