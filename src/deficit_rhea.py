"""Rhea curated reactions linked to Swiss-Prot, preregistration_rhea.md."""
import json,requests,pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact
from deficit_goa import COLS
URL='https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_sprot.tsv'

def assess(deficits,goa,rhea):
    human=goa[(goa.DB=='UniProtKB') & (goa.Taxon.fillna('').str.split('|').str[0]=='taxon:9606')]
    mapping=human.groupby('DB_Object_Symbol').DB_Object_ID.agg(lambda a:set(a)).to_dict()
    rmap=rhea.groupby('ID').RHEA_ID.agg(lambda a:set(a.astype(str))).to_dict()
    out={'source':URL,'preregistration':'results/preregistration_rhea.md (commit 9c649a5)',
         'n_rhea_rows':len(rhea),'n_rhea_accessions':len(rmap),'per_organ':{}}
    for o,g in deficits.groupby('organ'):
        genes=set(g.symbol.dropna());hit={x for x in genes if mapping.get(x,set())&rmap.keys()}
        ex={x:sorted(set().union(*(rmap.get(a,set()) for a in mapping.get(x,set()))))[:3] for x in sorted(hit)}
        out['per_organ'][o]={'n_input':len(genes),'n_goa_uniprot_mapped':sum(bool(mapping.get(x)) for x in genes),
                            'unmapped_symbols':sorted(x for x in genes if not mapping.get(x)),
                            'n_rhea_positive':len(hit),'positive_symbols_rhea_examples':ex}
    L=out['per_organ']['Liver'];B=out['per_organ']['Brain - Cortex']
    if min(L['n_goa_uniprot_mapped'],B['n_goa_uniprot_mapped'])<80 or sum(v['n_rhea_positive'] for v in out['per_organ'].values())<20:
        out['H1']={'verdict':'UNINTERPRETABLE','reason':'registered mapping/positive coverage floor'}
    else:
        table=[[L['n_rhea_positive'],L['n_input']-L['n_rhea_positive']],[B['n_rhea_positive'],B['n_input']-B['n_rhea_positive']]]
        OR,p=fisher_exact(table,alternative='greater')
        out['H1']={'table':table,'odds_ratio':OR,'p_one_sided':p,'verdict':'PASS' if p<0.05 else 'FAIL'}
    return out

if __name__=='__main__':
    path=Path('data/rhea/rhea2uniprot_sprot.tsv');path.parent.mkdir(exist_ok=True)
    if not path.exists():
        resp=requests.get(URL,timeout=35);resp.raise_for_status();path.write_bytes(resp.content)
    d=pd.read_csv('results/strict_organ_deficits.csv')
    g=pd.read_csv('data/goa/goa_human.gaf.gz',sep='\t',comment='!',header=None,names=COLS,dtype=str)
    r=pd.read_csv(path,sep='\t',dtype=str)
    out=assess(d,g,r);Path('results/deficit_rhea.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['H1']);print({o:(v['n_goa_uniprot_mapped'],v['n_rhea_positive']) for o,v in out['per_organ'].items()})
