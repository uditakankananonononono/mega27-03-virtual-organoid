"""Top-20 IntAct interaction count audit, preregistration_intact.md."""
import json, requests, pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from scipy.stats import mannwhitneyu
from numpy import log1p
from deficit_goa import COLS
URL='https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/'
ORGANS=('Liver','Brain - Cortex')

def fixed_rows(deficits,goa):
    human=goa[(goa.DB=='UniProtKB') & (goa.Taxon.fillna('').str.split('|').str[0]=='taxon:9606')]
    acc=human.groupby('DB_Object_Symbol').DB_Object_ID.agg(lambda a:sorted(set(a))).to_dict()
    rows=[]
    for organ in ORGANS:
        d=deficits[deficits.organ==organ].head(20)
        for rank,(_,x) in enumerate(d.iterrows(),1):
            a=acc.get(x.symbol,[])
            rows.append({'organ':organ,'rank':rank,'symbol':x.symbol,'accessions':a,
                         'accession':a[0] if len(a)==1 else None,'count':None,'status':'pending' if len(a)==1 else 'ambiguous_or_unmapped'})
    return rows

def query(row):
    if row['accession'] is None:return row
    url=URL+'identifier:'+row['accession']+'?format=count'
    row['url']=url
    try:
        r=requests.get(url,timeout=10);r.raise_for_status();row['count']=int(r.text.strip());row['status']='ok'
    except Exception as e:row['status']='error: '+type(e).__name__+' '+str(e)[:100]
    return row

def assess(rows):
    by={o:[x['count'] for x in rows if x['organ']==o and x['status']=='ok'] for o in ORGANS}
    out={'source':URL,'preregistration':'results/preregistration_intact.md (commit 1c7c0c8)', 'rows':rows,
         'per_organ':{o:{'n_measured':len(v),'median_count':float(pd.Series(v).median()) if v else None,
                         'median_log1p_count':float(pd.Series(log1p(v)).median()) if v else None} for o,v in by.items()}}
    if min(map(len,by.values()))<15:out['H1']={'verdict':'UNINTERPRETABLE','reason':'<15 measured genes in one group'}
    else:
        u,p=mannwhitneyu(log1p(by['Liver']),log1p(by['Brain - Cortex']),alternative='less')
        out['H1']={'u':float(u),'p_one_sided':float(p),'verdict':'PASS' if p<.05 else 'FAIL'}
    return out

if __name__=='__main__':
    d=pd.read_csv('results/strict_organ_deficits.csv')
    g=pd.read_csv('data/goa/goa_human.gaf.gz',sep='\t',comment='!',header=None,names=COLS,dtype=str)
    rows=fixed_rows(d,g)
    with ThreadPoolExecutor(max_workers=4) as pool: rows=list(pool.map(query,rows))
    out=assess(rows);Path('results/deficit_intact.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['per_organ'],out['H1'],flush=True)
    print([(r['organ'],r['symbol'],r['status']) for r in rows if r['status']!='ok'])
