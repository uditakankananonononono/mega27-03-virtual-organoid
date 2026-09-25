"""Fixed 12-target InterPro domain audit, preregistration_interpro.md."""
import io,json,time,requests,pandas as pd
from pathlib import Path
TARGETS={'liver':['CYP2D6','CYP2C9','CYP2C19','UGT1A1','SLCO1B1','VKORC1'],
         'brain':['OPALIN','MAG','C1QB','HBB','HBA2','S100A9']}
URL='https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/'
CATEGORIES=('cytochrome p450','udp-glucuronosyltransferase','organic anion transport','solute carrier','vitamin k epoxide reductase')

def classify(names):
    """Conservative keyword match from actual InterPro names, not a biological truth test."""
    return sorted({c for c in CATEGORIES if any(c in n.lower() for n in names)})

def run(session=requests):
    results={}
    for group, genes in TARGETS.items():
        for g in genes:
            q={'query':f'gene_exact:{g} AND organism_id:9606 AND reviewed:true',
               'fields':'accession,gene_primary','format':'tsv','size':10}
            u=session.get('https://rest.uniprot.org/uniprotkb/search',params=q,timeout=20);u.raise_for_status()
            tab=pd.read_csv(io.StringIO(u.text),sep='\t')
            acc=tab['Entry'].tolist() if 'Entry' in tab else []
            d={'group':group,'gene':g,'uniprot_query_url':u.url,'reviewed_accessions':acc}
            if len(acc)==1:
                endpoint=URL+acc[0]+'/?page_size=100';r=session.get(endpoint,timeout=25);r.raise_for_status();j=r.json()
                entries=[{'accession':x['metadata']['accession'],'name':x['metadata']['name'],'type':x['metadata'].get('type')} for x in j['results']]
                d.update({'interpro_url':r.url,'entry_count':j['count'],'entries':entries,
                          'classifications':classify([x['name'] for x in entries]),
                          'truncated':bool(j.get('next'))})
            else:d['unresolved']='reviewed mapping not one-to-one'
            results[g]=d
            print(g,acc,[(x['accession'],x['name']) for x in d.get('entries',[])][:5],flush=True)
            time.sleep(.15)
    deficit=pd.read_csv('results/strict_organ_deficits.csv')
    membership={group:[g for g in genes if g in set(deficit[deficit.organ==('Liver' if group=='liver' else 'Brain - Cortex')].symbol)] for group,genes in TARGETS.items()}
    hit=[g for g in membership['liver'] if results[g].get('classifications')]
    brain=[g for g in membership['brain'] if results[g].get('classifications')]
    out={'source':URL,'preregistration':'results/preregistration_interpro.md (commit cdb1051)',
         'targets':results,'deficit_top100_membership':membership,
         'H1':{'liver_matching_targets':hit,'brain_matching_targets':brain,
               'verdict':'INVALID - pre-registration incorrectly called 5 of 6 liver targets top-100 deficits'},
         'caveat':'The preregistration mistakenly designated six drug-metabolism genes as top-100 liver deficits; only CYP2D6 is in that top 100. The six gene domain annotations are correct but do not test H1 on top deficits. No success claim.'}
    Path('results/deficit_interpro.json').write_text(json.dumps(out,indent=2)+'\n')
    return out

if __name__=='__main__':run()
