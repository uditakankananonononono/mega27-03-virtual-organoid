"""Fixed Europe PMC searches preregistered in results/preregistration_epmc_priorart.md."""
import json, requests
from pathlib import Path
URL='https://www.ebi.ac.uk/europepmc/webservices/rest/search'
Q=[
'TITLE_ABS:"organoid" AND TITLE_ABS:"forskolin" AND TITLE_ABS:"size" AND TITLE_ABS:"swelling"',
'TITLE_ABS:"organoid" AND TITLE_ABS:"forskolin" AND TITLE_ABS:"lumen" AND TITLE_ABS:"swelling"',
'TITLE_ABS:"organoid" AND TITLE_ABS:"CFTR" AND TITLE_ABS:"size" AND TITLE_ABS:"response"']

def parse(res,query):
    out=[]
    for r in res['resultList']['result']:
        out.append({'id':r.get('id'),'source':r.get('source'), 'pmid':r.get('pmid'), 'doi':r.get('doi'),
                    'title':r.get('title'),'abstract':r.get('abstractText'),
                    'journal':r.get('journalTitle'),'year':r.get('pubYear'),
                    'url':f"https://europepmc.org/article/{r.get('source')}/{r.get('id')}",
                    'retrieved_by':query})
    return out

if __name__=='__main__':
    searched=[]; dedup={}
    for q in Q:
        r=requests.get(URL,params={'query':q,'format':'json','pageSize':50,'resultType':'core'},timeout=30);r.raise_for_status();j=r.json()
        docs=[]
        for x in j['resultList']['result']:
            d={'id':x.get('id'),'source':x.get('source'),'pmid':x.get('pmid'),'doi':x.get('doi'),
               'title':x.get('title'),'abstract':x.get('abstractText'),'year':x.get('pubYear'),
               'url':f"https://europepmc.org/article/{x.get('source')}/{x.get('id')}", 'query':q}
            docs.append(d);dedup.setdefault(x.get('pmid') or x.get('doi') or x.get('id'),d)
        searched.append({'query':q,'hitCount':j['hitCount'],'returned':len(docs),'results':docs})
        print('query hits',j['hitCount'],'returned',len(docs),flush=True)
    out={'source':URL,'preregistration':'results/preregistration_epmc_priorart.md (commit cdb1051)',
         'searches':searched,'n_unique_returned':len(dedup),
         'interpretation_status':'manual title/abstract review needed; search absence cannot prove novelty'}
    Path('results/epmc_priorart.json').write_text(json.dumps(out,indent=2)+'\n')
