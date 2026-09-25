"""JASPAR exact-name motif coverage for pre-fixed CollecTRI gaps."""
import json,requests
from pathlib import Path
URL='https://jaspar.elixir.no/api/v1/matrix/'
NAMES=['NEUROD2','NEUROD6','TBR1','PAX2']

def pick(name,records):
    exact=[r for r in records if r['name'].upper()==name.upper()]
    if not exact:return None
    return sorted(exact,key=lambda r:(int(r.get('version') or 0),r['matrix_id']))[-1]

def assess(found):
    hits={x:found[x]['matrix_id'] if found[x] else None for x in NAMES}
    yes=sum(hits[x] is not None for x in NAMES[:3])>=2 and hits['PAX2'] is not None
    return {'source':URL,'preregistration':'results/preregistration_jaspar.md (commit 0ebc540)',
            'selected':found,'exact_motifs':hits,'H1':{'verdict':'PASS' if yes else 'FAIL','n_brain':sum(hits[x] is not None for x in NAMES[:3]),'PAX2':hits['PAX2'] is not None}}

if __name__=='__main__':
    found={}
    for name in NAMES:
        url=URL+'?collection=CORE&tax_group=vertebrates&search='+name+'&page_size=100'
        try:
            r=requests.get(url,timeout=12);r.raise_for_status();response=r.json()
            if response.get('next'):
                raise ValueError('query has more than 100 matches; do not classify as absence')
            hit=pick(name,response['results']);found[name]={'query_url':url,'n_search_hits':response['count'],
                'matrix_id':hit['matrix_id'],'matrix_url':hit['url'],'name':hit['name']} if hit else None
        except Exception as e:
            found[name]={'query_url':url,'error':type(e).__name__+': '+str(e)[:150]}
    # Errors are not absence; avoid accidental PASS when an exact record is unavailable.
    if any(x and 'error' in x for x in found.values()):
        out={'source':URL,'preregistration':'results/preregistration_jaspar.md (commit 0ebc540)','selected':found,'H1':{'verdict':'UNINTERPRETABLE','reason':'source query error'}}
    else:out=assess(found)
    Path('results/deficit_jaspar.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out)
