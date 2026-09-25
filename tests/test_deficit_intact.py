import sys
sys.path.insert(0,'src')
import pandas as pd
from deficit_intact import fixed_rows,assess

def test_fixed_selection_ambiguous_and_floor():
    d=pd.DataFrame({'organ':['Liver','Brain - Cortex'],'symbol':['A','B']})
    g=pd.DataFrame({'DB':['UniProtKB']*3,'Taxon':['taxon:9606']*3,'DB_Object_Symbol':['A','A','B'],'DB_Object_ID':['P1','P2','P3']})
    r=fixed_rows(d,g)
    assert r[0]['status']=='ambiguous_or_unmapped' and r[1]['accession']=='P3'
    assert assess(r)['H1']['verdict']=='UNINTERPRETABLE'
