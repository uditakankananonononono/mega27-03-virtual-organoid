import sys
sys.path.insert(0,'src')
import pandas as pd
from deficit_rhea import assess

def test_symbol_mapping_and_coverage():
    d=pd.DataFrame({'organ':['Liver','Brain - Cortex'],'symbol':['A','B']})
    g=pd.DataFrame({'DB':['UniProtKB']*2,'Taxon':['taxon:9606']*2,'DB_Object_Symbol':['A','B'],'DB_Object_ID':['P1','P2']})
    r=pd.DataFrame({'ID':['P1'],'RHEA_ID':['123']})
    a=assess(d,g,r)
    assert a['per_organ']['Liver']['n_rhea_positive']==1
    assert a['per_organ']['Brain - Cortex']['n_rhea_positive']==0
    assert a['H1']['verdict']=='UNINTERPRETABLE'
