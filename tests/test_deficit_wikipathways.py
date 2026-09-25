import sys
import pandas as pd
sys.path.insert(0,'src')
from deficit_wikipathways import select,assess

def test_fixed_name_filter_and_low_coverage():
    p=[['Fatty acid oxidation','https://example.test','1'],['Cell cycle','https://example.test','2']]
    assert len(select(p))==1
    d=pd.DataFrame({'organ':['Liver','Brain - Cortex'],'symbol':['A','B']})
    m=pd.DataFrame({'GeneID':['1'],'Symbol':['A']})
    assert assess(d,p,m)['H1']['verdict']=='UNINTERPRETABLE'
