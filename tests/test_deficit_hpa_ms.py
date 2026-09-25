import sys
import pandas as pd
sys.path.insert(0,'src')
from deficit_hpa_ms import assess

def test_missing_reference_tissue_uninterpretable():
    d=pd.DataFrame({'organ':['Liver','Brain - Cortex'],'ens':['E1','E2'],'symbol':['A','B']})
    x=pd.DataFrame({'Gene':['E1','E2'],'Tissue':['liver','liver'],'Intensity':[1,2]})
    assert assess(d,x)['H1']['verdict']=='UNINTERPRETABLE'
