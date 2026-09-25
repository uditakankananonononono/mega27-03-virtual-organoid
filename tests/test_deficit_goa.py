import sys
import pandas as pd
sys.path.insert(0,'src')
from deficit_goa import assess

def test_not_qualifier_excluded_and_low_coverage():
    d=pd.DataFrame({'organ':['Liver','Brain - Cortex'],'symbol':['A','B']})
    g=pd.DataFrame({'DB_Object_Symbol':['A','B'],'Qualifier':['NOT|involved_in','involved_in'],
                    'GO_ID':['GO:0006631','GO:0006631'],'Aspect':['P','P'],'Evidence_Code':['IDA','IEA']})
    x=assess(d,g)
    assert x['per_organ']['Liver']['n_direct_metabolic_union']==0
    assert x['H1']['verdict']=='UNINTERPRETABLE'
