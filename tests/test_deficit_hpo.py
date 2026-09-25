import sys
import pandas as pd
sys.path.insert(0,'src')
from deficit_hpo import assess

def test_low_annotation_coverage_uninterpretable():
    d=pd.DataFrame({'organ':['Liver','Brain - Cortex'],'symbol':['A','B']})
    h=pd.DataFrame({'gene_symbol':['A'],'hpo_id':['HP:0001939']})
    assert assess(d,h)['H1']['verdict']=='UNINTERPRETABLE'
