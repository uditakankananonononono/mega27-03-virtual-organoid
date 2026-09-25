import sys
import pandas as pd
sys.path.insert(0,'src')
from deficit_kegg import paths,assess

def test_fixed_pathway_name_filter():
    x=paths('hsa0001\tFatty acid oxidation - Homo sapiens\nhsa0002\tCell cycle - Homo sapiens')
    assert len(x)==1 and x[0]['id']=='hsa0001'
    d=pd.DataFrame({'organ':['Liver','Brain - Cortex'],'symbol':['A','B']})
    m=pd.DataFrame({'GeneID':['1'],'Symbol':['A']})
    x[0]['gene_ids']=['1']
    assert assess(d,x,m)['H1']['verdict']=='UNINTERPRETABLE'
