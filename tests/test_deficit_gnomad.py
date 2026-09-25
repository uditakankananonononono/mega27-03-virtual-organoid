import sys
import pandas as pd
sys.path.insert(0, 'src')
from deficit_gnomad import unique_gene_constraint, analyze


def test_ambiguous_symbols_and_missing_scores_excluded():
    t = pd.DataFrame({'gene':['A','A','B','C'], 'oe_lof_upper':[0.1,0.2,None,0.5]})
    assert unique_gene_constraint(t).to_dict() == {'C': 0.5}


def test_small_mapping_uninterpretable():
    d = pd.DataFrame({'organ':['Liver','Brain - Cortex','Kidney - Cortex'], 'symbol':['A','B','C']})
    out = analyze(d, pd.Series({'A':1.,'B':0.2,'C':0.3}))
    assert out['H1']['verdict'] == 'UNINTERPRETABLE'
