import sys
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
sys.path.insert(0,'src')
from sc_kidney_markers import summarize

def test_sparse_coexpression_and_qc():
    g=pd.DataFrame({'symbol':['ALDOB','SLC17A3','OTHER']})
    X=csr_matrix(np.vstack([np.tile([1,1,0],40), np.tile([1,0,1],40), np.full(120,20)]))
    r=summarize(X,g,min_umis=1,min_genes=1)
    assert r['n_qc']==120 and r['n_both']==40 and r['n_either']==120
    assert r['H1_ge_1_percent']=='PASS' and r['H2_lt_20_percent']=='FAIL'
    assert summarize(X,g,min_umis=100,min_genes=1)['verdict']=='UNINTERPRETABLE'
