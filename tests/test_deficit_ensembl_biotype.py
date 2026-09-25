import sys
import pandas as pd
sys.path.insert(0, 'src')
from deficit_ensembl_biotype import analyze


def test_unmapped_does_not_enter_denominator():
    d = pd.DataFrame({'organ':['Liver','Liver','Brain - Cortex','Brain - Cortex'], 'ens':['L1','L2','B1','B2']})
    a = {'L1':{'biotype':'protein_coding','assembly_name':'GRCh38'},'B1':{'biotype':'lncRNA','assembly_name':'GRCh38'}}
    r = analyze(d,a)
    assert r['per_organ']['Liver']['n_mapped'] == 1
    assert r['per_organ']['Liver']['coding_fraction'] == 1
    assert r['H1']['verdict'] == 'UNINTERPRETABLE'
