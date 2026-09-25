import sys
sys.path.insert(0,'src')
import pandas as pd
from blocked_size import assess

def test_matched_block_and_floor():
    rows=[]
    for c in ('DMSO','VX445_VX661_VX770'):
        for sz in (500,1500):
            for _ in range(3):
                rows.append({'donor':'D','experiment':'plate1','forskolin_concentration_µM':1.0,'condition':c,'A0':sz,'swelling':2 if c=='DMSO' or sz==500 else 4})
    x=assess(pd.DataFrame(rows));assert x['n_eligible_blocks']==1
    assert x['n_positive_donors']==1 and x['H1']['verdict']=='UNINTERPRETABLE'
