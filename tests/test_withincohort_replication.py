import sys
sys.path.insert(0,'src')
import pandas as pd
from withincohort_replication import date_key,assess

def test_date_and_min_cell_floor():
    assert date_key('HitCF_A1_DIS_20230207_plate001')=='20230207'
    rows=[]
    for exp in ('20220911_p1','20221011_p2'):
        for condition in ('DMSO','VX445_VX661_VX770'):
            for area in (500,2000):
                rows += [dict(experiment=exp,donor='D',condition=condition,A0=area,swelling=2,
                               **{'forskolin_concentration_µM':0.128})]*5
    x=assess(pd.DataFrame(rows));assert x['portions']['earliest']['n_donors']==1
    assert x['portions']['later']['n_donors']==1
    assert x['H1']['verdict']=='UNINTERPRETABLE'
