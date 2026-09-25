import sys
sys.path.insert(0,'src')
import pandas as pd
from external_fis_robustness import assess

def test_well_floor_and_retention():
    rows=[]
    for t in [0,6]:
        for i,area in [(1,1000),(2,9000)]:
            rows.append({'Metadata_compound':'fsk','Metadata_concentration':.1,'Metadata_timeNum':t,
                         'Metadata_wellNum':1,'TrackObjects_Label_4':i,'Math_area_micronsq':area})
    a=assess([('W0001--fsk--0.1/objects.csv',pd.DataFrame(rows))],boots=10)
    assert a['R1']['verdict']=='UNINTERPRETABLE'
    assert a['R2']['by_condition_size']['fsk']['small']['rate']==1
