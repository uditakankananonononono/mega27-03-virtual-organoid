import sys
sys.path.insert(0,'src')
import pandas as pd
from external_fis_missingtracks import assess

def test_unobserved_endpoint_audit():
    t=pd.DataFrame([{'Metadata_compound':'fsk','Metadata_concentration':.1,'Metadata_timeNum':time,
     'Metadata_wellNum':1,'TrackObjects_Label_4':track,'Math_area_micronsq':area} for time,track,area in [(0,1,1000),(0,2,9000),(6,1,1200)]])
    x=assess([('W0001--fsk--0.1/objects.csv',t)])
    assert x['M1']['fsk']['large']['paired']==0
    assert x['M2']['verdict']=='UNINTERPRETABLE'
