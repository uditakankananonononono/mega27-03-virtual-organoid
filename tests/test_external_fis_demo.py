import sys
sys.path.insert(0,'src')
import pandas as pd
from external_fis_demo import selected,tracks,assess,PREFIX

def test_path_filter_and_pairs():
    p=PREFIX+'W0001--fsk--0.008/P001--fsk--0.008/objects.csv'
    assert selected([p,PREFIX+'W0002--fsk_809--0.008/x/objects.csv'])==[p]
    a=pd.DataFrame([{'Metadata_compound':'fsk','Metadata_concentration':1,'Metadata_timeNum':time,
                     'Metadata_wellNum':1,'TrackObjects_Label_4':1,'Math_area_micronsq':area} for time,area in [(0,100),(6,120)]])
    z,audit=tracks([(p,a)]);assert len(z)==1 and abs(z.log_swelling.iloc[0]-.1823215568)<1e-8
    assert assess(z,audit,[p])['H1']['verdict']=='UNINTERPRETABLE'
