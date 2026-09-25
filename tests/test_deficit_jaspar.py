import sys
sys.path.insert(0,'src')
from deficit_jaspar import pick,assess,NAMES

def test_exact_name_and_latest_version():
    data=[{'name':'NEUROD2-like','matrix_id':'MA1.1','version':'1'},
          {'name':'NEUROD2','matrix_id':'MA2.1','version':'1'},
          {'name':'Neurod2','matrix_id':'MA2.3','version':'3'}]
    assert pick('NEUROD2',data)['matrix_id']=='MA2.3'
    x={k:None for k in NAMES}; assert assess(x)['H1']['verdict']=='FAIL'
