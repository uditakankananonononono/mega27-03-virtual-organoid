import sys
sys.path.insert(0,'src')
from external_fis_fiji import selected,PREFIX

def test_exact_condition_filter():
    p=PREFIX+'W0001--fsk--0.008/P001--fsk--0.008/objects.csv'
    q=PREFIX+'W0003--fsk_809--0.008/P001--fsk_809--0.008/objects.csv'
    assert selected([p,q])==[p]
