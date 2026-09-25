import sys,json
sys.path.insert(0,'src')
from deficit_interpro import classify

def test_domain_keyword_classification_and_invalid_hypothesis():
    assert classify(['Cytochrome P450']) == ['cytochrome p450']
    assert classify(['Globin-like superfamily']) == []
    j=json.load(open('results/deficit_interpro.json'))
    assert j['deficit_top100_membership']['liver']==['CYP2D6']
    assert j['H1']['verdict'].startswith('INVALID')
