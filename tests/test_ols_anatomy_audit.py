import sys
sys.path.insert(0, 'src')
from ols_anatomy_audit import FIXED, assess

def test_missing_crucial_term_unresolved():
    r = {q:{'url':'https://example.invalid', 'matches':[{'label':label,'obo_id':obo}]}
         for q,(label,obo) in FIXED.items()}
    r['cl:proximal tubule cell'] = {'url':'https://example.invalid','matches':[]}
    assert assess(r)['primary_H1'] == 'PASS'
    r['cl:cholangiocyte']['matches'] = []
    assert assess(r)['primary_H1'] == 'UNRESOLVED'
