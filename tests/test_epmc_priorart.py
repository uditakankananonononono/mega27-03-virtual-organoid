import sys
sys.path.insert(0,'src')
from epmc_priorart import parse

def test_parse_preserves_abstract_not_author():
    row={'id':'42','source':'MED','pmid':'42','title':'Test','abstractText':'Size matters','authorString':'Wrong'}
    out=parse({'resultList':{'result':[row]}},'fixed query')
    assert out[0]['abstract']=='Size matters' and out[0]['url']=='https://europepmc.org/article/MED/42'
