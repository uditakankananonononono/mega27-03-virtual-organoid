from pathlib import Path
import json
import subprocess
import sys
import pandas as pd


def test_record_is_multi_patient_and_preregistered():
    r=json.loads(Path('results/drevinek_fis.json').read_text())
    assert r['preregistration_commits']==['c953790','edd828b']
    assert r['patients_source']==20 and r['plates_source']==54
    assert len(r['patient_results'])==20
    assert r['blocks_eligible']+sum('fewer_than_two_eligible_wells_per_arm' in b['reasons'] for b in r['blocks'])==r['blocks_total']
    assert sum(x['n_blocks'] for x in r['patient_results'])==r['blocks_eligible']


def test_auc_and_sign_p():
    from scripts.analyze_drevinek_fis import sign_p
    assert sign_p([1,1,1,1])==1/16
    assert sign_p([1,-1,0])==.75
    assert sign_p([]) is None
    r=json.loads(Path('results/drevinek_fis.json').read_text())
    b=next(b for b in r['blocks'] if not b['reasons'])
    for arm in ('TEZ/IVA','ELX/TEZ/IVA'):
        assert len(b['eligible_wells'][arm])>=2
        for w in b['eligible_wells'][arm]:
            assert w['horizon']==60 and w['n_times']==7
    assert r['primary']['exact_one_sided_sign_p']==2**-20


def test_recompute_does_not_change_result():
    path=Path('results/drevinek_fis.json')
    before=path.read_bytes()
    subprocess.run([sys.executable,'scripts/analyze_drevinek_fis.py'],check=True,capture_output=True)
    assert path.read_bytes()==before
