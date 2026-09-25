from pathlib import Path
import json
import subprocess
import sys

def test_preregistered_vti_eti_negative():
    r=json.loads(Path('results/vti_eti.json').read_text())
    assert r['preregistration_commits']==['08c5e22']
    assert r['patients_source']==12 and r['plates_source']==33
    assert r['blocks_eligible']==208 and r['blocks_total']==264
    assert r['primary']['eligible_patients']==9
    assert r['primary']['positive_patients']==6
    assert r['primary']['exact_one_sided_sign_p']==.25390625
    assert r['primary']['H1_pass'] is False
    assert sum(p['n_blocks'] for p in r['patient_results'])==208
    assert r['dose_optimization_descriptive']['complete_plate_dose_blocks']==48

def test_vti_recompute_stable():
    path=Path('results/vti_eti.json')
    before=path.read_bytes()
    subprocess.run([sys.executable,'scripts/analyze_vti_eti.py'],check=True,capture_output=True)
    assert path.read_bytes()==before
