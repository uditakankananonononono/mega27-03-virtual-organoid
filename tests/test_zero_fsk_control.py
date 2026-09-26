"""The frozen 0 vs .128 uM contrast is reproducible from source rows."""
import json
import subprocess


def test_recompute_negative_without_changing_verdict():
    subprocess.run(['python3','scripts/analyze_zero_fsk_control.py'],check=True,capture_output=True,text=True)
    r=json.load(open('results/zero_fsk_control.json'))
    assert (r['donor_experiments_source'],r['paired_experiments'],r['eligible_donors'],
            r['positive_paired_delta_donors'],r['verdict'])==(50,24,11,3,'FAIL')
    assert abs(r['exact_one_sided_sign_p']-0.96728515625)<1e-12
    assert all(abs(e['delta']-(e['effect_0128']-e['effect_0']))<1e-12 for e in r['experiments'])
