import json
from pathlib import Path
from math import exp
ROOT=Path(__file__).resolve().parents[1]/'results'
def test_neighbor_observability_full_frame_and_caveats():
    x=json.loads((ROOT/'postresult_neighbor_observability.json').read_text())
    p=json.loads((ROOT/'endpoint_observability.json').read_text())
    assert x['source_sha256']==p['source_sha256']
    assert (x['n_original_blocks'],x['n_donors'])==(31,12)
    assert (x['n_eligible_baseline'],x['n_observed_endpoint'])==(6964,4964)
    assert x['n_missing_nn1']==x['n_missing_nn3']==0
    assert x['n_eligible_fields']==143
    assert len(x['quartile_rates'])==16
    assert sum(z['n'] for z in x['quartile_rates'])==6964
    assert sum(z['observed'] for z in x['quartile_rates'])==4964
    assert all(0<=z['rate']<=1 for z in x['quartile_rates'])
    assert x['model_converged']
    assert abs(exp(x['model_terms']['log_nn1']['beta'])-x['model_terms']['log_nn1']['odds_ratio'])<1e-12
    assert x['within_field_close_minus_far']['n_usable_fields']==143
    assert x['within_field_close_minus_far']['n_permutations']==100
    assert x['within_field_close_minus_far']['observed_mean_field_difference']<0
