import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'results'
def test_stratified_neighbor_frame():
 x=json.loads((ROOT/'postresult_neighbor_stratified.json').read_text());p=json.loads((ROOT/'postresult_neighbor_observability.json').read_text())
 assert x['source_sha256']==p['source_sha256']
 assert (x['n_frame'],x['n_positive_A1'],x['n_images_frame'])==(6964,4964,143)
 assert (x['n_eligible_strata'],x['n_excluded_strata'])==(279,5)
 assert x['n_eligible_objects']+x['n_excluded_objects']==x['n_frame']
 assert sum(r['close_n']+r['far_n'] for r in x['eligible_strata'])==x['n_eligible_objects']
 assert sum(r['n'] for r in x['excluded_strata'])==x['n_excluded_objects']
 assert all(r['close_n']>0 and r['far_n']>0 for r in x['eligible_strata'])
 assert sum(r['n_strata'] for r in x['per_arm_size'])==279
 assert x['donor_n']==len(x['donor_differences'])==12
 assert x['permutation_n']==100
 assert x['equal_stratum_close_minus_far']<0
