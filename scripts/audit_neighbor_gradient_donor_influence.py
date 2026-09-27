"""Reconcile and influence-check already scored within-image/arm/size gradient."""
import hashlib,json,statistics
from pathlib import Path
src=Path('data/raw/orgasegment/ResearchData/DIS_database.csv')
x=json.load(open('results/postresult_neighbor_stratified.json'))
assert hashlib.sha256(src.read_bytes()).hexdigest()==x['source_sha256']=='0e0e7913947a29c7b32ce87ca73a98614033d41ce63063a1d65053e766b2f315'
assert x['n_eligible_strata']+x['n_excluded_strata']==284
assert x['n_eligible_objects']+x['n_excluded_objects']==x['n_frame']==6964
assert x['n_positive_A1']==4964 and x['donor_n']==12 and x['donor_negative_count']==12
vals=x['donor_differences'];assert len(vals)==12
assert abs(statistics.mean(vals.values())-x['donor_equal_close_minus_far'])<1e-12
loo={k:statistics.mean(v for d,v in vals.items() if d!=k) for k in vals}
out={'status':'post-result accounting of known gradient, no additional scientific validation','protocol':'notes/postresult_neighbor_gradient_decomposition_protocol.md','source_url':x['source_url'],'source_sha256':x['source_sha256'],'n_source_baseline_detected':x['n_frame'],'n_eligible_strata':x['n_eligible_strata'],'n_excluded_strata':x['n_excluded_strata'],'n_eligible_objects':x['n_eligible_objects'],'n_excluded_objects':x['n_excluded_objects'],'n_positive_A1':x['n_positive_A1'],'donor_equal_close_minus_far':x['donor_equal_close_minus_far'],'donors_sorted':[{'donor':k,'close_minus_far':v,'leave_this_donor_out_equal_mean':loo[k]} for k,v in sorted(vals.items(),key=lambda z:z[1])],'donor_median':statistics.median(vals.values()),'donor_min':min(vals.values()),'donor_max':max(vals.values()),'leave_one_out_equal_mean_range':[min(loo.values()),max(loo.values())],'limits':'Previously seen source association; no post-hoc p-value or causal crowding claim. Absent endpoints and baseline-undetected objects remain unobserved. Original sign/specificity FAIL, benchmark win withdrawn.'}
Path('results/postresult_neighbor_gradient_donor_influence.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('strata',out['n_eligible_strata'],out['n_excluded_strata'],'objects',out['n_eligible_objects'],out['n_excluded_objects'],'donor median/range',out['donor_median'],out['donor_min'],out['donor_max'],'LOO',out['leave_one_out_equal_mean_range'])
