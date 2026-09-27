"""Source-frame denominator and analysis-scope checks, no outcome imputation."""
import hashlib,json
from pathlib import Path
import pandas as pd
s=Path('data/raw/orgasegment/ResearchData/DIS_database.csv')
sha=hashlib.sha256(s.read_bytes()).hexdigest();assert sha=='0e0e7913947a29c7b32ce87ca73a98614033d41ce63063a1d65053e766b2f315'
source=pd.read_csv(s,nrows=2)
ret=json.load(open('results/tabular_track_retention.json'));mb=json.load(open('results/missing_endpoint_bounds.json'));cf=json.load(open('results/crossfit_inclusion_sensitivity.json'))
assert mb['source_sha256']==cf['source_sha256']==ret['source_sha256']==sha
cells={r['arm']+':'+r['size_bin']:r for r in ret['retention_by_arm_size'] if r['size_bin'] in ['small','large']}
ref=mb['percentile_reference'];assert set(cells)==set(ref)
rows=[]
for k,r in cells.items():
 z=ref[k];assert z['n_observed']==r['n_both_positive_area'] and z['n_missing']==r['n_baseline_detected']-r['n_both_positive_area']
 rows.append({'arm_size':k,'baseline_detected':r['n_baseline_detected'],'positive_A1':r['n_both_positive_area'],'no_positive_A1':z['n_missing'],'positive_A1_rate':r['n_both_positive_area']/r['n_baseline_detected']})
assert sum(r['baseline_detected'] for r in rows)==mb['n_baseline_detected']==6964
assert sum(r['positive_A1'] for r in rows)==cf['n_positive_source_pairs_in_original_small_large_cells']==4964
assert sum(r['no_positive_A1'] for r in rows)==2000
assert cf['n_derived_included_in_frame']<4964
out={'status':'post-result source-denominator and inclusion-scope audit, not missing-outcome correction',
 'protocol':'notes/postresult_attrition_model_limit_audit.md','source_url':'https://zenodo.org/records/10610438','source_sha256':sha,
 'available_source_columns':list(source.columns),'cells':rows,'n_baseline_detected':6964,'n_positive_A1':4964,
 'n_without_positive_A1':2000,'n_derived_included_given_positive_A1':cf['n_derived_included_in_frame'],
 'n_positive_pairs_excluded_from_derived':4964-cf['n_derived_included_in_frame'],
 'crossfit_uses_endpoint_features':['log_A1','score1'], 'crossfit_target':'derived inclusion only among positive A0 and A1 source pairs',
 'crossfit_cannot_weight_missing_A1_count':2000,'hypothetical_response_scenarios':{'adverse_equal_donor_mean':mb['minimize_contrast']['equal_donor_mean_median'],'favorable_equal_donor_mean':mb['maximize_contrast']['equal_donor_mean_median'],'zero_contained':mb['zero_contained']},
 'unavailable_from_tabular_source':['raw t0/t1 image pixels and masks for these source rows','reasons for absent or nonpositive A1','objects never detected at t0'],
 'limits':'Denominators conditional on source-detected positive A0 in the fixed 31 blocks/4 cells. Missing response scenario endpoints are assumed from observed 5th/95th percentiles, not empirical support or confidence bounds. Crossfit includes endpoint features and does not identify a full-population causal effect. Sign-only/specificity gates remain failed.'}
Path('results/postresult_attrition_model_limit.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:out[k] for k in ['n_baseline_detected','n_positive_A1','n_without_positive_A1','n_derived_included_given_positive_A1','n_positive_pairs_excluded_from_derived']})
print(rows)
