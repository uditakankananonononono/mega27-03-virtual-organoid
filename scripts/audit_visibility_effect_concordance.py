"""Post-result, same-accession donor concordance; not biological validation."""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr,spearmanr
v=json.load(open('results/postresult_donor_attrition_decomposition.json'))
y=json.load(open('results/tabular_track_retention.json'))
assert v['source_sha256']==y['source_sha256']=='0e0e7913947a29c7b32ce87ca73a98614033d41ce63063a1d65053e766b2f315'
a=v['donors'];b=y['source_paired_positive_area']['donor_medians'];assert set(a)==set(b) and len(a)==12
ids=sorted(a);X=np.array([a[g]['gap_interaction_large_minus_small'] for g in ids]);Y=np.array([b[g] for g in ids]);
def stat(x,y):return {'pearson_r':float(pearsonr(x,y).statistic),'spearman_rho':float(spearmanr(x,y).statistic)}
s=stat(X,Y);leave={g:stat(np.delete(X,i),np.delete(Y,i)) for i,g in enumerate(ids)}
out={'status':'exploratory post-result donor observability and observed-size-contrast concordance; not mechanism or gate','protocol':'notes/postresult_visibility_effect_concordance_protocol.md','source_url':v['source_url'],'source_sha256':v['source_sha256'],'n_donors':len(ids),'donors':[{'id':g,'baseline_detected':a[g]['n_baseline_detected'],'positive_A1':a[g]['n_positive_A1'],'visibility_gap_large_minus_small':float(X[i]),'observed_positive_pair_log_swelling_contrast':float(Y[i])} for i,g in enumerate(ids)],'associations':s,'leave_one_donor_out':leave,'leave_one_out_ranges':{k:[min(q[k] for q in leave.values()),max(q[k] for q in leave.values())] for k in s},'sign_agreement':int(sum(np.sign(X)==np.sign(Y))),'limits':'No inferential p-values: outcomes and the selection question were known, n=12, unequal donor counts and different within-donor aggregation rules. Neither positive nor null association identifies 2,000 absent swelling responses, selection causality or original biological gate.'}
Path('results/postresult_visibility_effect_concordance.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('n',len(ids),'associations',s,'leave-one ranges',out['leave_one_out_ranges'],'sign agreement',out['sign_agreement'])
for row in out['donors']:print(row)
