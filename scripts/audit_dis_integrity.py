"""Audit predeclared OrgaSegment identity/arithmetic checks; never alter source rows."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

source = Path('data/raw/orgasegment/dis_merged_A0.csv')
df = pd.read_csv(source)
key = ['donor','experiment','well','condition','forskolin_concentration_µM','particle']
metrics = {'source': 'https://zenodo.org/records/10610438',
           'preregistration': 'results/preregistration_dis_integrity_audit.md (commit 12b6575)',
           'rows': len(df), 'donors': int(df.donor.nunique()),
           'full_row_duplicates': int(df.duplicated().sum()),
           'duplicate_ID_rows': int(df.duplicated(['ID'],keep=False).sum()),
           'duplicate_diagnostic_key_rows': int(df.duplicated(key,keep=False).sum())}
metrics['columns'] = {}
for col in ('A0','A1','swelling'):
    v = pd.to_numeric(df[col], errors='coerce').to_numpy()
    metrics['columns'][col] = {'missing': int(np.isnan(v).sum()),
                                'nonfinite': int((~np.isfinite(v)).sum()),
                                'nonpositive': int((v<=0).sum())}
a0 = pd.to_numeric(df.A0, errors='coerce').to_numpy()
a1 = pd.to_numeric(df.A1, errors='coerce').to_numpy()
s = pd.to_numeric(df.swelling, errors='coerce').to_numpy()
with np.errstate(divide='ignore',invalid='ignore'):
    rel = np.abs(s-a1/a0)/np.maximum(np.abs(s),1e-12)
metrics['swelling_A1_over_A0_relative_error'] = {
    'max': float(np.nanmax(rel)), 'over_1e-6': int(np.sum(rel>1e-6)),
    'over_1e-3': int(np.sum(rel>1e-3)),
    'five_highest_source_IDs': df.iloc[np.argsort(rel)[-5:][::-1]].ID.tolist()}
metrics['by_condition_and_dose'] = [
    {'condition': str(cond), 'dose_uM': float(dose), 'rows': int(len(g)),
     'donors': int(g.donor.nunique())}
    for (cond,dose),g in df.groupby(['condition','forskolin_concentration_µM'])]
Path('results/dis_integrity.json').write_text(json.dumps(metrics,indent=2)+'\n')
print(json.dumps({k:metrics[k] for k in ('rows','donors','full_row_duplicates','duplicate_ID_rows','duplicate_diagnostic_key_rows','swelling_A1_over_A0_relative_error')},indent=2))
