"""Preregistered well-level multi-patient FIS comparison; see results/preregistration_drevinek_fis.md."""
import json
from pathlib import Path
from math import comb
import numpy as np
import pandas as pd

SRC = Path('data/drevinek_fis/area_data.csv')
OUT = Path('results/drevinek_fis.json')
ARMS = ('TEZ/IVA', 'ELX/TEZ/IVA')

def sign_p(effects):
    vals = [v for v in effects if v != 0]
    n, k = len(vals), sum(v > 0 for v in vals)
    return sum(comb(n, j) for j in range(k, n + 1)) / 2**n if n else None

def main():
    df = pd.read_csv(SRC)
    total = len(df)
    valid = np.isfinite(df.area) & (df.area > 0) & np.isfinite(df.time) & (df.time >= 0) & np.isfinite(df.fsk_concentration)
    for c in ['patient','date','filename','well','type','replicate']:
        valid &= df[c].notna()
    invalid_rows = int((~valid).sum()); df = df[valid].copy()
    patient = df[df.type.isin(ARMS)].copy()
    # Well is unique within source plate, but keep all plate map fields to audit ambiguities.
    groupcols = ['patient','date','filename','well','type','fsk_concentration','replicate']
    duplicate_key = groupcols + ['time']
    dup = patient.groupby(duplicate_key, dropna=False).area.agg(['size','nunique']).reset_index()
    conflicting = dup[(dup['size'] > 1) & (dup['nunique'] > 1)]
    ambiguous_keys = set(tuple(x) for x in conflicting[duplicate_key].itertuples(index=False,name=None))
    if ambiguous_keys:
        patient = patient[~patient[duplicate_key].apply(tuple,axis=1).isin(ambiguous_keys)]
    patient = patient.drop_duplicates(duplicate_key)
    well_map = {}
    for key,g in patient.groupby(groupcols, dropna=False):
        well_map[key] = {int(t):float(a) for t,a in zip(g.time,g.area)}
    blockcols = ['patient','date','filename','fsk_concentration']
    blocks = []
    for bkey,g in patient.groupby(blockcols):
        wells = {arm:[] for arm in ARMS}
        for key, times in well_map.items():
            if (key[0],key[1],key[2],key[5]) == bkey and key[4] in ARMS:
                wells[key[4]].append((key,times))
        # Intersection of recorded times across the two arms: present in each arm, with baseline.
        common = set.intersection(*(set().union(*(set(t) for _,t in wells[a])) for a in ARMS)) if all(wells.values()) else set()
        # For each well, only use actual common times where that well has measurements; require 0 and >=3 points.
        eligible = {a:[] for a in ARMS}
        reasons = []
        if 0 not in common or len(common)<3:
            reasons.append('fewer_than_three_shared_treatment_times_including_zero')
        else:
            for arm in ARMS:
                for key,t in wells[arm]:
                    points = sorted(common.intersection(t))
                    if 0 not in points or len(points)<3: continue
                    horizon = max(points)
                    if not horizon: continue
                    areas=np.array([t[p] for p in points]); times=np.array(points)
                    norm=float(np.trapezoid(areas/areas[0],times)/horizon)
                    raw=float(np.trapezoid(areas,times)/horizon)
                    eligible[arm].append({'well':int(key[3]),'replicate':str(key[6]),'n_times':len(points),'horizon':int(horizon),'baseline':float(areas[0]),'auc_normalized':norm,'auc_raw':raw})
            if any(len(eligible[a])<2 for a in ARMS):reasons.append('fewer_than_two_eligible_wells_per_arm')
        record={'patient':bkey[0],'date':bkey[1],'filename':bkey[2],'dose':float(bkey[3]),'available_wells':{a:len(wells[a]) for a in ARMS},'eligible_wells':eligible,'shared_times':sorted(common),'reasons':reasons}
        if not reasons:
            means={a:float(np.mean([w['auc_normalized'] for w in eligible[a]])) for a in ARMS}
            medians={a:float(np.median([w['auc_normalized'] for w in eligible[a]])) for a in ARMS}
            raws={a:float(np.mean([w['auc_raw'] for w in eligible[a]])) for a in ARMS}
            record.update({'mean_normalized':means,'effect':means[ARMS[1]]-means[ARMS[0]],'median_effect':medians[ARMS[1]]-medians[ARMS[0]],'raw_effect':raws[ARMS[1]]-raws[ARMS[0]]})
        blocks.append(record)
    eligible_blocks=[b for b in blocks if not b['reasons']]
    pats=[]
    for p,g in pd.DataFrame([{'patient':b['patient'],'effect':b['effect'],'median_effect':b['median_effect'],'raw_effect':b['raw_effect']} for b in eligible_blocks]).groupby('patient'):
        pats.append({'patient':p,'n_blocks':len(g),'effect':float(g.effect.median()),'median_well_effect':float(g.median_effect.median()),'raw_effect':float(g.raw_effect.median()),'doses':sorted(set(b['dose'] for b in eligible_blocks if b['patient']==p))})
    effects=[p['effect'] for p in pats]
    passed=len(pats)>=6 and sum(e>0 for e in effects)/len(effects)>=.75 and sign_p(effects)<.05
    leave_one=[{'omitted':p['patient'],'n':len(pats)-1,'positive':sum(q['effect']>0 for q in pats if q!=p),'p':sign_p([q['effect'] for q in pats if q!=p]),'median_effect':float(np.median([q['effect'] for q in pats if q!=p]))} for p in pats]
    reversals=[b for b in eligible_blocks if b['effect']*b['raw_effect']<0]
    result={'source':'https://zenodo.org/records/4771466','preregistration_commits':['c953790','edd828b'],'source_zip_sha256':'366e4da2e931e96051ad69ac3862a83cd44dff275c97b77e607de6c9622285ae','total_rows':total,'invalid_rows':invalid_rows,'conflicting_duplicate_keys':len(conflicting),'patients_source':int(df.patient.nunique()),'plates_source':int(df[['patient','date','filename']].drop_duplicates().shape[0]),'blocks_total':len(blocks),'blocks_eligible':len(eligible_blocks),'block_exclusions':pd.Series([r for b in blocks for r in b['reasons']]).value_counts().to_dict(),'patient_results':pats,'primary':{'eligible_patients':len(pats),'positive_patients':sum(e>0 for e in effects),'median_patient_effect':float(np.median(effects)),'exact_one_sided_sign_p':sign_p(effects),'H1_pass':passed},'sensitivity':{'median_well_positive_patients':sum(p['median_well_effect']>0 for p in pats),'raw_area_positive_patients':sum(p['raw_effect']>0 for p in pats),'raw_vs_normalized_block_reversals':len(reversals),'leave_one_patient_out':leave_one},'blocks':blocks}
    OUT.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:result[k] for k in ['total_rows','invalid_rows','patients_source','plates_source','blocks_total','blocks_eligible','block_exclusions','primary','sensitivity'] if k!='sensitivity'},indent=2))
    print('sensitivities', {k:v for k,v in result['sensitivity'].items() if k!='leave_one_patient_out'})
if __name__=='__main__':main()
