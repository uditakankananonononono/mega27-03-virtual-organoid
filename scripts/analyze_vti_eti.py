"""Preregistered VTI versus ETI comparison; see results/preregistration_vti_eti.md."""
import json
from pathlib import Path
from math import comb
import numpy as np
import pandas as pd

SRC = Path('data/vti_eti/organoids_data.csv')
OUT = Path('results/vti_eti.json')
ARMS = ('ETI', 'VTI')

def sign_p(effects):
    vals = [v for v in effects if v != 0]
    n, k = len(vals), sum(v > 0 for v in vals)
    return sum(comb(n, j) for j in range(k, n + 1)) / 2**n if n else None

def main():
    df = pd.read_csv(SRC, dtype={'patient':str})
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
        record={'patient':str(bkey[0]),'date':str(bkey[1]),'filename':str(bkey[2]),'dose':float(bkey[3]),'available_wells':{a:len(wells[a]) for a in ARMS},'eligible_wells':eligible,'shared_times':sorted(common),'reasons':reasons}
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
    # Dose-optimization plates are a distinct descriptive subset, never substituted into H1.
    opt = df[df.type.isin(('ETI','VTI_002','VTI_02','VTI_2'))].copy()
    dose_cols = ['patient','date','filename','fsk_concentration','type','well']
    opt_auc = []
    for key, g in opt.groupby(dose_cols):
        g=g.sort_values('time')
        if len(g)!=7 or list(g.time)!=[0,10,20,30,40,50,60]: continue
        opt_auc.append((*key,float(np.trapezoid(g.area.to_numpy()/g.area.iloc[0], g.time.to_numpy())/60)))
    o = pd.DataFrame(opt_auc, columns=dose_cols+['auc'])
    o = o.groupby(dose_cols[:-1]).auc.mean().unstack('type')
    o = o.dropna(subset=['VTI_002','VTI_02','VTI_2'])
    doseopt = {'complete_plate_dose_blocks':len(o),'patients':sorted(set(str(v) for v in o.index.get_level_values('patient'))),
               'monotone_002_lt_02_lt_2':int(((o.VTI_002<o.VTI_02)&(o.VTI_02<o.VTI_2)).sum()),
               'variant_minus_ETI_block_medians':{k:float((o[k]-o.ETI).median()) for k in ['VTI_002','VTI_02','VTI_2']},
               'patient_VTI_02_minus_ETI_medians':{str(k):float(v) for k,v in (o.VTI_02-o.ETI).groupby(level='patient').median().items()}}
    result={'source':'https://zenodo.org/records/15754800','preregistration_commits':['08c5e22'],'source_csv_sha256':'46dbbf719d9c4fe4d195508fd81446f4e6e0e03e5f79bc45735851b0fae1ffc4','total_rows':total,'invalid_rows':invalid_rows,'conflicting_duplicate_keys':len(conflicting),'patients_source':int(df.patient.nunique()),'plates_source':int(df[['patient','date','filename']].drop_duplicates().shape[0]),'blocks_total':len(blocks),'blocks_eligible':len(eligible_blocks),'block_exclusions':pd.Series([r for b in blocks for r in b['reasons']]).value_counts().to_dict(),'patient_results':pats,'primary':{'eligible_patients':len(pats),'positive_patients':sum(e>0 for e in effects),'median_patient_effect':float(np.median(effects)),'exact_one_sided_sign_p':sign_p(effects),'H1_pass':passed},'sensitivity':{'median_well_positive_patients':sum(p['median_well_effect']>0 for p in pats),'raw_area_positive_patients':sum(p['raw_effect']>0 for p in pats),'raw_vs_normalized_block_reversals':len(reversals),'leave_one_patient_out':leave_one},'dose_optimization_descriptive':doseopt,'blocks':blocks}
    OUT.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:result[k] for k in ['total_rows','invalid_rows','patients_source','plates_source','blocks_total','blocks_eligible','block_exclusions','primary','sensitivity'] if k!='sensitivity'},indent=2))
    print('sensitivities', {k:v for k,v in result['sensitivity'].items() if k!='leave_one_patient_out'})
if __name__=='__main__':main()
