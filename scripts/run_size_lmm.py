"""VO3-R2 + amendment A1: donor-level size-dependence LMM + exact donor sign-flip test.
Blocks/cells byte-identical to src/blocked_size.py assess(); LMM on block-cell mean log swelling."""
import json, hashlib, sys, itertools
import numpy as np, pandas as pd
sys.path.insert(0, 'src')
from blocked_size import assess
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import chi2

df = pd.read_csv('data/raw/orgasegment/dis_merged_A0.csv')
base = assess(df)  # identical filters/eligibility to locked analysis
blocks = base['blocks']
keys = {(b['donor'], b['experiment'], str(b['forskolin_concentration_µM'])) for b in blocks}

y = df[(df['forskolin_concentration_µM'] > 0) & df.condition.isin(['DMSO', 'VX445_VX661_VX770'])
       & (df.swelling > 0) & np.isfinite(df.swelling) & (df.A0 > 0) & np.isfinite(df.A0)].copy()
y['size'] = np.select([y.A0 < 722, y.A0 >= 1400], ['small', 'large'], default='middle')
y = y[y['size'] != 'middle'].copy()
y['log_swelling'] = np.log(y.swelling)
y['_k'] = list(zip(y['donor'].astype(str), y['experiment'].astype(str), y['forskolin_concentration_µM'].astype(str)))
y = y[y['_k'].isin(keys)].drop(columns='_k')
cm = (y.groupby(['donor', 'experiment', 'forskolin_concentration_µM', 'condition', 'size'])
        .log_swelling.mean().reset_index())
cm['treatment'] = (cm.condition == 'VX445_VX661_VX770').astype(int)
cm['large'] = (cm['size'] == 'large').astype(int)
assert len(cm) == 4 * len(blocks), (len(cm), len(blocks))

def fit(form):
    md = smf.mixedlm(form, cm, groups=cm['donor'],
                     vc_formula={'experiment': '0 + C(experiment)'})
    return md.fit(reml=False, method='lbfgs')

full = fit('log_swelling ~ treatment * large')
red = fit('log_swelling ~ treatment + large')
stat = 2 * (full.llf - red.llf)
beta = full.params.get('treatment:large', np.nan)
p2 = chi2.sf(stat, 1)
p1 = p2 / 2 if beta > 0 else 1 - p2 / 2

# exact donor sign-flip on donor median block effects (2^12 enumeration)
dm = pd.DataFrame(blocks).groupby('donor').effect.median()
vals = dm.values; obs = vals.mean()
stats = np.array([np.mean(vals * np.array(s)) for s in itertools.product([1, -1], repeat=len(vals))])
p_flip = float((stats >= obs - 1e-12).mean())

lrt_verdict = 'DISCOVERY' if (p1 < 0.05 and p_flip < 0.05) else \
              'INCONCLUSIVE' if (0.05 <= p1 < 0.10 or 0.05 <= p_flip < 0.10) else 'NEGATIVE'
out = {'n_blocks': len(blocks), 'n_donors': int(dm.size), 'n_cell_rows': len(cm),
       'lmm': {'formula': 'log_swelling ~ treatment * large + (1|donor) + (1|experiment)',
               'fit': 'ML, lbfgs', 'interaction_beta': float(beta), 'lrt_stat': float(stat),
               'lrt_p_two_sided': float(p2), 'lrt_p_one_sided_positive': float(p1),
               'converged_full': bool(full.converged), 'converged_reduced': bool(red.converged)},
       'exact_donor_signflip': {'n_donors': int(dm.size), 'observed_mean_median_effect': float(obs),
                                'n_assignments': int(len(stats)), 'p_one_sided_positive': p_flip},
       'verdict': lrt_verdict,
       'rule': 'discovery requires LRT p<0.05 AND sign-flip p<0.05; [0.05,0.10] => INCONCLUSIVE',
       'input_sha256': {'data/raw/orgasegment/dis_merged_A0.csv':
                        hashlib.sha256(open('data/raw/orgasegment/dis_merged_A0.csv', 'rb').read()).hexdigest()}}
json.dump(out, open('results/size_lmm.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
