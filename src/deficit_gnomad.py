"""Pre-registered population-constraint comparison, see results/preregistration_gnomad.md.

Uses gnomAD v2.1.1 gene-level pLoF metrics (LOEUF, oe_lof_upper).
"""
import gzip
import json
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu

SOURCE = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz'
PATH = Path('data/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')


def unique_gene_constraint(table):
    """Keep only non-ambiguous, finite gene-level constraint scores."""
    table = table.loc[~table.gene.duplicated(keep=False), ['gene', 'oe_lof_upper']].copy()
    table['oe_lof_upper'] = pd.to_numeric(table.oe_lof_upper, errors='coerce')
    return table.dropna().set_index('gene').oe_lof_upper


def analyze(deficits, constraint):
    organ_scores = {}
    for organ, grp in deficits.groupby('organ'):
        # The held-fixed input has 100 genes per organ. Avoid duplicating a gene.
        syms = grp.symbol.dropna().drop_duplicates()
        organ_scores[organ] = constraint.reindex(syms).dropna()
    L, B = organ_scores['Liver'], organ_scores['Brain - Cortex']
    out = {
        'source': SOURCE,
        'preregistration': 'results/preregistration_gnomad.md (commit 58a422c)',
        'mapping': 'unique gene symbol; ambiguous gnomAD gene symbols and missing LOEUF excluded',
        'per_organ': {o: {'n_mapped': len(v), 'median_loeuf': float(v.median())} for o, v in organ_scores.items()},
    }
    if min(len(L), len(B)) < 30:
        out['H1'] = {'verdict': 'UNINTERPRETABLE', 'reason': 'Fewer than 30 mapped genes in liver or brain'}
    else:
        U, p = mannwhitneyu(L, B, alternative='greater')
        out['H1'] = {'U': float(U), 'p_one_sided': float(p), 'verdict': 'PASS' if p < 0.05 else 'FAIL',
                     'direction': 'liver LOEUF greater than brain'}
    K = organ_scores['Kidney - Cortex']
    U, p = mannwhitneyu(K, B, alternative='greater')
    out['kidney_vs_brain_descriptive'] = {'U': float(U), 'p_one_sided': float(p)}
    return out


if __name__ == '__main__':
    with gzip.open(PATH, 'rt') as f:
        g = pd.read_csv(f, sep='\t', usecols=['gene', 'oe_lof_upper'])
    d = pd.read_csv('results/strict_organ_deficits.csv')
    out = analyze(d, unique_gene_constraint(g))
    Path('results/deficit_gnomad.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))
