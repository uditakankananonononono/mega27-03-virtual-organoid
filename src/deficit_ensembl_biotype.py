"""Ensembl REST coding-composition check, preregistration_ensembl_biotype.md."""
import json
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact

URL = 'https://rest.ensembl.org/lookup/id'


def analyze(deficits, annotations):
    per = {}
    for organ, grp in deficits.groupby('organ'):
        ids = grp.ens.dropna().drop_duplicates().tolist()
        matched = [annotations[i] for i in ids if annotations.get(i) and annotations[i].get('biotype')]
        n = len(matched)
        coding = sum(x['biotype'] == 'protein_coding' for x in matched)
        per[organ] = {'n_input': len(ids), 'n_mapped': n, 'n_protein_coding': coding,
                      'coding_fraction': coding/n if n else None}
    l, b = per['Liver'], per['Brain - Cortex']
    result = {'source': URL, 'documentation': 'https://rest.ensembl.org/documentation/info/lookup_post',
              'preregistration': 'results/preregistration_ensembl_biotype.md (commit 8ea5c50)',
              'assembly': sorted(set(v.get('assembly_name') for v in annotations.values() if v and v.get('assembly_name'))),
              'per_organ': per}
    if min(l['n_mapped'], b['n_mapped']) < 80:
        result['H1'] = {'verdict': 'UNINTERPRETABLE', 'reason': 'mapping <80 in one organ'}
    else:
        table = [[l['n_protein_coding'], l['n_mapped']-l['n_protein_coding']],
                 [b['n_protein_coding'], b['n_mapped']-b['n_protein_coding']]]
        odds, p = fisher_exact(table, alternative='greater')
        result['H1'] = {'table': table, 'odds_ratio': odds, 'p_one_sided': p,
                        'verdict': 'PASS' if p < 0.05 else 'FAIL'}
    return result


if __name__ == '__main__':
    d = pd.read_csv('results/strict_organ_deficits.csv')
    annotations = json.loads(Path('data/ref/ensembl_biotypes_lookup.json').read_text())
    result = analyze(d, annotations)
    Path('results/deficit_ensembl_biotype.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
