"""Summarize fixed OLS4 ontology search audit from preregistration_ols_anatomy.md.

The live query snapshot is results/ols_anatomy.json; no network request is made on import.
"""
import json
from pathlib import Path

FIXED = {'uberon:kidney cortex':('cortex of kidney','UBERON:0001225'),
         'uberon:kidney':('kidney','UBERON:0002113'),
         'uberon:liver':('liver','UBERON:0002107'),
         'uberon:bile duct':('bile duct','UBERON:0002394'),
         'uberon:cerebral cortex':('cerebral cortex','UBERON:0000956'),
         'uberon:colon':('colon','UBERON:0001155'),
         'uberon:lung':('lung','UBERON:0002048'),
         'cl:cholangiocyte':('cholangiocyte','CL:1000488'),
         'cl:hepatocyte':('hepatocyte','CL:0000182')}


def assess(results):
    resolved = {}
    for query, (label, obo) in FIXED.items():
        matches = [m for m in results.get(query, {}).get('matches', [])
                   if m.get('label', '').lower() == label and m.get('obo_id') == obo]
        resolved[query] = {'url': results.get(query, {}).get('url'),
                           'obo_id': obo if len(matches) == 1 else None,
                           'matches_exact_label_and_id': len(matches)}
    primary = all(resolved[q]['matches_exact_label_and_id'] == 1 for q in
                  ('uberon:kidney cortex','uberon:kidney','cl:cholangiocyte','cl:hepatocyte'))
    return {'source': 'https://www.ebi.ac.uk/ols4/api/search',
            'preregistration': 'results/preregistration_ols_anatomy.md (commit fecf4f7)',
            'resolved': resolved,
            'primary_H1': 'PASS' if primary else 'UNRESOLVED',
            'proximal_tubule_exact_query': {'url': results['cl:proximal tubule cell']['url'],
              'n_returned':len(results['cl:proximal tubule cell']['matches']), 'verdict':'UNRESOLVED (no exact query result)'},
            'interpretation': 'Kidney cortex and kidney have distinct UBERON IDs. Cholangiocyte and hepatocyte have distinct CL IDs. Distinct labels do not prove a specific part-of hierarchy or molecular similarity. Comparing biliary organoids to whole liver GTEx and kidney organoids to adult kidney cortex requires cell-type/compartment-aware validation; ontology cannot rescue or refute the replication outcome.'}

if __name__ == '__main__':
    out = assess(json.loads(Path('results/ols_anatomy.json').read_text()))
    Path('results/ols_anatomy_audit.json').write_text(json.dumps(out, indent=2)+'\n')
    print(out['primary_H1'],out['proximal_tubule_exact_query']['verdict'])
