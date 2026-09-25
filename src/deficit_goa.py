"""Direct EBI GOA human BP annotation of top deficits, preregistration_goa.md."""
import json
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact
PATH='data/goa/goa_human.gaf.gz'
URL='https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz'
TERMS={'GO:0006631':'fatty acid metabolic process','GO:0006805':'xenobiotic metabolic process','GO:0006699':'bile acid biosynthetic process'}
COLS=['DB','DB_Object_ID','DB_Object_Symbol','Qualifier','GO_ID','DB_Reference','Evidence_Code','With_From','Aspect','DB_Object_Name','DB_Object_Synonym','DB_Object_Type','Taxon','Date','Assigned_By','Annotation_Extension','Gene_Product_Form_ID']

def assess(deficits,gaf):
    a=gaf[(gaf.Aspect=='P') & (~gaf.Qualifier.fillna('').str.contains('NOT'))]
    all_annotated=set(a.DB_Object_Symbol.dropna())
    selected=a[a.GO_ID.isin(TERMS)]
    genes=set(selected.DB_Object_Symbol.dropna())
    out={'source':URL,'preregistration':'results/preregistration_goa.md (commit fc39e36)',
         'terms':TERMS,'n_selected_rows':len(selected),'n_selected_symbols':len(genes),
         'evidence_codes_selected':selected.Evidence_Code.value_counts().to_dict(),'per_organ':{}}
    for o,g in deficits.groupby('organ'):
        s=set(g.symbol.dropna());out['per_organ'][o]={'n_input':len(s),'n_any_BP_annotation':len(s&all_annotated),
            'n_direct_metabolic_union':len(s&genes),'direct_metabolic_symbols':sorted(s&genes)}
    L,B=out['per_organ']['Liver'],out['per_organ']['Brain - Cortex']
    if min(L['n_any_BP_annotation'],B['n_any_BP_annotation'])<80:
        out['H1']={'verdict':'UNINTERPRETABLE','reason':'<80 genes with any BP annotation in one list'}
    else:
        tab=[[L['n_direct_metabolic_union'],L['n_input']-L['n_direct_metabolic_union']],
             [B['n_direct_metabolic_union'],B['n_input']-B['n_direct_metabolic_union']]]
        OR,p=fisher_exact(tab,alternative='greater')
        out['H1']={'table':tab,'odds_ratio':OR,'p_one_sided':p,'verdict':'PASS' if p<.05 else 'FAIL'}
    return out

if __name__=='__main__':
    d=pd.read_csv('results/strict_organ_deficits.csv')
    a=pd.read_csv(PATH,sep='\t',comment='!',header=None,names=COLS,dtype=str)
    out=assess(d,a);Path('results/deficit_goa.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['H1']);print({k:(v['n_any_BP_annotation'],v['n_direct_metabolic_union'],v['direct_metabolic_symbols']) for k,v in out['per_organ'].items()})
