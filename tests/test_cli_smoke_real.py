"""Read-only real-data smoke tests of the shipped commands (no download)."""
import json
from pathlib import Path
import pandas as pd
from vorganoid.cli import main


def test_real_sizeaware_smoke(capsys):
    assert main(['sizeaware','data/raw/orgasegment/dis_merged_A0.csv',
                 '--drug','VX445_VX661_VX770','--forskolin-col','forskolin_concentration_µM'])==0
    out=capsys.readouterr().out
    d=json.loads(out[out.index('{'):])
    assert d['n_donors']==14 and .24<d['mean_attenuation']<.26


def test_real_fidelity_profile_smoke(capsys):
    f=Path('data/geo/profiles/GSE186249.csv.gz')
    if not f.exists(): return  # large GEO profiles are intentionally not in a fresh clone
    assert main(['fidelity',str(f),'--profile','--reference','data/ref/gtex_v8_median_tpm.gct.gz',
                 '--organ-tissue','Liver','--purity','data/ref/hpa_single_cell_type_v_api.tsv.gz'])==0
    d=json.loads(capsys.readouterr().out)
    assert d['organ_tissue_rank']==1 and d['top']['Liver']>.65
