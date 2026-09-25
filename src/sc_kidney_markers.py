"""Real GSE108291 organoid 10x matrices, preregistration_sc_kidney.md.

Gene x barcode sparse matrix. Cell QC is >=500 total UMIs and >=200 detected genes;
no mitochondrial exclusion (not fixed in preregistration, so report as a deviation).
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import isspmatrix
from statsmodels.stats.proportion import proportion_confint
import scanpy as sc
from importlib.metadata import version

ROOT=Path('data/sc_kidney')
TARGETS=('ALDOB','SLC17A3')
URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108291/suppl/'


def summarize(X, genes, min_umis=500, min_genes=200):
    """Never densify the full expression matrix."""
    if X.shape[0] != len(genes): raise ValueError('gene dimension mismatch')
    if not isspmatrix(X): raise TypeError('sparse expression required')
    X=X.tocsr()
    idx={g:int(np.flatnonzero(genes.symbol.to_numpy()==g)[0]) for g in TARGETS
         if (genes.symbol==g).sum()==1}
    if len(idx)!=2: return {'verdict':'UNINTERPRETABLE','reason':'markers not uniquely mapped','marker_indices':idx}
    totals=np.asarray(X.sum(axis=0)).ravel(); detected=np.asarray(X.getnnz(axis=0)).ravel()
    candidates=np.flatnonzero(totals>=min_umis)
    # Scanpy computes QC metrics on the retained sparse cells; retain only the
    # small subset rather than materializing 2.2 million empty droplets.
    adata=sc.AnnData(X[:,candidates].T.tocsr())
    qc=sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=False)[0]
    if not np.array_equal(qc['total_counts'].to_numpy(),totals[candidates]):
        raise AssertionError('Scanpy QC totals disagree with sparse sums')
    if not np.array_equal(qc['n_genes_by_counts'].to_numpy(),detected[candidates]):
        raise AssertionError('Scanpy QC gene detection disagrees with sparse counts')
    mask=(totals>=min_umis)&(detected>=min_genes); n=int(mask.sum())
    a=np.asarray(X[idx['ALDOB'],:].toarray()).ravel()>0
    b=np.asarray(X[idx['SLC17A3'],:].toarray()).ravel()>0
    if n<100: return {'verdict':'UNINTERPRETABLE','reason':'<100 QC cells','n_qc':n}
    both=int((a&b&mask).sum()); either=int(((a|b)&mask).sum()); n_a=int((a&mask).sum()); n_b=int((b&mask).sum())
    lo,hi=proportion_confint(both,n,method='wilson')
    return {'verdict':'INTERPRETABLE','n_barcodes':X.shape[1], 'n_qc':n,
            'qc_min_umis':min_umis,'qc_min_genes':min_genes,
            'n_aldob':n_a,'n_slc17a3':n_b,'n_either':either,'n_both':both,
            'fraction_both':both/n,'fraction_either':either/n,'wilson95_both':[lo,hi],
            'H1_ge_1_percent':'PASS' if both/n>=0.01 else 'FAIL',
            'H2_lt_20_percent':'PASS' if both/n<0.20 else 'FAIL'}


if __name__=='__main__':
    import scanpy as sc
    runs={}
    for run in ('org4','org'):
        p=ROOT/run
        genes=pd.read_csv(str(p)+'_genes.tsv.gz',sep='\t',header=None,names=['ens','symbol'])
        print('Reading',run,flush=True)
        X=mmread(str(p)+'_matrix.mtx.gz').tocsr()
        res=summarize(X,genes)
        res['matrix_url']=URL+f'GSE108291_{run}_matrix.mtx.gz'
        res['genes_url']=URL+f'GSE108291_{run}_genes.tsv.gz'
        runs[run]=res
        del X
        print(run, {k:v for k,v in res.items() if k.startswith('n_') or k.startswith('H') or k=='fraction_both'}, flush=True)
    out={'accession':'GSE108291','source_project':'https://explore.data.humancellatlas.org/projects/7b947aa2-43a7-4082-afff-222a3e3a4635',
         'preregistration':'results/preregistration_sc_kidney.md (commit 3330723)',
         'analysis_library':'scanpy '+version('scanpy'), 'qc_deviation':'QC thresholds 500 UMIs and 200 genes fixed after preregistration, before seeing coexpression; no mitochondrial filter',
         'runs':runs}
    Path('results/sc_kidney_markers.json').write_text(json.dumps(out,indent=2)+'\n')
