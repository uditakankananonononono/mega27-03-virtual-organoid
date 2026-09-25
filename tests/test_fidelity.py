import numpy as np, pandas as pd, pytest
from vorganoid.fidelity import mean_log_cpm, fidelity_scores

def _ref(n=800, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.lognormal(1, 2, size=(n, 4)), index=[f"ENSG{i:011d}" for i in range(n)], columns=list("ABCD"))

def test_recovers_source_tissue():
    ref = _ref(); rng = np.random.default_rng(1)
    counts = pd.DataFrame({f"s{j}": rng.poisson(ref["C"] * 200) for j in range(3)}, index=ref.index)
    sc = fidelity_scores(mean_log_cpm(counts), ref)
    assert sc.index[0] == "C" and sc["C"] > 0.8

def test_drops_non_library_columns():
    ref = _ref(); counts = pd.DataFrame({"s": np.full(800, 500.0), "log2fc": np.ones(800)}, index=ref.index)
    p = mean_log_cpm(counts); assert np.isfinite(p).all()

def test_too_few_genes_raises():
    ref = _ref(); counts = pd.DataFrame({"s": np.full(100, 5000.0)}, index=[f"X{i}" for i in range(100)])
    with pytest.raises(ValueError): fidelity_scores(mean_log_cpm(counts), ref)

def test_htseq_summary_rows_excluded():
    ref = _ref(); counts = pd.DataFrame({"s": np.full(800, 500.0)}, index=ref.index)
    extra = pd.DataFrame({"s": [1e9]}, index=["__no_feature"])
    a = mean_log_cpm(counts); b = mean_log_cpm(pd.concat([counts, extra]))
    assert np.allclose(a.values, b.loc[a.index].values)

def test_parenchymal_genes_and_restriction(tmp_path):
    from vorganoid.fidelity import parenchymal_genes
    ref = _ref(); ids = list(ref.index)
    rows = ["Gene\tEnsembl\tspec\ttypes"] + [f"G{i}\t{ids[i]}\tenh\t{'Hepatocytes: 5.0' if i % 2 else 'T-cells: 3.0;Hepatocytes: 1.0'}" for i in range(800)]
    f = tmp_path / "hpa.tsv"; f.write_text("\n".join(rows))
    par = parenchymal_genes(str(f)); assert len(par) == 400 and ids[1] in par and ids[0] not in par
    rng = np.random.default_rng(2); counts = pd.DataFrame({"s": rng.poisson(ref["B"] * 200)}, index=ref.index)
    sc = fidelity_scores(mean_log_cpm(counts), ref, genes=par, min_genes=300); assert sc.index[0] == "B"

def test_cli_precomputed_profile_requires_flag(tmp_path, capsys):
    from vorganoid.cli import main
    ref = _ref(); rng = np.random.default_rng(9)
    f = tmp_path / 'ref.gct'; cols = list(ref.columns)
    # GTEx loader consumes a GCT-like table with two header rows.
    with open(f,'w') as fh:
        fh.write('#1.2\n800\t4\nName\tDescription\t'+'\t'.join(cols)+'\n')
        for i,(idx,row) in enumerate(ref.iterrows()):
            fh.write(idx+'\tG'+str(i)+'\t'+'\t'.join(map(str,row.values))+'\n')
    prof = tmp_path / 'profile.csv'
    pd.Series(np.log1p(ref['C'])+rng.normal(0,.001,len(ref)),index=ref.index,name='lcpm').to_csv(prof)
    assert main(['fidelity',str(prof),'--reference',str(f),'--organ-tissue','C','--profile'])==0
    assert '"organ_tissue_rank": 1' in capsys.readouterr().out
    with pytest.raises(SystemExit, match='use --profile'):
        main(['fidelity',str(prof),'--reference',str(f)])
