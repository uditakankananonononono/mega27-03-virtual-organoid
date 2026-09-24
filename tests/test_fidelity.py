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
