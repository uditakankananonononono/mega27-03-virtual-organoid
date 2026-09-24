"""Organoid-to-tissue fidelity scoring against a tissue reference (e.g. GTEx v8 median TPM).
fidelity_scores: all-gene Spearman (SALL) of mean log1p(CPM) profile vs each reference tissue (log1p TPM)."""
from __future__ import annotations
import numpy as np, pandas as pd
from scipy.stats import rankdata

def mean_log_cpm(counts: pd.DataFrame, min_lib: float = 1e5) -> pd.Series:
    X = counts[~counts.index.astype(str).str.startswith("__")]  # drop HTSeq summary rows
    X = X[~X.index.duplicated()]
    X = X.select_dtypes(include=[np.number]).clip(lower=0)
    X = X.loc[:, X.sum() > min_lib]
    if X.shape[1] == 0: raise ValueError("no library-scale sample columns")
    return np.log1p(X / X.sum() * 1e6).mean(axis=1)

def fidelity_scores(profile: pd.Series, reference: pd.DataFrame, min_genes: int = 500) -> pd.Series:
    profile = profile[~profile.index.duplicated()]; reference = reference[~reference.index.duplicated()]
    common = profile.index.intersection(reference.index)
    if len(common) < min_genes: raise ValueError(f"only {len(common)} shared genes")
    r = rankdata(profile[common]); G = np.log1p(reference.loc[common])
    out = {t: np.corrcoef(r, rankdata(G[t]))[0, 1] for t in G.columns}
    return pd.Series(out).sort_values(ascending=False)

def load_gtex(path: str) -> pd.DataFrame:
    g = pd.read_csv(path, sep="\t", skiprows=2)
    g.index = g.Name.str.split(".").str[0]
    g = g[~g.index.duplicated()]
    return g.drop(columns=["Name", "Description"])
