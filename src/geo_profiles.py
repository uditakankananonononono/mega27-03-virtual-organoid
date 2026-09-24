"""Cache per-GSE mean log1p(CPM) profiles (library-scale columns only) for offline method checks.
Re-downloads the exact file recorded in results/geo_fidelity.csv. Output: data/geo/profiles/<GSE>.csv.gz (gitignored)."""
import os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from geo_fidelity import load_gtex, list_suppl, fetch, parse, match_genes
os.makedirs("data/geo/profiles", exist_ok=True)
g = load_gtex(); d = pd.read_csv("results/geo_fidelity.csv")
for _, r in d.iterrows():
    out = f"data/geo/profiles/{r.gse}.csv.gz"
    if os.path.exists(out): continue
    try:
        url, _ = list_suppl(r.gse); df = parse(fetch(url + r.file), r.file); df, key = match_genes(df, g)
        X = df.clip(lower=0); X = X.loc[:, X.sum() > 1e5]
        prof = np.log1p(X / X.sum() * 1e6).mean(axis=1)
        if key == "Description":  # map symbols to Ensembl for a common index
            m = g.drop_duplicates("Description").set_index("Description").ens
            prof = prof[prof.index.isin(m.index)]; prof.index = m.loc[prof.index].values
        prof = prof[~prof.index.duplicated()].astype("float32")
        prof.to_csv(out, header=["lcpm"])
        print(r.gse, len(prof), flush=True)
    except Exception as e: print(r.gse, "fail", type(e).__name__, flush=True)
    time.sleep(0.2)
