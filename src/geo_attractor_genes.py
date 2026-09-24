"""Which genes drive the culture-fibroblast attractor, and what are they (g:Profiler)?
For clean-subset series whose best match (SALL) is GTEx cultured fibroblasts, score each gene by
d_g = mean_s [ z(p_sg) * ( z(G_g,fib) - z(G_g,organ_s) ) ], where z = rank-based standard scores within the profile.
Top 300 genes by d_g are tested with g:Profiler (GO:BP, REAC, KEGG; g:SCS correction; background = genes scored).
Outputs: results/attractor_genes.csv, results/attractor_gprofiler.csv."""
import json, urllib.request, glob
import numpy as np, pandas as pd
from scipy.stats import rankdata
import sys; sys.path.insert(0, "src"); from geo_fidelity import load_gtex, ORGANS
g = load_gtex(); T = [c for c in g.columns if c not in ("Name", "Description", "ens")]
G = np.log1p(g.drop_duplicates("ens").set_index("ens")[T]); sym = g.drop_duplicates("ens").set_index("ens").Description
a = pd.read_csv("results/geo_label_audit.csv"); R = pd.read_csv("results/geo_fidelity_methods_ranks.csv")
clean = set(a.gse[a.label_consistent & a.label_unambiguous & a.organoid_in_summary])
S = R[(R.method == "SALL") & R.gse.isin(clean) & (R.best_tissue == "Cells - Cultured fibroblasts")]
match = {o[1][0]: o[1][0] for o in ORGANS}
def z(v): r = rankdata(v); return (r - r.mean()) / r.std()
acc = None; n = 0
for s, organ in zip(S.gse, S.organ):
    p = pd.read_csv(f"data/geo/profiles/{s}.csv.gz", index_col=0).lcpm; c = p.index.intersection(G.index)
    d = pd.Series(z(p[c].values) * (z(G.loc[c, "Cells - Cultured fibroblasts"].values) - z(G.loc[c, organ].values)), index=c)
    acc = d if acc is None else acc.add(d, fill_value=0); n += 1
acc = (acc / n).sort_values(ascending=False)
top = acc.index[:300]
pd.DataFrame({"ens": acc.index[:1000], "symbol": sym.reindex(acc.index[:1000]).values, "d": acc.values[:1000]}).to_csv("results/attractor_genes.csv", index=False)
body = {"organism": "hsapiens", "query": list(top), "sources": ["GO:BP", "REAC", "KEGG"], "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS", "domain_scope": "custom", "background": list(acc.index), "no_evidences": True}
req = urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", data=json.dumps(body).encode(), headers={"Content-Type": "application/json"})
res = json.load(urllib.request.urlopen(req, timeout=100))["result"]
E = pd.DataFrame([{"source": r["source"], "term_id": r["native"], "name": r["name"], "p_adj": r["p_value"], "term_size": r["term_size"],
                   "intersection": r["intersection_size"]} for r in res]).sort_values("p_adj")
E.to_csv("results/attractor_gprofiler.csv", index=False)
print(n, "series;", len(E), "terms"); print(E.head(15).to_string(index=False)); print(list(sym.reindex(top[:30]).values))
