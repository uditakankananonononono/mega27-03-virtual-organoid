"""Method checks and nulls for the GEO organoid-fidelity scan (offline, from cached profiles).
Variants: S3000 = Spearman on 3000 most tissue-variable genes (original); SALL = Spearman on all shared genes;
CENT = Pearson after centring both organoid and tissue profiles on the GTEx across-tissue mean (tissue-specific signal only).
Nulls: (1) analytic chance k/54 per series for top-1 (k = GTEx columns counted as a match);
(2) permutation of organ labels across series (10,000 draws), preserving the organ mix."""
import glob, json, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
sys.path.insert(0, "src")
from geo_fidelity import load_gtex, ORGANS
g = load_gtex(); T = [c for c in g.columns if c not in ("Name", "Description", "ens")]
G = np.log1p(g.drop_duplicates("ens").set_index("ens")[T])
d = pd.read_csv("results/geo_fidelity.csv")
match = {o[1][0]: o[1] for o in ORGANS}
P = {}
for f in glob.glob("data/geo/profiles/*.csv.gz"):
    P[f.split("/")[-1][:-7]] = pd.read_csv(f, index_col=0).lcpm
d = d[d.gse.isin(P)].reset_index(drop=True)
var_top = G.var(axis=1).sort_values().index[-3000:]
def rhos(p, method):
    c = p.index.intersection(G.index)
    if method == "S3000":
        c = c.intersection(var_top); return np.array([spearmanr(p[c], G.loc[c, t])[0] for t in T])
    if method == "SALL":
        R = rankdata(p[c]); GR = G.loc[c].rank().values
        return np.array([np.corrcoef(R, GR[:, j])[0, 1] for j in range(len(T))])
    if method == "CENT":
        mu = G.loc[c].mean(axis=1); x = p[c] - p[c].mean() - (mu - mu.mean())
        Y = G.loc[c].sub(mu, axis=0); Y = Y - Y.mean()
        return (Y.values * x.values[:, None]).sum(0) / (np.linalg.norm(x) * np.linalg.norm(Y.values, axis=0))
res = {"n_series": int(len(d)), "variants": {}}
rng = np.random.default_rng(0)
rows = []
for m in ["S3000", "SALL", "CENT"]:
    R = np.array([rhos(P[s], m) for s in d.gse])
    order = np.argsort(-R, axis=1)
    pos = np.argsort(order, axis=1) + 1  # pos[i, t] = rank of tissue t in series i
    organs = sorted(match); oi = {o: k for k, o in enumerate(organs)}
    M = np.stack([pos[:, [T.index(t) for t in match[o]]].min(1) for o in organs], 1)
    lab = np.array([oi[o] for o in d.organ]); idx = np.arange(len(d))
    rk = M[idx, lab]
    chance1 = float(np.mean([len(match[o]) / len(T) for o in d.organ]))
    perm1 = np.array([(M[idx, rng.permutation(lab)] == 1).mean() for _ in range(10000)])
    best = pd.Series([T[o[0]] for o in order]).value_counts()
    res["variants"][m] = {"top1": float((rk == 1).mean()), "top5": float((rk <= 5).mean()), "median_rank": float(np.median(rk)),
        "analytic_chance_top1": chance1, "perm_top1_mean": float(perm1.mean()), "perm_top1_p": float((1 + (perm1 >= (rk == 1).mean()).sum()) / (1 + len(perm1))),
        "top_best_tissues": best.head(5).to_dict(),
        "per_organ_top1": pd.DataFrame({"o": d.organ, "t": rk == 1}).groupby("o").t.agg(["count", "mean"]).round(3).to_dict("index")}
    for s, o, r in zip(d.gse, d.organ, rk): rows.append({"gse": s, "organ": o, "method": m, "rank_match": int(r), "best_tissue": T[order[list(d.gse).index(s)][0]]})
    print(m, res["variants"][m]["top1"], res["variants"][m]["perm_top1_mean"], res["variants"][m]["perm_top1_p"], flush=True)
json.dump(res, open("results/geo_fidelity_methods.json", "w"), indent=1)
pd.DataFrame(rows).to_csv("results/geo_fidelity_methods_ranks.csv", index=False)
