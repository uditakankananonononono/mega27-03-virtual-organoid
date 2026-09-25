"""Pre-registered (results/preregistration_decoupler.md): decoupler ulm TF (CollecTRI) and pathway (PROGENy) activity of organoid deficits.
Nets fetched with decoupler.op.collectri / op.progeny(top=500) from OmniPath, cached in data/ref/. Outputs results/deficit_decoupler*.{csv,json}"""
import json, sys, numpy as np, pandas as pd, decoupler as dc
from scipy.stats import rankdata, binomtest
sys.path.insert(0, "src"); from geo_fidelity import load_gtex
MASTER = {"Liver": ["HNF4A", "HNF1A", "NR1H4", "CEBPA"], "Kidney - Cortex": ["HNF1B", "HNF4A", "PAX2", "PAX8"], "Lung": ["NKX2-1", "FOXA2"],
          "Colon - Transverse": ["CDX2", "CDX1", "HNF4A"], "Brain - Cortex": ["NEUROD2", "NEUROD6", "TBR1"]}

def top_frac_hits(act, organ_tfs, q=0.10):
    """act: Series of TF activities (higher = more deficient). Returns (hits, n_tested)."""
    cut = act.quantile(1 - q); t = [x for x in organ_tfs if x in act.index]
    return sum(act[x] >= cut for x in t), len(t)

def deficit_matrix():
    m = pd.read_csv("results/geo_sample_meta.csv"); R = pd.read_csv("results/geo_fidelity_methods_ranks.csv")
    r = R[(R.method == "SALL_PARONLY") & R.gse.isin(set(m.gse[m.frac_organoid == 1]))]
    g = load_gtex(); T = [c for c in g.columns if c not in ("Name", "Description", "ens")]
    G = np.log1p(g.drop_duplicates("ens").set_index("ens")[T]); sym = g.drop_duplicates("ens").set_index("ens").Description
    z = lambda v: (rankdata(v) - rankdata(v).mean()) / rankdata(v).std(); rows = {}
    for organ, grp in r.groupby("organ"):
        if organ not in MASTER or len(grp) < 4: continue
        acc = []
        for s in grp.gse:
            p = pd.read_csv(f"data/geo/profiles/{s}.csv.gz", index_col=0).lcpm; c = p.index.intersection(G.index)
            acc.append(pd.Series(z(G.loc[c, organ].values) - z(p[c].values), index=c))
        d = pd.concat(acc, axis=1, join="inner").mean(1); d.index = sym.reindex(d.index).values
        rows[organ] = d[pd.notna(d.index)].groupby(level=0).mean()
    return pd.DataFrame(rows).T.dropna(axis=1)

if __name__ == "__main__":
    X = deficit_matrix(); col = pd.read_csv("data/ref/collectri_human.csv"); prg = pd.read_csv("data/ref/progeny_human_top500.csv")
    tf, tfp = dc.mt.ulm(data=X, net=col[["source", "target", "weight"]], tmin=5)
    pw, pwp = dc.mt.ulm(data=X, net=prg[["source", "target", "weight"]], tmin=5)
    tf.to_csv("results/deficit_decoupler_tf.csv"); pw.to_csv("results/deficit_decoupler_progeny.csv")
    per = {}; H = N = 0
    for o in X.index:
        h, n = top_frac_hits(tf.loc[o], MASTER[o]); H += h; N += n
        a = tf.loc[o]; per[o] = {"hits": h, "tested": n, "ranks_of_master": {x: int((a > a[x]).sum()) + 1 for x in MASTER[o] if x in a.index}, "n_tfs": int(a.size),
                                 "top10_deficit_TFs": a.sort_values(ascending=False).head(10).round(2).to_dict(),
                                 "progeny": pw.loc[o].round(2).to_dict()}
    bt = binomtest(H, N, 0.10, alternative="greater")
    out = {"organs": list(X.index), "n_genes": int(X.shape[1]), "pooled_hits": H, "pooled_tested": N, "binom_p": float(bt.pvalue), "verdict": "PASS" if bt.pvalue < 0.05 else "FAIL", "per_organ": per}
    json.dump(out, open("results/deficit_decoupler.json", "w"), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
