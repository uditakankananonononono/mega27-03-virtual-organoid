"""Attractor and recovery on the strict subset: series where every sample is annotated as organoid (GEOparse sample metadata),
vs series with no organoid-annotated samples. Output: results/geo_fidelity_strict.json."""
import json, numpy as np, pandas as pd
from scipy.stats import fisher_exact
m = pd.read_csv("results/geo_sample_meta.csv"); R = pd.read_csv("results/geo_fidelity_methods_ranks.csv")
strict = set(m.gse[m.frac_organoid == 1]); none_ = set(m.gse[m.frac_organoid == 0])
rng = np.random.default_rng(3); out = {"n_strict": len(strict), "n_no_organoid_samples": len(none_)}
for meth, g in R.groupby("method"):
    res = {}
    for name, S in [("strict", strict), ("no_organoid", none_)]:
        c = g[g.gse.isin(S)]; fib = (c.best_tissue == "Cells - Cultured fibroblasts")
        bs = [fib.sample(len(fib), replace=True, random_state=int(rng.integers(1e9))).mean() for _ in range(2000)]
        res[name] = {"n": int(len(c)), "top1": float((c.rank_match == 1).mean()), "fibroblast_best_frac": float(fib.mean()),
                     "fibroblast_ci95": [float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))], "fib_count": int(fib.sum())}
    a, b = res["strict"]["fib_count"], res["strict"]["n"] - res["strict"]["fib_count"]
    c_, d = res["no_organoid"]["fib_count"], res["no_organoid"]["n"] - res["no_organoid"]["fib_count"]
    res["fisher_p_strict_vs_none"] = float(fisher_exact([[a, b], [c_, d]])[1]); out[meth] = res
json.dump(out, open("results/geo_fidelity_strict.json", "w"), indent=1)
for k, v in out.items():
    if isinstance(v, dict): print(k, v["strict"]["n"], round(v["strict"]["top1"], 2), round(v["strict"]["fibroblast_best_frac"], 2), [round(x, 2) for x in v["strict"]["fibroblast_ci95"]], "| none:", round(v["no_organoid"]["fibroblast_best_frac"], 2), "p", round(v["fisher_p_strict_vs_none"], 3))
