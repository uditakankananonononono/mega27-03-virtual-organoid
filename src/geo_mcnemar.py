"""Paired McNemar tests (statsmodels) of top-1 organ recovery: base method vs gene-set-restricted variant, per series."""
import json, pandas as pd
from statsmodels.stats.contingency_tables import mcnemar
R = pd.read_csv("results/geo_fidelity_methods_ranks.csv"); W = R.pivot(index="gse", columns="method", values="rank_match") == 1
out = {}
for base in ["SALL", "CENT"]:
    for v in ["NOCULT", "NONPAR", "PARONLY", "PC"]:
        a, b = W[base], W[f"{base}_{v}"]
        t = [[int((a & b).sum()), int((a & ~b).sum())], [int((~a & b).sum()), int((~a & ~b).sum())]]
        r = mcnemar(t, exact=True)
        out[f"{base}_vs_{base}_{v}"] = {"table": t, "top1_base": float(a.mean()), "top1_variant": float(b.mean()), "p_exact": float(r.pvalue)}
json.dump(out, open("results/geo_mcnemar.json", "w"), indent=1)
for k, v in out.items(): print(k, round(v["top1_base"], 3), round(v["top1_variant"], 3), v["table"][0][1], v["table"][1][0], f"{v['p_exact']:.2g}")
