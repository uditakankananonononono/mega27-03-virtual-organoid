"""Share of each organ's top-100 deficit genes that are HPA parenchymal-only vs non-parenchymal (immune/blood/vascular/stromal) enhanced.
Uses the same HPA split as geo_fidelity_methods.py. Fisher test: parenchymal share in unfaithful (liver, kidney, lung) vs faithful (colon, brain) organs.
Output: results/strict_organ_parenchymal.json"""
import json, pandas as pd
from scipy.stats import fisher_exact
src = open("src/geo_fidelity_methods.py").read().splitlines()
ns = {"pd": pd}; exec("\n".join(src[27:37]), ns)  # lines 28-37: NONPAR set and HPA split
D = pd.read_csv("results/strict_organ_deficits.csv"); out = {}
for o, g in D.groupby("organ"):
    par = g.ens.isin(ns["par_ens"]).sum(); non = g.ens.isin(ns["nonpar_ens"]).sum(); out[o] = {"parenchymal_only": int(par), "nonparenchymal": int(non), "unannotated": int(len(g) - par - non), "par_share_of_annotated": par / max(par + non, 1)}
bad = ["Liver", "Kidney - Cortex", "Lung"]; good = ["Colon - Transverse", "Brain - Cortex"]
t = [[sum(out[o]["parenchymal_only"] for o in bad), sum(out[o]["nonparenchymal"] for o in bad)], [sum(out[o]["parenchymal_only"] for o in good), sum(out[o]["nonparenchymal"] for o in good)]]
out["fisher_unfaithful_vs_faithful"] = {"table": t, "odds_ratio": float(fisher_exact(t)[0]), "p": float(fisher_exact(t)[1])}
json.dump(out, open("results/strict_organ_parenchymal.json", "w"), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
