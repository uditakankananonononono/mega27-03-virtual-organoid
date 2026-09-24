"""Re-test tissue recovery and the fibroblast attractor on audited series only
(title organ confirmed by summary, single organ in summary, 'organoid' in summary). Output: results/geo_fidelity_clean.json."""
import json, numpy as np, pandas as pd
a = pd.read_csv("results/geo_label_audit.csv"); R = pd.read_csv("results/geo_fidelity_methods_ranks.csv")
clean = set(a.gse[a.label_consistent & a.label_unambiguous & a.organoid_in_summary])
rng = np.random.default_rng(1); out = {"n_clean": len(clean)}
for m, g in R.groupby("method"):
    c = g[g.gse.isin(clean)]; fib = (c.best_tissue == "Cells - Cultured fibroblasts")
    # binomial-style bootstrap CI over series
    bs = [c.sample(len(c), replace=True, random_state=int(rng.integers(1e9))) for _ in range(2000)]
    out[m] = {"n": int(len(c)), "top1": float((c.rank_match == 1).mean()), "top5": float((c.rank_match <= 5).mean()),
              "fibroblast_best_frac": float(fib.mean()),
              "fibroblast_best_ci95": [float(np.quantile([(b.best_tissue == "Cells - Cultured fibroblasts").mean() for b in bs], q)) for q in (0.025, 0.975)],
              "fibroblast_best_by_organ": c.assign(f=fib).groupby("organ").f.agg(["count", "sum"]).to_dict("index")}
json.dump(out, open("results/geo_fidelity_clean.json", "w"), indent=1)
print(json.dumps({m: {k: v for k, v in out[m].items() if k != "fibroblast_best_by_organ"} for m in out if m != "n_clean"}), out["n_clean"])
