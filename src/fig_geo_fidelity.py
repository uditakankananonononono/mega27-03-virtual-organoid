"""Figure: per-organ top-1 tissue recovery and fibroblast best-match share (clean subset, SALL and CENT_NOCULT)."""
import json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, numpy as np, pandas as pd
a = pd.read_csv("results/geo_label_audit.csv"); R = pd.read_csv("results/geo_fidelity_methods_ranks.csv")
clean = set(a.gse[a.label_consistent & a.label_unambiguous & a.organoid_in_summary])
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for k, m in enumerate(["SALL", "CENT_NOCULT"]):
    c = R[(R.method == m) & R.gse.isin(clean)].copy(); c["short"] = c.organ.str.split(" - ").str[0]
    t = c.groupby("short").agg(n=("gse", "size"), top1=("rank_match", lambda r: (r == 1).mean()),
                               fib=("best_tissue", lambda b: (b == "Cells - Cultured fibroblasts").mean()))
    t = t[t.n >= 3].sort_values("n", ascending=False); x = np.arange(len(t))
    ax[k].bar(x - 0.2, t.top1, 0.4, label="organ of origin is top-1"); ax[k].bar(x + 0.2, t.fib, 0.4, label="cultured fibroblasts top-1")
    ax[k].set_xticks(x); ax[k].set_xticklabels([f"{i}\n(n={n})" for i, n in zip(t.index, t.n)], fontsize=7, rotation=45)
    ax[k].set_title({"SALL": "All-gene Spearman", "CENT_NOCULT": "Centred Pearson, culture genes removed"}[m], fontsize=9); ax[k].set_ylim(0, 1)
ax[0].set_ylabel("fraction of series"); ax[0].legend(fontsize=7); plt.tight_layout(); plt.savefig("results/figures/fig_geo_fidelity.png", dpi=200)
