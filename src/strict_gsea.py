"""GSEApy prerank GSEA of the organ-consensus deficit ranking (delta_g = z(GTEx organ) - z(organoid), strict series, as in
strict_organ_deficits.py but over all genes) against MSigDB Hallmark v2023.2. Tests whether liver/kidney organoids specifically
lack metabolic hallmarks (xenobiotic, fatty acid, bile acid, oxidative phosphorylation) relative to colon/brain organoids.
Outputs: results/strict_gsea_hallmark.csv, results/strict_gsea_summary.json"""
import json, sys, numpy as np, pandas as pd, gseapy as gp
from scipy.stats import rankdata, mannwhitneyu
sys.path.insert(0, "src"); from geo_fidelity import load_gtex
m = pd.read_csv("results/geo_sample_meta.csv"); R = pd.read_csv("results/geo_fidelity_methods_ranks.csv")
r = R[(R.method == "SALL_PARONLY") & R.gse.isin(set(m.gse[m.frac_organoid == 1]))]
g = load_gtex(); T = [c for c in g.columns if c not in ("Name", "Description", "ens")]
G = np.log1p(g.drop_duplicates("ens").set_index("ens")[T]); sym = g.drop_duplicates("ens").set_index("ens").Description
def z(v): q = rankdata(v); return (q - q.mean()) / q.std()
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open("data/ref/h.all.v2023.2.Hs.symbols.gmt")}
rows = []
for organ, grp in r.groupby("organ"):
    if len(grp) < 4: continue
    acc = []
    for s in grp.gse:
        p = pd.read_csv(f"data/geo/profiles/{s}.csv.gz", index_col=0).lcpm; c = p.index.intersection(G.index)
        acc.append(pd.Series(z(G.loc[c, organ].values) - z(p[c].values), index=c))
    d = pd.concat(acc, axis=1, join="inner").mean(1); d.index = sym.reindex(d.index).values
    d = d[pd.notna(d.index)].groupby(level=0).mean().sort_values(ascending=False)
    res = gp.prerank(rnk=d, gene_sets=gmt, permutation_num=500, min_size=15, max_size=500, seed=1, threads=1, outdir=None, verbose=False).res2d
    res["organ"] = organ; rows.append(res[["organ", "Term", "NES", "NOM p-val", "FDR q-val"]])
D = pd.concat(rows); D.columns = ["organ", "term", "NES", "p", "fdr"]; D.to_csv("results/strict_gsea_hallmark.csv", index=False)
MET = ["HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_PEROXISOME", "HALLMARK_ADIPOGENESIS", "HALLMARK_COAGULATION", "HALLMARK_COMPLEMENT"]
W = D.pivot(index="term", columns="organ", values="NES").astype(float)
unf = [o for o in ["Liver", "Kidney - Cortex", "Lung"] if o in W]; fa = [o for o in ["Colon - Transverse", "Brain - Cortex"] if o in W]
met = W.loc[[t for t in MET[:6] if t in W.index]]
out = {"metabolic_NES": met.round(2).to_dict(), "mean_metabolic_NES": met.mean().round(3).to_dict(),
       "mwu_metabolic_NES_liver_kidney_vs_colon_brain": float(mannwhitneyu(met[["Liver", "Kidney - Cortex"]].values.ravel(), met[fa].values.ravel(), alternative="greater").pvalue),
       "top5_per_organ": {o: D[D.organ == o].sort_values("NES", ascending=False).head(5)[["term", "NES", "fdr"]].values.tolist() for o in W.columns}}
json.dump(out, open("results/strict_gsea_summary.json", "w"), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
