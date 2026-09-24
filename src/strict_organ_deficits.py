"""Which organs do strict organoid series recover, and what do the unfaithful ones lack?
Strict series: every GSM annotated organoid (results/geo_sample_meta.csv). Recovery = GTEx rank of the labelled organ (SALL_PARONLY).
Deficit score per gene and series: delta_g = z(G_g,organ) - z(p_g), rank-based z within profile; consensus per organ = mean delta.
Top-100 deficit genes per organ are annotated with the STRING functional-enrichment API (string-db.org) and Enrichr PanglaoDB_Augmented_2021.
Outputs: results/strict_organ_recovery.csv, results/strict_organ_deficits.csv, results/strict_organ_deficit_enrichment.csv, results/strict_organ_summary.json"""
import json, sys, time, urllib.request, urllib.parse, numpy as np, pandas as pd
from scipy.stats import rankdata, fisher_exact
sys.path.insert(0, "src"); from geo_fidelity import load_gtex
m = pd.read_csv("results/geo_sample_meta.csv"); R = pd.read_csv("results/geo_fidelity_methods_ranks.csv")
S = set(m.gse[m.frac_organoid == 1]); r = R[(R.method == "SALL_PARONLY") & R.gse.isin(S)].copy()
rec = r.groupby("organ").agg(n=("gse", "size"), top1=("rank_match", lambda x: (x == 1).mean()), median_rank=("rank_match", "median")).reset_index()
rec.to_csv("results/strict_organ_recovery.csv", index=False)
good = r[r.organ.isin(["Colon - Transverse", "Brain - Cortex"])]; bad = r[r.organ.isin(["Liver", "Kidney - Cortex", "Lung"])]
tab = [[int((good.rank_match == 1).sum()), int((good.rank_match != 1).sum())], [int((bad.rank_match == 1).sum()), int((bad.rank_match != 1).sum())]]
out = {"n_strict": len(r), "fisher_colon_brain_vs_liver_kidney_lung": {"table": tab, "p": float(fisher_exact(tab)[1])}}
g = load_gtex(); T = [c for c in g.columns if c not in ("Name", "Description", "ens")]
G = np.log1p(g.drop_duplicates("ens").set_index("ens")[T]); sym = g.drop_duplicates("ens").set_index("ens").Description
def z(v): q = rankdata(v); return (q - q.mean()) / q.std()
rows = []
for organ, grp in r.groupby("organ"):
    if len(grp) < 4: continue
    acc = None
    for s in grp.gse:
        p = pd.read_csv(f"data/geo/profiles/{s}.csv.gz", index_col=0).lcpm; c = p.index.intersection(G.index)
        d = pd.Series(z(G.loc[c, organ].values) - z(p[c].values), index=c); acc = d.to_frame(s) if acc is None else acc.join(d.to_frame(s), how="inner")
    cons = acc.mean(1).sort_values(ascending=False); frac = (acc > 1).mean(1)
    for e in cons.index[:100]: rows.append({"organ": organ, "ens": e, "symbol": sym.get(e), "delta": cons[e], "frac_series_delta_gt1": frac[e], "n_series": acc.shape[1]})
D = pd.DataFrame(rows); D.to_csv("results/strict_organ_deficits.csv", index=False)
def string_enrich(genes):
    data = urllib.parse.urlencode({"identifiers": "%0d".join(genes), "species": 9606, "caller_identity": "mega27"}).encode()
    q = urllib.request.Request("https://string-db.org/api/json/enrichment", data=data, headers={"User-Agent": "mega27/1.0"})
    return json.load(urllib.request.urlopen(q, timeout=90))
def enrichr(genes, lib="PanglaoDB_Augmented_2021"):
    b = "----b"; body = (f"--{b}\r\nContent-Disposition: form-data; name=\"list\"\r\n\r\n" + "\n".join(genes) + f"\r\n--{b}--\r\n").encode()
    uid = json.load(urllib.request.urlopen(urllib.request.Request("https://maayanlab.cloud/Enrichr/addList", data=body, headers={"Content-Type": f"multipart/form-data; boundary={b}"}), timeout=60))["userListId"]
    return json.load(urllib.request.urlopen(f"https://maayanlab.cloud/Enrichr/enrich?userListId={uid}&backgroundType={lib}", timeout=60))[lib]
E = []
for organ, grp in D.groupby("organ"):
    genes = [x for x in grp.symbol.dropna()]
    for x in string_enrich(genes):
        if x["category"] in ("Process", "KEGG", "Component") and x["fdr"] < 0.05: E.append({"organ": organ, "source": "STRING:" + x["category"], "term": x["description"], "fdr": x["fdr"], "n_genes": x["number_of_genes"], "genes": x["inputGenes"] if isinstance(x["inputGenes"], str) else ",".join(x["inputGenes"])})
    for x in enrichr(genes)[:8]: E.append({"organ": organ, "source": "Enrichr:PanglaoDB", "term": x[1], "fdr": x[6], "n_genes": len(x[5]), "genes": ",".join(x[5])})
    time.sleep(1)
EE = pd.DataFrame(E); EE.to_csv("results/strict_organ_deficit_enrichment.csv", index=False)
out["top_deficit_genes"] = {o: list(grp.symbol[:15]) for o, grp in D.groupby("organ")}
json.dump(out, open("results/strict_organ_summary.json", "w"), indent=1, default=str)
print(rec.to_string()); print(out["fisher_colon_brain_vs_liver_kidney_lung"])
for o, grp in EE.groupby("organ"):
    print("==", o, out["top_deficit_genes"][o][:12])
    for src in grp.source.unique(): print("  ", src, grp[grp.source == src].sort_values("fdr").head(3)[["term", "fdr"]].values.tolist())
