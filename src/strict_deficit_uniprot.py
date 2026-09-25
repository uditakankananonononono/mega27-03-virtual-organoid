"""UniProt (rest.uniprot.org) annotation of strict-subset deficit genes: share of reviewed human proteins with a signal peptide
or 'Secreted' subcellular location, per organ, vs a background of 2,000 random GTEx-expressed genes (same query). Fisher test per organ vs background.
Output: results/strict_deficit_uniprot.csv, results/strict_deficit_uniprot.json"""
import io, json, time, urllib.parse, urllib.request, numpy as np, pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
D = pd.read_csv("results/strict_organ_deficits.csv")
def uni(genes):
    rows = []
    for i in range(0, len(genes), 100):
        q = " OR ".join(f"gene_exact:{g}" for g in genes[i:i + 100])
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({"query": f"({q}) AND organism_id:9606 AND reviewed:true", "fields": "gene_primary,ft_signal,cc_subcellular_location", "format": "tsv", "size": 500})
        rows.append(pd.read_csv(io.StringIO(urllib.request.urlopen(url, timeout=90).read().decode()), sep="\t")); time.sleep(0.3)
    U = pd.concat(rows); U.columns = ["gene", "signal", "loc"]; U = U.drop_duplicates("gene")
    U["secreted"] = U.signal.notna() | U["loc"].fillna("").str.contains("Secreted"); return U
sys_ = __import__("sys"); sys_.path.insert(0, "src"); from geo_fidelity import load_gtex
g = load_gtex(); expr = g[g.drop(columns=["Name", "Description", "ens"]).max(1) > 5].Description.dropna().unique()
bg = list(np.random.default_rng(0).choice(expr, 2000, replace=False))
B = uni(bg); out = {"background": {"n": len(B), "secreted": int(B.secreted.sum())}}; rows = []
for o, grp in D.groupby("organ"):
    U = uni(list(grp.symbol.dropna())); a, b = int(U.secreted.sum()), len(U) - int(U.secreted.sum())
    OR, p = fisher_exact([[a, b], [out["background"]["secreted"], len(B) - out["background"]["secreted"]]])
    rows.append({"organ": o, "n_mapped": len(U), "secreted": a, "share": a / len(U), "bg_share": B.secreted.mean(), "OR": OR, "p": p})
R = pd.DataFrame(rows); R["q_bh"] = multipletests(R.p, method="fdr_bh")[1]; R.to_csv("results/strict_deficit_uniprot.csv", index=False)
out["organs"] = R.to_dict("records"); json.dump(out, open("results/strict_deficit_uniprot.json", "w"), indent=1, default=float); print(R.round(4).to_string())
