"""Pre-registered (results/preregistration_clinpgx.md): are core pharmacogenes (ClinPGx/PharmGKB VIP or CPIC-guideline genes;
https://api.clinpgx.org/v1/download/file/data/genes.zip, release 2026-09-05) deficient in liver organoids?
H2 permutation: within each pharmacogene, its per-organ standardized deficits are shuffled across the 5 organs."""
import json, sys, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
sys.path.insert(0, "src"); from geo_fidelity import load_gtex
from deficit_bodymap import deficits, BM

def core_pharmacogenes(df, use_vip=False):
    """Deviation from pre-registration (documented): in release 2026-09-05 the 'Is VIP' column is 'Yes' for all 25,041 genes,
    so the VIP flag is uninformative. Default set = CPIC-dosing-guideline genes only."""
    m = df["Has CPIC Dosing Guideline"].astype(str).str.lower() == "yes"
    if use_vip: m |= df["Is VIP"].astype(str).str.lower() == "yes"
    return set(df.loc[m, "Ensembl Id"].dropna().astype(str).str.split(",").explode().str.strip())

def h2_stat(Z, organs, target="Liver"):
    med = np.median(Z, axis=0); i = organs.index(target)
    return med[i] - np.max(np.delete(med, i))

if __name__ == "__main__":
    pg = core_pharmacogenes(pd.read_csv("data/clinpgx/genes.tsv", sep="\t", dtype=str))
    g = load_gtex(); G = np.log1p(g.drop_duplicates("ens").set_index("ens")[list(BM)]); D = deficits(G)
    sym = g.drop_duplicates("ens").set_index("ens").Description
    L = D["Liver"]; inset = L.index.isin(pg)
    u = mannwhitneyu(L[inset], L[~inset], alternative="greater")
    top = L[inset].sort_values(ascending=False).head(15); top.index = sym.reindex(top.index).values
    organs = list(D); Zs = {o: (D[o] - D[o].mean()) / D[o].std() for o in organs}
    common = sorted(set.intersection(*[set(Zs[o].index) for o in organs]) & pg)
    Z = np.column_stack([Zs[o][common].values for o in organs]); T = h2_stat(Z, organs)
    rng = np.random.default_rng(0); null = np.array([h2_stat(rng.permuted(Z, axis=1), organs) for _ in range(10000)])
    p2 = (1 + (null >= T).sum()) / 10001
    out = {"n_core_pharmacogenes": len(pg), "liver": {"n_pg_tested": int(inset.sum()), "n_other": int((~inset).sum()), "median_deficit_pg": float(L[inset].median()),
           "median_deficit_other": float(L[~inset].median()), "mwu_p_greater": float(u.pvalue), "H1": "PASS" if u.pvalue < 0.01 else "FAIL", "top15_pg_deficits": top.round(2).to_dict()},
           "H2": {"n_pg_all_organs": len(common), "median_std_deficit_by_organ": dict(zip(organs, np.median(Z, axis=0).round(3).tolist())), "stat_liver_minus_max_other": float(T), "perm_p": float(p2), "verdict": "PASS" if p2 < 0.05 else "FAIL"}}
    out["deviation"] = "Is VIP = Yes for all 25,041 genes in the release; set restricted to CPIC-guideline genes (34 in file)"
    json.dump(out, open("results/deficit_clinpgx.json", "w"), indent=1); print(json.dumps(out, indent=1))
