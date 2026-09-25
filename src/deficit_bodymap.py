"""Pre-registered (results/preregistration_bodymap.md) reference-robustness test: GTEx vs Illumina Body Map (Expression Atlas
E-MTAB-513 baseline TPMs, https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-513/E-MTAB-513-tpms.tsv;
cells are 5 comma-separated quantiles, the middle one = median is used)."""
import json, sys, numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr, hypergeom
sys.path.insert(0, "src"); from geo_fidelity import load_gtex
import deficit_opentargets as ot
BM = {"Liver": "g6", "Kidney - Cortex": "g5", "Lung": "g11", "Brain - Cortex": "g9", "Colon - Transverse": "g13"}
z = lambda v: (rankdata(v) - rankdata(v).mean()) / rankdata(v).std()

def parse_cell(x):
    return float(str(x).split(",")[2]) if "," in str(x) else float(x)

def load_bodymap(path="data/atlas/E-MTAB-513-tpms.tsv"):
    b = pd.read_csv(path, sep="\t", comment="#").set_index("GeneID")
    return np.log1p(b[list(BM.values())].map(parse_cell)).rename(columns={v: k for k, v in BM.items()})

def deficits(ref):
    m = pd.read_csv("results/geo_sample_meta.csv"); R = pd.read_csv("results/geo_fidelity_methods_ranks.csv")
    r = R[(R.method == "SALL_PARONLY") & R.gse.isin(set(m.gse[m.frac_organoid == 1]))]; out = {}
    for organ, grp in r.groupby("organ"):
        if organ not in BM or len(grp) < 4: continue
        acc = []
        for s in grp.gse:
            p = pd.read_csv(f"data/geo/profiles/{s}.csv.gz", index_col=0).lcpm; c = p.index.intersection(ref.index)
            acc.append(pd.Series(z(ref.loc[c, organ].values) - z(p[c].values), index=c))
        out[organ] = pd.concat(acc, axis=1, join="inner").mean(1)
    return out

if __name__ == "__main__":
    g = load_gtex(); G = np.log1p(g.drop_duplicates("ens").set_index("ens")[list(BM)]); B = load_bodymap()
    common = G.index.intersection(B.index); dG = deficits(G.loc[common]); dB = deficits(B.loc[common]); res = {}; lists = {}
    for o in dG:
        c = dG[o].index.intersection(dB[o].index); a, b = dG[o][c], dB[o][c]; rho = spearmanr(a, b).correlation
        tA, tB = set(a.nlargest(100).index), set(b.nlargest(100).index); k = len(tA & tB)
        res[o] = {"n_genes": len(c), "spearman": float(rho), "top100_overlap": k, "hypergeom_p": float(hypergeom.sf(k - 1, len(c), 100, 100))}; lists[o] = tB
    organs = list(ot.DIS); sets = {o: ot.top_targets(ot.DIS[o])[2] for o in organs}
    M = ot.overlap_matrix(lists, sets, organs); D = ot.diag_dominance(M)
    genes = pd.DataFrame([(o, x) for o in organs for x in lists[o]], columns=["organ", "ens"]); rng = np.random.default_rng(0); null = []
    for _ in range(10000):
        lab = rng.permutation(genes.organ.values); null.append(ot.diag_dominance(ot.overlap_matrix({o: set(genes.ens[lab == o]) for o in organs}, sets, organs)))
    p = (1 + (np.array(null) >= D).sum()) / 10001
    primary = all(v["spearman"] >= 0.5 and v["hypergeom_p"] < 1e-3 for v in res.values())
    out = {"n_common_genes": len(common), "per_organ": res, "primary_verdict": "PASS" if primary else "FAIL",
           "opentargets_rerun": {"M": M.tolist(), "D": float(D), "perm_p": float(p), "verdict": "PASS" if p < 0.05 else "FAIL"}}
    json.dump(out, open("results/deficit_bodymap.json", "w"), indent=1); print(json.dumps(out, indent=1))
