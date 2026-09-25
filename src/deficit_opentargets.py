"""Pre-registered test (results/preregistration_opentargets.md): organoid deficit genes vs Open Targets organ-disease targets."""
import json, requests, numpy as np, pandas as pd
from scipy.stats import fisher_exact
URL = "https://api.platform.opentargets.org/api/v4/graphql"
DIS = {"Liver": "MONDO_0005154", "Kidney - Cortex": "MONDO_0005240", "Lung": "MONDO_0005275",
       "Brain - Cortex": "MONDO_0005560", "Colon - Transverse": "MONDO_0003409"}
Q = "query($id:String!,$i:Int!){disease(efoId:$id){name associatedTargets(page:{index:$i,size:500}){count rows{score target{id approvedSymbol}}}}}"

def top_targets(efo, n=500):
    r = requests.post(URL, json={"query": Q, "variables": {"id": efo, "i": 0}}, timeout=60).json()["data"]["disease"]
    rows = sorted(r["associatedTargets"]["rows"], key=lambda x: -x["score"])[:n]
    return r["name"], r["associatedTargets"]["count"], {x["target"]["id"] for x in rows}

def diag_dominance(M):
    k = M.shape[0]; off = (M.sum() - np.trace(M)) / (k * k - k)
    return np.trace(M) / k - off

def overlap_matrix(lists, sets, organs):
    return np.array([[len(lists[a] & sets[b]) for b in organs] for a in organs], float)

if __name__ == "__main__":
    d = pd.read_csv("results/strict_organ_deficits.csv"); organs = list(DIS)
    lists = {o: set(d[d.organ == o].ens) for o in organs}; sets = {}; meta = {}
    for o, efo in DIS.items():
        name, cnt, s = top_targets(efo); sets[o] = s; meta[o] = {"efo": efo, "name": name, "n_assoc_total": cnt, "n_top": len(s)}
    M = overlap_matrix(lists, sets, organs); D = diag_dominance(M)
    genes = d[d.organ.isin(organs)][["organ", "ens"]].reset_index(drop=True); rng = np.random.default_rng(0); null = []
    for _ in range(10000):
        lab = rng.permutation(genes.organ.values)
        L = {o: set(genes.ens[lab == o]) for o in organs}; null.append(diag_dominance(overlap_matrix(L, sets, organs)))
    null = np.array(null); p = (1 + (null >= D).sum()) / (1 + len(null))
    per = {}
    for o in organs:
        own = lists[o]; oth = set().union(*[lists[x] for x in organs if x != o]) - own
        a = len(own & sets[o]); b = len(own) - a; c = len(oth & sets[o]); e = len(oth) - c
        OR, pf = fisher_exact([[a, b], [c, e]], alternative="greater"); per[o] = {"own_hits": a, "own_n": len(own), "other_hits": c, "other_n": len(oth), "OR": float(OR), "p": float(pf)}
    out = {"diseases": meta, "overlap_matrix": {"rows_deficit_lists": organs, "cols_disease_sets": organs, "M": M.tolist()},
           "D": float(D), "null_mean": float(null.mean()), "perm_p": float(p), "per_organ": per, "verdict": "PASS" if p < 0.05 else "FAIL"}
    json.dump(out, open("results/deficit_opentargets.json", "w"), indent=1); print(json.dumps(out, indent=1))
