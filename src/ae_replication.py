"""Pre-registered out-of-sample test (results/preregistration_organ_split.md) on EBI ArrayExpress/BioStudies organoid datasets.
For each curated accession (data/ae/ae_curated.csv): download processed expression files (<40 MB; skip mtx/ATAC), parse numeric tables,
map genes (Ensembl or symbol, via the first column with >=50% GTEx match), average samples into one profile (log1p CPM for counts, as-is for
log-scale tables), score with SALL_PARONLY (identical code path to geo_fidelity_methods.rhos) and take the best rank over the organ's GTEx columns.
Outputs: data/ae/profiles/<acc>.csv.gz, results/ae_replication.csv, results/ae_replication.json"""
import io, os, re, gzip, json, urllib.request, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
src = open("src/geo_fidelity_methods.py").read(); exec(src[:src.index("res = {")])
os.makedirs("data/ae/profiles", exist_ok=True)
ORG = {"liver": "Liver", "kidney": "Kidney - Cortex", "lung": "Lung", "intestine": "Colon - Transverse", "brain": "Brain - Cortex"}
sym2ens = g.drop_duplicates("Description").set_index("Description").ens
def fetch(u):
    loc = "data/ae/raw/" + u.split("/")[-1]
    if os.path.exists(loc): b = open(loc, "rb").read(); return gzip.decompress(b) if u.endswith(".gz") else b  # local cache (curl-downloaded)
    req = urllib.request.Request(u, method="HEAD", headers={"User-Agent": "mega27/1.0"})
    n = int(urllib.request.urlopen(req, timeout=30).headers.get("Content-Length", 0))
    if n > 40_000_000: raise ValueError("too big")
    b = urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "mega27/1.0"}), timeout=120).read()
    return gzip.decompress(b) if u.endswith(".gz") else b
def table(b, name):
    first = b[:5000].decode("utf8", "ignore").lstrip("\ufeff").splitlines()
    first = [l for l in first if not l.startswith("#")][:2]
    sep = "," if re.search(r"\.csv", name) else ("\t" if "\t" in first[0] else r"\s+")
    df = pd.read_csv(io.BytesIO(b), sep=sep, comment="#", **({"engine": "python"} if sep == r"\s+" else {"low_memory": False}))
    ens = set(G.index); best, key = None, 0
    df = df.loc[:, ~df.columns.astype(str).str.contains("Detection|Pval", case=False)]
    if not isinstance(df.index, pd.RangeIndex): df = df.reset_index()  # header shorter than rows: gene ids became the index
    for c in [c for c in df.columns if df[c].dtype == object][:8]:
        v = df[c].astype(str).str.split(".").str[0]
        fe = v.isin(ens).mean(); fs = v.isin(sym2ens.index).mean()
        if max(fe, fs) > key: key, best = max(fe, fs), (c, "ens" if fe >= fs else "sym")
    if key < 0.5: return None
    c, kind = best; ids = df[c].astype(str).str.split(".").str[0]
    if kind == "sym": ids = ids.map(sym2ens)
    X = df.select_dtypes(include=[np.number]); X = X.loc[:, ~X.columns.str.contains("length|start|end|chr|gc", case=False)]
    X.index = ids; X = X[X.index.notna()]; X = X.groupby(level=0).sum() if kind == "ens" else X.groupby(level=0).mean()
    return X
cur = pd.read_csv("data/ae/ae_curated.csv"); rows = []
for r in cur.itertuples():
    prof = []; fc = 0
    if os.path.exists(f"data/ae/profiles/{r.accession}.csv.gz"):
        prof = [pd.read_csv(f"data/ae/profiles/{r.accession}.csv.gz", index_col=0).lcpm]; print(r.accession, "cached profile", flush=True)
    for u in ([] if prof else r.files.split(";")):
        nm = u.split("/")[-1]
        if re.search(r"\.mtx|\.h5$|atac|fragments|meta|index|\.tbi|barcode|feature_counts_GEX|assay_table", nm, re.I): continue
        if "feature_counts" in u:
            fc += 1
            if fc > 3: continue  # per-sample featureCounts files: first 3 samples only (download budget)
        try:
            X = table(fetch(u), nm)
        except Exception as e:
            print(r.accession, nm, "ERR", type(e).__name__, str(e)[:60], flush=True); continue
        if X is None or X.shape[1] == 0 or len(X) < 5000: print(r.accession, nm, "skip (no gene match / too few genes)", flush=True); continue
        X = X.clip(lower=0); logscale = X.values.max() < 50
        prof.append(X.mean(1) if logscale else np.log1p(X / X.sum() * 1e6).mean(1))
        print(r.accession, nm, X.shape, "log" if logscale else "counts", flush=True)
    if not prof: rows.append({"accession": r.accession, "organ": r.organ, "status": "no usable file"}); continue
    p = pd.concat(prof, axis=1).mean(1); p.to_frame("lcpm").to_csv(f"data/ae/profiles/{r.accession}.csv.gz")
    R_ = rhos(p, "SALL_PARONLY"); order = np.argsort(-R_); pos = np.empty(len(T), int); pos[order] = np.arange(1, len(T) + 1)
    o = ORG[r.organ]; rk = int(min(pos[T.index(t)] for t in match[o]))
    rows.append({"accession": r.accession, "organ": r.organ, "gtex_organ": o, "n_genes": int(len(p.index.intersection(G.index))), "rank_match": rk, "best_tissue": T[order[0]], "status": "ok"})
D = pd.DataFrame(rows); D.to_csv("results/ae_replication.csv", index=False)
ok = D[D.status == "ok"]; unf = ok[ok.organ.isin(["liver", "kidney", "lung"])].rank_match; fa = ok[ok.organ.isin(["intestine", "brain"])].rank_match
res = {"n_scored": len(ok), "n_unfaithful": len(unf), "n_faithful": len(fa), "median_rank_unfaithful": float(unf.median()) if len(unf) else None,
       "median_rank_faithful": float(fa.median()) if len(fa) else None, "top1_unfaithful": float((unf == 1).mean()) if len(unf) else None, "top1_faithful": float((fa == 1).mean()) if len(fa) else None}
if len(unf) and len(fa): res["mwu_one_sided_p"] = float(mannwhitneyu(unf, fa, alternative="greater").pvalue)
json.dump(res, open("results/ae_replication.json", "w"), indent=1); print(D.to_string()); print(json.dumps(res, indent=1))
