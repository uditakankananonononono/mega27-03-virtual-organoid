"""Organoid-to-tissue transcriptome fidelity across GEO series.
For each human organoid bulk RNA-seq GSE: download one supplementary gene-level matrix, convert each sample to
log1p(CPM), average samples, and compute Spearman rho against every GTEx v8 tissue median profile.
Records the expected organ (from the series title), the rank of the matched GTEx tissue among 54, and a z-score.
Resumable: appends to results/geo_fidelity.csv. Real downloads only; failures are logged, never imputed."""
import re, io, gzip, sys, os, csv, time, urllib.request
import numpy as np, pandas as pd
from scipy.stats import spearmanr
ORGANS = [  # (regex on title, GTEx tissue columns that count as a match)
 (r"colon|colorectal|intestin|gut|ileum|ileal|duoden|jejun|enteroid|colonoid|rectal", ["Colon - Transverse", "Small Intestine - Terminal Ileum", "Colon - Sigmoid"]),
 (r"liver|hepat|cholangio|biliary", ["Liver"]),
 (r"brain|cerebral|cortical|cortex|neural organoid|forebrain|midbrain|assembloid", ["Brain - Cortex", "Brain - Frontal Cortex (BA9)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Caudate (basal ganglia)", "Brain - Hippocampus", "Brain - Putamen (basal ganglia)", "Brain - Nucleus accumbens (basal ganglia)", "Brain - Hypothalamus", "Brain - Substantia nigra", "Brain - Amygdala", "Brain - Cerebellum", "Brain - Cerebellar Hemisphere", "Brain - Spinal cord (cervical c-1)"]),
 (r"kidney|renal|nephron|tubul", ["Kidney - Cortex", "Kidney - Medulla"]),
 (r"lung|airway|alveol|bronch", ["Lung"]),
 (r"pancrea", ["Pancreas"]),
 (r"stomach|gastric", ["Stomach"]),
 (r"prostat", ["Prostate"]),
 (r"breast|mammary", ["Breast - Mammary Tissue"]),
 (r"cardiac|heart|cardio", ["Heart - Left Ventricle", "Heart - Atrial Appendage"]),
 (r"esophag|oesophag", ["Esophagus - Mucosa"]),
 (r"skin|epiderm|keratinocyte", ["Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"]),
 (r"thyroid", ["Thyroid"]), (r"ovar", ["Ovary"]), (r"endometri|uter", ["Uterus"]), (r"bladder|urothel", ["Bladder"]),
]
def organ_of(title):
    hits = [i for i, (rx, _) in enumerate(ORGANS) if re.search(rx, title, re.I)]
    return hits[0] if len(hits) == 1 else None
def load_gtex(path="data/ref/gtex_v8_median_tpm.gct.gz"):
    g = pd.read_csv(path, sep="\t", skiprows=2)
    g["ens"] = g.Name.str.split(".").str[0]
    return g
def fetch(url, maxb=40_000_000):
    req = urllib.request.urlopen(url, timeout=90)
    b = req.read(maxb + 1)
    if len(b) > maxb: raise ValueError("too big")
    return b
def list_suppl(gse):
    stub = gse[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/"
    html = fetch(url, 2_000_000).decode("utf8", "ignore")
    files = re.findall(r'href="([^"]+\.(?:txt|tsv|csv)(?:\.gz)?)"', html, re.I)
    sizes = re.findall(r'href="[^"]+"[^\n]*?\s(\d+(?:\.\d+)?[KMG]?)\s*$', html, re.M)
    return url, [f for f in files if not re.search(r"filelist|readme|meta|annot|sample|barcode|peak|diff|deg|de_|result", f, re.I)]
def parse(b, name):
    if name.endswith(".gz"):
        with gzip.open(io.BytesIO(b)) as fh: b = fh.read(150_000_001)
        if len(b) > 150_000_000: raise ValueError("decompressed too big")
    sep = "," if re.search(r"\.csv", name, re.I) else "\t"
    df = pd.read_csv(io.BytesIO(b), sep=sep, index_col=0, low_memory=False)
    df = df.select_dtypes(include=[np.number])
    return df
def match_genes(df, g):
    idx = df.index.astype(str).str.split(".").str[0]
    ens = set(g.ens); sym = set(g.Description)
    fe = np.mean([i in ens for i in idx[:2000]]); fs = np.mean([i in sym for i in idx[:2000]])
    key = "ens" if fe >= fs else "Description"
    if max(fe, fs) < 0.5: return None
    df = df.copy(); df.index = idx if key == "ens" else df.index.astype(str)
    df = df[~df.index.duplicated()]
    return df, key
def score(df, key, g, tissues):
    X = df.clip(lower=0)
    X = X.loc[:, X.sum() > 1e5]  # require library-scale columns (drops fold-change / p-value tables)
    if X.shape[1] == 0 or X.shape[0] < 5000: return None
    lcpm = np.log1p(X / X.sum() * 1e6).mean(axis=1)
    G = g.drop_duplicates(key).set_index(key)[tissues]
    common = lcpm.index.intersection(G.index)
    if len(common) < 5000: return None
    L = np.log1p(G.loc[common])
    var = L.var(axis=1); top = var.sort_values().index[-3000:]  # tissue-informative genes
    rho = np.array([spearmanr(lcpm[top], L.loc[top, t])[0] for t in tissues])
    return rho, X.shape[1], len(common)
def main(limit):
    g = load_gtex(); tissues = [c for c in g.columns if c not in ("Name", "Description", "ens")]
    cand = pd.read_csv("data/geo/candidates.csv").fillna("")
    s = cand.suppfile.str.upper()
    cand = cand[~s.str.contains("MTX|H5|RDS") & s.str.contains("TXT|TSV|CSV")]
    out = "results/geo_fidelity.csv"; log = "results/geo_fidelity_log.csv"
    done = set(pd.read_csv(out).gse) if os.path.exists(out) else set()
    tried = set(pd.read_csv(log).gse) if os.path.exists(log) else set()
    n = 0
    for _, r in cand.iterrows():
        if n >= limit: break
        if r.gse in done or r.gse in tried: continue
        oi = organ_of(r.title)
        if oi is None: continue
        n += 1; status = "ok"
        try:
            url, files = list_suppl(r.gse); res = None
            for f in files[:3]:
                try:
                    df = parse(fetch(url + f), f); m = match_genes(df, g)
                    if m is None: continue
                    res = score(m[0], m[1], g, tissues)
                    if res is not None: break
                except Exception as e: status = f"parse:{type(e).__name__}"
            if res is None:
                status = status if status != "ok" else "no_usable_matrix"
            else:
                rho, ns, ng = res
                match = ORGANS[oi][1]; mi = [tissues.index(t) for t in match]
                order = np.argsort(-rho); best = tissues[order[0]]
                rank = int(min(np.where(np.isin(order, mi))[0])) + 1
                z = float((rho[mi].max() - rho.mean()) / rho.std())
                row = {"gse": r.gse, "title": r.title, "organ": ORGANS[oi][1][0], "file": f, "n_samples": ns, "n_genes": ng,
                       "rho_match": float(rho[mi].max()), "best_tissue": best, "rho_best": float(rho.max()), "rank_match": rank, "z_match": z, "pdat": r.pdat}
                new = not os.path.exists(out)
                with open(out, "a", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=list(row)); new and w.writeheader(); w.writerow(row)
        except Exception as e: status = f"list:{type(e).__name__}"
        newl = not os.path.exists(log)
        with open(log, "a", newline="") as fh:
            w = csv.writer(fh); newl and w.writerow(["gse", "status"]); w.writerow([r.gse, status])
        print(r.gse, status, flush=True); time.sleep(0.3)
if __name__ == "__main__": main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
