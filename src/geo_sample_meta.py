"""Sample-level annotation of the clean-subset GEO series with GEOparse (SOFT family files).
Per series: number of samples, fraction whose title/source/characteristics mention organoids, fraction that look like non-organoid
material (tissue, biopsy, 2D/monolayer, cell line, fibroblast) without an organoid mention. Resumable. Output: results/geo_sample_meta.csv."""
import os, re, shutil, sys, logging
import pandas as pd, GEOparse
logging.getLogger("GEOparse").setLevel(logging.ERROR)
a = pd.read_csv("results/geo_label_audit.csv"); clean = a.gse[a.label_consistent & a.label_unambiguous & a.organoid_in_summary].tolist()
out = "results/geo_sample_meta.csv"; done = set(pd.read_csv(out).gse) if os.path.exists(out) else set()
NON = re.compile(r"tissue|biops|\b2d\b|monolayer|cell line|fibroblast|primary cell|ipsc(?!.*organoid)|blood|pbmc", re.I)
for gse in clean:
    if gse in done: continue
    d = f"/tmp/geo_{gse}"; os.makedirs(d, exist_ok=True)
    try:
        g = GEOparse.get_GEO(geo=gse, destdir=d, silent=True, how="full")
        n = org = non = 0
        for gsm in g.gsms.values():
            m = gsm.metadata; txt = " ".join(m.get("title", []) + m.get("source_name_ch1", []) + m.get("characteristics_ch1", []))
            n += 1; o = "organoid" in txt.lower() or "enteroid" in txt.lower() or "colonoid" in txt.lower() or "spheroid" in txt.lower()
            org += o; non += (not o) and bool(NON.search(txt))
        row = {"gse": gse, "n_gsm": n, "frac_organoid": org / n if n else None, "frac_non_organoid": non / n if n else None, "status": "ok"}
    except Exception as e: row = {"gse": gse, "n_gsm": None, "frac_organoid": None, "frac_non_organoid": None, "status": type(e).__name__}
    shutil.rmtree(d, ignore_errors=True)
    pd.DataFrame([row]).to_csv(out, mode="a", header=not os.path.exists(out), index=False); print(row, flush=True)
