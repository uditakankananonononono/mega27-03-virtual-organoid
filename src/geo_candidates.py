"""List human organoid RNA-seq GEO series with supplementary files (NCBI E-utilities)."""
import json, urllib.request, urllib.parse, csv, time, sys
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
term = 'organoid[All] AND "Homo sapiens"[orgn] AND gse[etyp] AND "expression profiling by high throughput sequencing"[gdstype]'
def get(url):
    for k in range(4):
        try: return json.load(urllib.request.urlopen(url, timeout=60))
        except Exception as e: time.sleep(2 + 3 * k)
    raise RuntimeError(url)
ids = []
for start in range(0, 3000, 1000):
    r = get(E + "esearch.fcgi?" + urllib.parse.urlencode({"db": "gds", "term": term, "retmax": 1000, "retstart": start, "retmode": "json"}))
    ids += r["esearchresult"]["idlist"]
rows = []
for i in range(0, len(ids), 200):
    r = get(E + "esummary.fcgi?" + urllib.parse.urlencode({"db": "gds", "id": ",".join(ids[i:i+200]), "retmode": "json"}))["result"]
    for u in r["uids"]:
        d = r[u]
        rows.append({"gse": "GSE" + d["gse"], "title": d["title"], "n_samples": d["n_samples"], "suppfile": d.get("suppfile", ""),
                     "pdat": d["pdat"], "ftplink": d.get("ftplink", "")})
    time.sleep(0.4)
with open("data/geo/candidates.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(len(ids), len(rows))
