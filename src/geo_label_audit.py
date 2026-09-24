"""Audit title-derived organ labels against GSE summaries (NCBI esummary). Output: results/geo_label_audit.csv."""
import json, re, sys, time, urllib.request, urllib.parse
import pandas as pd
sys.path.insert(0, "src"); from geo_fidelity import ORGANS
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
d = pd.read_csv("results/geo_fidelity.csv")
term = " OR ".join(f"{g}[ACCN]" for g in d.gse) 
rows = []
for i in range(0, len(d), 50):
    t = " OR ".join(f"{g}[ACCN]" for g in d.gse[i:i+50]) + " AND gse[etyp]"
    ids = json.load(urllib.request.urlopen(E + "esearch.fcgi?" + urllib.parse.urlencode({"db": "gds", "term": t, "retmax": 200, "retmode": "json"})))["esearchresult"]["idlist"]
    r = json.load(urllib.request.urlopen(E + "esummary.fcgi?" + urllib.parse.urlencode({"db": "gds", "id": ",".join(ids), "retmode": "json"})))["result"]
    for u in r["uids"]:
        x = r[u]; rows.append({"gse": "GSE" + x["gse"], "summary": x["summary"]})
    time.sleep(0.4)
S = pd.DataFrame(rows).drop_duplicates("gse")
d = d.merge(S, on="gse", how="left")
def organs(text): return sorted({ORGANS[i][1][0] for i, (rx, _) in enumerate(ORGANS) if re.search(rx, str(text), re.I)})
d["summary_organs"] = d.summary.map(lambda s: ";".join(organs(s)))
d["organoid_in_summary"] = d.summary.str.contains("organoid", case=False, na=False)
d["label_consistent"] = [o in so.split(";") for o, so in zip(d.organ, d.summary_organs)]
d["label_unambiguous"] = d.summary_organs.map(lambda s: len([x for x in s.split(";") if x]) == 1)
d[["gse", "organ", "summary_organs", "organoid_in_summary", "label_consistent", "label_unambiguous", "summary"]].to_csv("results/geo_label_audit.csv", index=False)
print(len(d), d.label_consistent.sum(), (d.label_consistent & d.label_unambiguous & d.organoid_in_summary).sum())
