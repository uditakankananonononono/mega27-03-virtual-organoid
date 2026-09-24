"""Cell-identity and pathway annotation of culture-fibroblast attractor driver genes with two independent services:
Enrichr (maayanlab.cloud; libraries CellMarker_2024, PanglaoDB_Augmented_2021, MSigDB_Hallmark_2020) and the
Reactome AnalysisService (reactome.org/AnalysisService, over-representation, projected to human).
Input: results/attractor_genes.csv (top-200 by d, and ranks 801-1000 as a within-file contrast). Output: results/attractor_enrichr.csv, results/attractor_reactome.csv."""
import json, time, urllib.request, urllib.parse, pandas as pd
G = pd.read_csv("results/attractor_genes.csv").dropna(subset=["symbol"]).sort_values("d")
sets = {"top200": list(G.symbol[::-1][:200]), "rank801_1000": list(G.symbol[:200])}  # file holds only top-1000 d>0 genes; the second set is a within-file contrast, not down-regulated genes
LIBS = ["CellMarker_2024", "PanglaoDB_Augmented_2021", "MSigDB_Hallmark_2020"]
def enrichr(genes):
    b = "----b"; body = (f"--{b}\r\nContent-Disposition: form-data; name=\"list\"\r\n\r\n" + "\n".join(genes) + f"\r\n--{b}\r\nContent-Disposition: form-data; name=\"description\"\r\n\r\nattractor\r\n--{b}--\r\n").encode()
    r = urllib.request.Request("https://maayanlab.cloud/Enrichr/addList", data=body, headers={"Content-Type": f"multipart/form-data; boundary={b}"})
    return json.load(urllib.request.urlopen(r, timeout=60))["userListId"]
rows, rrows = [], []
for dirn, genes in sets.items():
    uid = enrichr(genes)
    for lib in LIBS:
        res = json.load(urllib.request.urlopen(f"https://maayanlab.cloud/Enrichr/enrich?userListId={uid}&backgroundType={lib}", timeout=60))[lib]
        for x in res[:10]: rows.append({"direction": dirn, "library": lib, "term": x[1], "p": x[2], "combined_score": x[4], "q_bh": x[6], "overlap": ";".join(x[5])})
        time.sleep(0.5)
    r = urllib.request.Request("https://reactome.org/AnalysisService/identifiers/projection?pageSize=15&page=1", data="\n".join(genes).encode(), headers={"Content-Type": "text/plain", "User-Agent": "mega27-vorganoid/1.0", "Accept": "application/json"})
    for p in json.load(urllib.request.urlopen(r, timeout=90))["pathways"]:
        rrows.append({"direction": dirn, "stId": p["stId"], "name": p["name"], "found": p["entities"]["found"], "total": p["entities"]["total"], "p": p["entities"]["pValue"], "fdr": p["entities"]["fdr"]})
pd.DataFrame(rows).to_csv("results/attractor_enrichr.csv", index=False); pd.DataFrame(rrows).to_csv("results/attractor_reactome.csv", index=False)
E = pd.DataFrame(rows); R = pd.DataFrame(rrows)
for d in sets:
    for lib in LIBS: print(d, lib, E[(E.direction == d) & (E.library == lib)].head(3)[["term", "q_bh"]].values.tolist())
    print(d, "Reactome", R[R.direction == d].head(4)[["name", "fdr"]].values.tolist())
