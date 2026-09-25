"""Scan EBI BioStudies/ArrayExpress for human organoid RNA-seq studies with processed expression files.
Writes data/ae/ae_candidates.csv (accession, title, organ guess, processed file URLs)."""
import json, re, urllib.request, urllib.parse, pandas as pd, time
Q = ["liver organoid", "hepatic organoid", "hepatocyte organoid", "cholangiocyte organoid", "kidney organoid", "tubuloid", "nephron organoid",
     "lung organoid", "alveolar organoid", "airway organoid", "bronchial organoid", "intestinal organoid", "colon organoid", "colonic organoid",
     "cerebral organoid", "cortical organoid", "brain organoid", "gastric organoid", "pancreas organoid"]
def get(u, t=30): return urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "mega27/1.0"}), timeout=t).read().decode("utf8", "ignore")
seen, rows = set(), []
for q in Q:
    for page in (1, 2):
        try: d = json.loads(get(f"https://www.ebi.ac.uk/biostudies/api/v1/arrayexpress/search?query={urllib.parse.quote(q)}&pageSize=50&page={page}"))
        except Exception: continue
        for h in d.get("hits", []):
            a = h["accession"]
            if a in seen or not a.startswith("E-MTAB"): continue
            seen.add(a); t = h.get("title", "")
            if re.search(r"single.?cell|scRNA|snRNA|ATAC|ChIP|CUT&|methyl|mouse|murine|CLIP|10x", t, re.I): continue
            try:
                ftp = json.loads(get(f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{a}/info")).get("ftpLink", "").replace("ftp://", "https://")
                files = re.findall(r'href="([^"/?][^"]*)"', get(ftp + "/Files/"))
            except Exception: continue
            proc = [f for f in files if re.search(r"count|tpm|fpkm|rpkm|expr|matrix|ReadsPerGene|processed|\.txt$|\.tsv|\.csv", f, re.I) and not re.search(r"fastq|sdrf|idf|md5|bam|cram", f, re.I)]
            if proc: rows.append({"accession": a, "title": t, "query": q, "n_files": len(proc), "files": ";".join(ftp + "/Files/" + f for f in proc[:60])}); print(a, q, len(proc), t[:80], flush=True)
            time.sleep(0.2)
pd.DataFrame(rows).to_csv("data/ae/ae_candidates.csv", index=False)
