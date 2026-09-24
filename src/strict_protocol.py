"""Classify each strict organoid series as PSC-derived (iPSC/ESC/hPSC/pluripotent) or adult/tissue-derived from GEO series and GSM
text (GEOparse), then test whether derivation explains organ recovery. Output: results/strict_protocol.csv, results/strict_protocol.json"""
import os, re, json, logging, pandas as pd, GEOparse
from scipy.stats import fisher_exact
import statsmodels.formula.api as smf
logging.getLogger("GEOparse").setLevel(logging.ERROR)
PSC = re.compile(r"\b(i?psc|hipsc|hesc|esc-derived|pluripotent|embryonic stem|induced pluripotent|h9|h1 cells|wtc-?11|directed differentiation)\b", re.I)
ADULT = re.compile(r"\b(adult stem|lgr5|biops(y|ies)|patient-derived|tissue-derived|primary tissue|resection|crypts?|donor tissue|surgical)\b", re.I)
m = pd.read_csv("results/geo_sample_meta.csv"); strict = list(m.gse[m.frac_organoid == 1])
out_csv = "results/strict_protocol.csv"; done = set(pd.read_csv(out_csv).gse) if os.path.exists(out_csv) else set()
for gse in strict:
    if gse in done: continue
    try:
        g = GEOparse.get_GEO(geo=gse, destdir="/tmp/geosoft", silent=True)
        txt = " ".join(" ".join(v) for k, v in g.metadata.items() if k in ("title", "summary", "overall_design"))
        txt += " " + " ".join(" ".join(" ".join(x.metadata.get(k, [])) for k in ("title", "source_name_ch1", "characteristics_ch1", "growth_protocol_ch1", "extract_protocol_ch1")) for x in g.gsms.values())
        np_, na = len(PSC.findall(txt)), len(ADULT.findall(txt))
        row = {"gse": gse, "psc_hits": np_, "adult_hits": na, "derivation": "PSC" if np_ > na else ("adult" if na > np_ else "unclear"), "status": "ok"}
    except Exception as e:
        row = {"gse": gse, "psc_hits": 0, "adult_hits": 0, "derivation": "unclear", "status": f"err:{type(e).__name__}"}
    pd.DataFrame([row]).to_csv(out_csv, mode="a", header=not os.path.exists(out_csv), index=False)
    for f in os.listdir("/tmp/geosoft"): os.remove(os.path.join("/tmp/geosoft", f))
    print(row, flush=True)
P = pd.read_csv(out_csv); R = pd.read_csv("results/geo_fidelity_methods_ranks.csv"); R = R[R.method == "SALL_PARONLY"]
D = P.merge(R[["gse", "organ", "rank_match"]], on="gse"); D["top1"] = (D.rank_match == 1).astype(int)
D["group"] = D.organ.map(lambda o: "faithful_organ" if o in ("Colon - Transverse", "Brain - Cortex") else ("unfaithful_organ" if o in ("Liver", "Kidney - Cortex", "Lung") else "other"))
res = {"counts": D.groupby(["organ", "derivation"]).size().unstack(fill_value=0).to_dict(),
       "top1_by_derivation": D.groupby("derivation").top1.agg(["mean", "size"]).to_dict()}
d2 = D[D.derivation != "unclear"]
t = [[int(((d2.derivation == "PSC") & (d2.top1 == 1)).sum()), int(((d2.derivation == "PSC") & (d2.top1 == 0)).sum())],
     [int(((d2.derivation == "adult") & (d2.top1 == 1)).sum()), int(((d2.derivation == "adult") & (d2.top1 == 0)).sum())]]
res["fisher_psc_vs_adult_top1"] = {"table": t, "p": float(fisher_exact(t)[1])}
for grp in ["faithful_organ", "unfaithful_organ"]:
    s = d2[d2.group == grp]; res[f"within_{grp}"] = s.groupby("derivation").top1.agg(["mean", "size"]).to_dict()
# Organ effect within each derivation stratum (Fisher) and pooled Cochran-Mantel-Haenszel; a logistic GLM separates perfectly here.
from statsmodels.stats.contingency_tables import StratifiedTable
tabs = []
for der in ["PSC", "adult"]:
    s_ = d2[(d2.derivation == der) & (d2.group != "other")]
    tt = [[int(((s_.group == "faithful_organ") & (s_.top1 == 1)).sum()), int(((s_.group == "faithful_organ") & (s_.top1 == 0)).sum())],
          [int(((s_.group == "unfaithful_organ") & (s_.top1 == 1)).sum()), int(((s_.group == "unfaithful_organ") & (s_.top1 == 0)).sum())]]
    res[f"organ_effect_within_{der}"] = {"table": tt, "fisher_p": float(fisher_exact(tt)[1])}; tabs.append(tt)
st = StratifiedTable([[[a + 0.5 for a in r] for r in t_] for t_ in tabs])
res["cmh_organ_effect_given_derivation"] = {"pooled_OR_haldane": float(st.oddsratio_pooled), "p": float(st.test_null_odds(correction=True).pvalue), "note": "0.5 Haldane correction for zero cells"}
json.dump(res, open("results/strict_protocol.json", "w"), indent=1, default=str); print(json.dumps(res, indent=1, default=str))
