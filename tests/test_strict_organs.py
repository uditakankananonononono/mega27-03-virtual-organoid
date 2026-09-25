import json, pandas as pd
def test_strict_organ_outputs():
    r = pd.read_csv("results/strict_organ_recovery.csv"); assert r.n.sum() == json.load(open("results/strict_organ_summary.json"))["n_strict"]
    D = pd.read_csv("results/strict_organ_deficits.csv"); assert (D.groupby("organ").size() == 100).all()
    p = json.load(open("results/strict_organ_parenchymal.json")); assert 0 < p["fisher_unfaithful_vs_faithful"]["p"] <= 1
def test_protocol_outputs():
    P = pd.read_csv("results/strict_protocol.csv"); assert set(P.derivation) <= {"PSC", "adult", "unclear"} and P.gse.is_unique
    d = json.load(open("results/strict_protocol.json")); assert 0 < d["cmh_organ_effect_given_derivation"]["p"] <= 1
def test_uniprot_outputs():
    U = pd.read_csv("results/strict_deficit_uniprot.csv"); assert ((U.share >= 0) & (U.share <= 1)).all() and (U.q_bh >= U.p - 1e-12).all()
def test_ae_replication():
    d = json.load(open("results/ae_replication.json")); D = pd.read_csv("results/ae_replication.csv")
    assert d["n_scored"] == int((D.status == "ok").sum()) and 0 < d["mwu_one_sided_p"] <= 1
def test_gsea_outputs():
    D = pd.read_csv("results/strict_gsea_hallmark.csv"); assert D.organ.nunique() >= 4 and D.NES.abs().max() < 10
