import json, pandas as pd
def test_strict_json_consistent():
    d = json.load(open("results/geo_fidelity_strict.json")); m = pd.read_csv("results/geo_sample_meta.csv")
    assert d["n_strict"] == int((m.frac_organoid == 1).sum())
    for k, v in d.items():
        if isinstance(v, dict):
            assert 0 <= v["strict"]["fibroblast_best_frac"] <= 1 and v["strict"]["n"] <= d["n_strict"]
