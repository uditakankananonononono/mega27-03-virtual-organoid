import pandas as pd
def test_annotation_outputs():
    E = pd.read_csv("results/attractor_enrichr.csv"); R = pd.read_csv("results/attractor_reactome.csv")
    assert set(E.direction) == {"top200", "rank801_1000"} and (E.q_bh >= 0).all() and (R.fdr <= 1).all()
