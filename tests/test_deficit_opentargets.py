import numpy as np, sys
sys.path.insert(0, "src")
from deficit_opentargets import diag_dominance, overlap_matrix

def test_diag_dominance_identity():
    assert diag_dominance(np.eye(3) * 3) == 3.0
    assert diag_dominance(np.ones((4, 4))) == 0.0

def test_overlap_matrix():
    L = {"a": {1, 2}, "b": {3}}; S = {"a": {2, 3}, "b": {3}}
    assert overlap_matrix(L, S, ["a", "b"]).tolist() == [[1, 0], [1, 1]]
