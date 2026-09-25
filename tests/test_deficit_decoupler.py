import sys, pandas as pd
sys.path.insert(0, "src")
from deficit_decoupler import top_frac_hits

def test_top_frac_hits():
    a = pd.Series(range(100), index=[f"T{i}" for i in range(100)], dtype=float)
    assert top_frac_hits(a, ["T99", "T95", "T0", "MISSING"]) == (2, 3)
