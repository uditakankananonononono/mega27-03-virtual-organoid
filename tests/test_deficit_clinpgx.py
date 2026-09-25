import sys, numpy as np, pandas as pd
sys.path.insert(0, "src")
from deficit_clinpgx import core_pharmacogenes, h2_stat

def test_core_pharmacogenes_cpic_only():
    df = pd.DataFrame({"Ensembl Id": ["E1", "E2", None], "Is VIP": ["Yes", "Yes", "Yes"], "Has CPIC Dosing Guideline": ["Yes", "No", "Yes"]})
    assert core_pharmacogenes(df) == {"E1"}
    assert core_pharmacogenes(df, use_vip=True) == {"E1", "E2"}

def test_h2_stat():
    Z = np.array([[3.0, 1.0, 0.0], [3.0, 1.0, 2.0]])
    assert h2_stat(Z, ["Liver", "A", "B"]) == 2.0
