import numpy as np, pandas as pd
from vorganoid.sizeaware import size_effects, attenuation
from vorganoid.cli import main


def _table(atten=True, seed=0):
    rng = np.random.default_rng(seed); rows = []
    for d in range(8):
        for cond in ("DMSO", "DRUG"):
            A0 = np.exp(rng.uniform(5, 9, 400))
            eff = (0.1 + 0.4 * (np.log(A0) - 5) / 4) if atten else 0.3
            ls = rng.normal(0.2 + (eff if cond == "DRUG" else 0), 0.1, size=400)
            rows += [{"donor": f"D{d}", "condition": cond, "A0": a, "swelling": np.exp(l)} for a, l in zip(A0, ls)]
    return pd.DataFrame(rows)


def test_detects_planted_attenuation():
    r = attenuation(size_effects(_table(True), "DRUG")[0])
    assert r["mean_attenuation"] > 0.2 and r["ci95"][0] > 0 and r["donors_positive"] == 8


def test_no_attenuation_when_size_invariant():
    r = attenuation(size_effects(_table(False), "DRUG")[0])
    assert abs(r["mean_attenuation"]) < 0.03


def test_cli(tmp_path, capsys):
    f = tmp_path / "t.csv"; _table(True).to_csv(f, index=False)
    assert main(["sizeaware", str(f), "--drug", "DRUG"]) == 0
    assert "mean_attenuation" in capsys.readouterr().out
