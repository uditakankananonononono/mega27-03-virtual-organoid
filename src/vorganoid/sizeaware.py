"""Size-aware FIS analysis: modulator-vs-control swelling effect as a function of starting organoid size.

Input: per-organoid table with columns donor, condition, A0 (starting area) and swelling (A1/A0), optionally
a forskolin column. For each donor, organoids are binned by starting-size quantiles computed on the pooled
table; the per-bin effect is Delta_k = mean log s (drug, bin k) - mean log s (control, bin k).
Attenuation = Delta_top - Delta_bottom; 95% CI by donor bootstrap.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def size_effects(df, drug, control="DMSO", n_bins=4, min_per_cell=5):
    d = df[df.condition.isin([drug, control])].copy()
    d = d[(d.swelling > 0) & (d.A0 > 0)]
    d["ls"] = np.log(d.swelling)
    edges = np.quantile(d.A0, np.linspace(0, 1, n_bins + 1)); edges[-1] += 1e-9
    d["bin"] = np.clip(np.searchsorted(edges, d.A0, side="right") - 1, 0, n_bins - 1)
    rows = []
    for don, g in d.groupby("donor"):
        eff = []
        for k in range(n_bins):
            a = g[(g.condition == drug) & (g.bin == k)].ls; c = g[(g.condition == control) & (g.bin == k)].ls
            eff.append(a.mean() - c.mean() if len(a) >= min_per_cell and len(c) >= min_per_cell else np.nan)
        rows.append([don] + eff)
    return pd.DataFrame(rows, columns=["donor"] + [f"bin{k}" for k in range(n_bins)]), edges


def attenuation(per_donor, n_boot=1000, seed=0):
    cols = [c for c in per_donor.columns if c.startswith("bin")]
    x = per_donor.dropna(subset=[cols[0], cols[-1]])
    att = (x[cols[-1]] - x[cols[0]]).values
    rng = np.random.default_rng(seed)
    boots = [att[rng.integers(0, len(att), len(att))].mean() for _ in range(n_boot)] if len(att) else [np.nan]
    return {"n_donors": int(len(att)), "mean_attenuation": float(np.mean(att)) if len(att) else float("nan"),
            "ci95": [float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))],
            "donors_positive": int((att > 0).sum())}
