"""Biophysical swelling twin for CFTR-driven organoid lumen expansion.

Model family (volume V ∝ A^{3/2} for a near-spherical organoid of projected area A):
    dV/dt = k V^alpha
    alpha = 2/3 : secretion limited by epithelial surface (constant flux per area)
    alpha = 1   : secretion proportional to volume/cell mass (size-invariant fold change)
    alpha > 1   : super-linear; larger organoids swell disproportionately more
Closed form for alpha != 1:  V1^{1-alpha} = V0^{1-alpha} + (1-alpha) k t
The observable is the area fold change s = A1/A0 over the assay window.

Size-dependence diagnostic: slope b of log s on log A0.
    alpha < 1 -> b < 0 ; alpha = 1 -> b = 0 ; alpha > 1 -> b > 0
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar


def predict_fold(A0: np.ndarray, alpha: float, c: float) -> np.ndarray:
    V0 = np.asarray(A0, float) ** 1.5
    if abs(alpha - 1.0) < 1e-6:
        V1 = V0 * np.exp(c)
    else:
        base = V0 ** (1 - alpha) + c
        base = np.maximum(base, 1e-12)
        V1 = base ** (1 / (1 - alpha))
    return (V1 / V0) ** (2.0 / 3.0)


def fit_alpha(A0: np.ndarray, s: np.ndarray, alphas=None) -> dict:
    """Least squares in log fold change over alpha grid, c profiled out."""
    A0 = np.asarray(A0, float); ls = np.log(np.asarray(s, float))
    alphas = np.linspace(0.0, 2.5, 51) if alphas is None else alphas
    best = None
    for a in alphas:
        V0 = A0 ** 1.5
        if abs(a - 1) < 1e-6:
            c0, lo, hi = np.log(np.median(s) ** 1.5), -10, 10
        else:
            scale = np.median(V0 ** (1 - a))
            c0, lo, hi = 0.0, -0.999 * scale if a < 1 else -10 * scale, 10 * scale
        f = lambda c: np.mean((np.log(predict_fold(A0, a, c)) - ls) ** 2)
        res = minimize_scalar(f, bounds=(lo, hi), method="bounded")
        if best is None or res.fun < best["mse"]:
            best = {"alpha": float(a), "c": float(res.x), "mse": float(res.fun)}
    best["mse_null"] = float(np.var(ls))  # size-invariant (alpha=1) equivalent
    return best


def size_slope(A0, s) -> float:
    return float(np.polyfit(np.log(A0), np.log(s), 1)[0])


def cluster_bootstrap_slope(df: pd.DataFrame, cluster: str = "well", n_boot: int = 1000,
                            seed: int = 0) -> tuple[float, float, float]:
    """Slope with CI from resampling wells (organoids in a well are not independent)."""
    rng = np.random.default_rng(seed)
    groups = [g for _, g in df.groupby(cluster)]
    est = size_slope(df.A0, df.swelling)
    boots = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        d = pd.concat([groups[i] for i in pick])
        boots.append(size_slope(d.A0, d.swelling))
    return est, float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))
