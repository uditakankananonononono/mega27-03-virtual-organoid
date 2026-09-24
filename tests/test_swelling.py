import numpy as np
import pandas as pd
import pytest

from vorganoid.swelling import cluster_bootstrap_slope, fit_alpha, predict_fold, size_slope


def test_alpha_one_is_size_invariant():
    A0 = np.array([100.0, 1000.0, 10000.0])
    s = predict_fold(A0, 1.0, 0.6)
    assert np.allclose(s, s[0])


def test_surface_limited_gives_negative_slope():
    A0 = np.geomspace(100, 20000, 50)
    s = predict_fold(A0, 2 / 3, 5.0)
    assert size_slope(A0, s) < 0


def test_superlinear_gives_positive_slope():
    A0 = np.geomspace(100, 20000, 50)
    s = predict_fold(A0, 1.3, -1e-3)
    assert size_slope(A0, s) > 0


def test_fold_is_one_when_no_secretion():
    A0 = np.geomspace(100, 5000, 10)
    assert np.allclose(predict_fold(A0, 2 / 3, 0.0), 1.0)


@pytest.mark.parametrize("alpha", [0.7, 1.0, 1.3])
def test_fit_alpha_recovers_truth(alpha):
    rng = np.random.default_rng(1)
    A0 = np.geomspace(200, 20000, 400)
    c = {0.7: 3.0, 1.0: 0.8, 1.3: -2e-2}[alpha]
    s = predict_fold(A0, alpha, c) * np.exp(rng.normal(0, 0.02, len(A0)))
    fit = fit_alpha(A0, s)
    assert abs(fit["alpha"] - alpha) <= 0.1


def test_cluster_bootstrap_ci_contains_estimate():
    rng = np.random.default_rng(0)
    A0 = np.geomspace(200, 20000, 300)
    df = pd.DataFrame({"A0": A0, "swelling": predict_fold(A0, 1.2, -0.01) * np.exp(rng.normal(0, 0.05, 300)),
                       "well": np.repeat(np.arange(30), 10)})
    est, lo, hi = cluster_bootstrap_slope(df, n_boot=200)
    assert lo <= est <= hi and est > 0
