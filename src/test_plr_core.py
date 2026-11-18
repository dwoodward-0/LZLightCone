"""Unit tests for plr_core (plain asserts; run: lc run local -- python src/test_plr_core.py ...).

1. Analytic gradients vs central finite differences (free and conditional).
2. One-bin, no-nuisance exact Neyman construction with the class's q_mu reproduces
   the colleague FC_ints upper limits on FC's own 0.005-event grid.
3. Asimov data at the nominal model: free fit returns mu = 0, theta = Table I means;
   q_0 = 0.
4. Optional benchmark on the real 40 GeV/c2 templates (--signal/--background).
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

import numpy as np
from scipy.stats import poisson

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from plr_core import BinnedLikelihood

ALPHA = 0.1
FC_STEP = 0.005


def test_gradients(rng):
    nb, J = 40, 3
    shapes = rng.random((J, nb))
    lk = BinnedLikelihood(rng.random(nb), shapes, means=[50.0, 20.0, 5.0], sigmas=[5.0, 3.0, 1.0], mu_scale=10.0)
    n = rng.poisson(lk.expected(7.0, lk.means)).astype(float)
    a = lk.means + 0.3
    for mu_fixed, dim in ((None, J + 1), (4.0, J)):
        z = rng.random(dim) + 0.5
        f, g = lk._nll_grad(z, n, a, mu_fixed)
        h = 1e-6
        num = np.array([(lk._nll_grad(z + h * e, n, a, mu_fixed)[0] - lk._nll_grad(z - h * e, n, a, mu_fixed)[0])
                        / (2 * h) for e in np.eye(dim)])
        assert np.allclose(g, num, rtol=1e-5, atol=1e-6), (mu_fixed, g, num)
    print("PASS gradients")


def one_bin_upper(lk, n_obs, b, mu_lo, mu_hi):
    """Largest mu on FC's 0.005 grid (scanned over [mu_lo, mu_hi]) whose 90% LR-ordered acceptance set contains n_obs."""
    best = None
    q_cache = {}
    for mu in np.arange(mu_lo, mu_hi + FC_STEP / 2, FC_STEP):
        mu = round(float(mu), 6)
        tot = mu + b
        nmax = int(tot + 12 * np.sqrt(tot) + 15)
        ns = np.arange(nmax + 1)
        q = np.empty(nmax + 1)
        for k in ns:
            key = (int(k), mu)
            if key not in q_cache:
                q_cache[key] = lk.q_mu(np.array([float(k)]), mu)[0]
            q[k] = q_cache[key]
        p = poisson.pmf(ns, tot)
        order = np.lexsort((-p, q))  # ascending q (= descending LR); ties by probability
        cum = np.cumsum(p[order])
        accepted = order[: int(np.searchsorted(cum, 1 - ALPHA)) + 1]
        if n_obs in accepted:
            best = mu
    return best


def test_one_bin_fc(fc_module):
    cases = [(9, 3.7), (0, 0.6879), (3, 2.0), (1, 0.5)]
    for n_obs, b in cases:
        lk = BinnedLikelihood([1.0], np.zeros((0, 1)), means=[], sigmas=[], fixed_background=[b])
        lo_fc, hi_fc = fc_module.FC_ints(n_obs, b, alpha=ALPHA)
        hi = one_bin_upper(lk, n_obs, b, max(hi_fc - 0.3, 0.0), hi_fc + 0.3)
        print(f"  n={n_obs} b={b}: PLR-construction upper {hi:.3f} vs FC_ints {hi_fc:.3f}")
        assert abs(hi - hi_fc) <= FC_STEP + 1e-9, (n_obs, b, hi, hi_fc)
    print("PASS one-bin FC reproduction")


def test_asimov(rng):
    nb, J = 200, 4
    shapes = rng.random((J, nb)) ** 3
    sig = np.exp(-np.linspace(0, 5, nb))
    lk = BinnedLikelihood(sig, shapes, means=[300.0, 40.0, 5.0, 0.5], sigmas=[30.0, 4.0, 1.0, 0.2])
    asimov = lk.expected(0.0, lk.means)
    fr = lk.fit(asimov)
    assert fr.mu < 1e-3 and np.allclose(fr.theta, lk.means, rtol=1e-4, atol=1e-4), (fr.mu, fr.theta)
    q0, _ = lk.q_mu(asimov, 0.0)
    assert q0 < 1e-8, q0
    print("PASS Asimov recovery")


def benchmark(signal_path, background_path, n_toys, rng):
    import h5py
    sys.path.insert(0, str(HERE))
    from validate_backgrounds import TOTAL_SOURCES
    import csv
    with h5py.File(signal_path, "r") as f:
        s = f["ws2024/mass_40gev/hist_counts"][:]
    means, sigmas, shapes = [], [], []
    with h5py.File(background_path, "r") as f:
        for name in TOTAL_SOURCES:
            g = f[name]
            h = g["hist_counts"][:]
            shapes.append(h if h.sum() > 0 else g["hist_mc_counts"][:].astype(float))
            means.append(float(g.attrs["table1_mean"]))
            sigmas.append(float(g.attrs["table1_sigma"]))
    s_ref = float(s.sum())
    lk = BinnedLikelihood(s, np.array([x.ravel() for x in shapes]), means, sigmas, mu_scale=5.0)
    print(f"  analysis bins: {lk.support.sum()} of {lk.support.size}; signal events at 1e-36: {s_ref:.4g}")
    t_free = t_cond = 0.0
    nfev = []
    for _ in range(n_toys):
        n, a = lk.toy(rng, 0.0)
        t0 = time.perf_counter()
        free = lk.fit(n, a)
        t1 = time.perf_counter()
        for mu in (2.0, 5.0, 10.0):
            c = lk.fit(n, a, mu_fixed=mu, start=(mu, free.theta))
            nfev.append(c.nfev)
        t2 = time.perf_counter()
        t_free += t1 - t0
        t_cond += (t2 - t1) / 3
        assert free.success
    print(f"  per fit: free {1e3 * t_free / n_toys:.2f} ms, conditional {1e3 * t_cond / n_toys:.2f} ms "
          f"(median nfev {np.median(nfev):.0f})")
    return t_free / n_toys, t_cond / n_toys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fc-source", required=True, type=Path)
    ap.add_argument("--signal", type=Path)
    ap.add_argument("--background", type=Path)
    ap.add_argument("--bench-toys", type=int, default=50)
    args = ap.parse_args()
    sys.path.insert(0, str(args.fc_source.resolve().parent))
    import FeldmanCousins
    rng = np.random.default_rng(20261118)
    test_gradients(rng)
    test_asimov(rng)
    test_one_bin_fc(FeldmanCousins)
    if args.signal and args.background:
        benchmark(args.signal, args.background, args.bench_toys, rng)


if __name__ == "__main__":
    main()
