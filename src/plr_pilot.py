"""Toy-MC Neyman construction pilot for the binned 2D PLR at one WIMP mass.

1. Critical values: at each mu on a geometric grid (signal events in the full ROI),
   toys drawn at (mu, Table-I theta) with fluctuated global observables; q_crit(mu)
   is the 90% quantile of the two-sided q_mu.
2. Background-only toys: q_mu(toy) on the same grid; the upper limit is the largest
   crossing of q_mu - q_crit (linear interpolation). Median and +-1/2 sigma quantiles
   of the toy limits give the sensitivity; the power constraint lifts limits below
   the -1 sigma quantile to it (applied_1sigma).

Read-only on the archives; writes one JSON to --output (scratch). No results/ writes.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

import h5py
import numpy as np
from scipy.stats import chi2, norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from plr_core import BinnedLikelihood, neyman_upper_limit
from validate_backgrounds import TOTAL_SOURCES

REFERENCE_SIGMA_CM2 = 1e-36
CL = 0.9
QUANT = {"minus2sigma": norm.cdf(-2), "minus1sigma": norm.cdf(-1), "median": 0.5,
         "plus1sigma": norm.cdf(1), "plus2sigma": norm.cdf(2)}

_PATHS = None
_LKS = {}


def build(signal_path, background_path, mass_key):
    with h5py.File(signal_path, "r") as f:
        g = f[f"ws2024/{mass_key}"]
        if g.attrs["reference_sigma_cm2"] != REFERENCE_SIGMA_CM2:
            raise ValueError("unexpected reference sigma")
        s = g["hist_counts"][:]
        mass = float(g.attrs["mass_gev"])
    shapes, means, sigmas = [], [], []
    with h5py.File(background_path, "r") as f:
        for name in TOTAL_SOURCES:
            g = f[name]
            h = g["hist_counts"][:]
            shapes.append((h if h.sum() > 0 else g["hist_mc_counts"][:].astype(float)).ravel())
            means.append(float(g.attrs["table1_mean"]))
            sigmas.append(float(g.attrs["table1_sigma"]))
    lk = BinnedLikelihood(s, np.array(shapes), means, sigmas, mu_scale=5.0)
    return lk, float(s.sum()), mass


def _init(signal_path, background_path):
    global _PATHS
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    _PATHS = (signal_path, background_path)


def _lk(mass_key):
    if mass_key not in _LKS:
        _LKS.clear()
        _LKS[mass_key] = build(*_PATHS, mass_key)[0]
    return _LKS[mass_key]


def _signal_task(args):
    mass_key, mu, n_toys, seed = args
    _LK = _lk(mass_key)
    rng = np.random.default_rng(seed)
    q, mu_hat, fails = np.empty(n_toys), np.empty(n_toys), 0
    for t in range(n_toys):
        n, a = _LK.toy(rng, mu)
        qq, free = _LK.q_mu(n, mu, a)
        q[t], mu_hat[t] = qq, free.mu
        fails += not free.success
    return mu, q, mu_hat, fails


def _bkg_task(args):
    mass_key, mu_grid, n_toys, seed = args
    _LK = _lk(mass_key)
    rng = np.random.default_rng(seed)
    out = np.empty((n_toys, len(mu_grid)))
    mu_hat, fails = np.empty(n_toys), 0
    for t in range(n_toys):
        n, a = _LK.toy(rng, 0.0)
        free = _LK.fit(n, a)
        mu_hat[t] = free.mu
        fails += not free.success
        for k, mu in enumerate(mu_grid):
            out[t, k] = _LK.q_mu(n, mu, a, free=free)[0]
    return out, mu_hat, fails


def chunks(total, size):
    return [min(size, total - i) for i in range(0, total, size)]


def run_mass(pool, args, mass_key, mass_index):
    """Full construction at one mass using an existing worker pool; returns the report dict.

    mass_index=None seeds with args.seed alone (the 2026-11-18 pilot's stream)."""
    t0 = time.perf_counter()
    lk, s_ref, mass = build(args.signal, args.background, mass_key)
    to_sigma = REFERENCE_SIGMA_CM2 / s_ref
    mu_grid = np.geomspace(args.mu_min, args.mu_max, args.mu_points)

    entropy = args.seed if mass_index is None else [args.seed, mass_index]
    ss = np.random.SeedSequence(entropy)
    sig_jobs = []
    for mu in mu_grid:
        for n in chunks(args.toys_per_mu, 250):
            sig_jobs.append((mass_key, float(mu), n, ss.spawn(1)[0]))
    bkg_jobs = [(mass_key, mu_grid, n, ss.spawn(1)[0]) for n in chunks(args.bkg_toys, 25)]

    sig_res = pool.map(_signal_task, sig_jobs, chunksize=1)
    t_sig = time.perf_counter()
    bkg_res = pool.map(_bkg_task, bkg_jobs, chunksize=1)
    t_bkg = time.perf_counter()

    q_by_mu = {float(mu): [] for mu in mu_grid}
    muhat_sig, sig_fails = {float(mu): [] for mu in mu_grid}, 0
    for mu, q, mh, f in sig_res:
        q_by_mu[mu].append(q)
        muhat_sig[mu].append(mh)
        sig_fails += f
    crit, crit_err, frac_zero = [], [], []
    for mu in mu_grid:
        q = np.concatenate(q_by_mu[float(mu)])
        c = float(np.quantile(q, CL))
        # binomial order-statistic 68% interval on the quantile
        n = q.size
        lo_i = int(max(0, math.floor(n * CL - math.sqrt(n * CL * (1 - CL)))))
        hi_i = int(min(n - 1, math.ceil(n * CL + math.sqrt(n * CL * (1 - CL)))))
        qs = np.sort(q)
        crit.append(c)
        crit_err.append(0.5 * float(qs[hi_i] - qs[lo_i]))
        frac_zero.append(float(np.mean(np.concatenate(muhat_sig[float(mu)]) == 0.0)))
    crit = np.array(crit)

    qb = np.concatenate([r[0] for r in bkg_res])
    muhat_b = np.concatenate([r[1] for r in bkg_res])
    bkg_fails = sum(r[2] for r in bkg_res)
    limits = np.array([neyman_upper_limit(mu_grid, row, crit) for row in qb])
    n_beyond = int(np.sum(~np.isfinite(limits)))
    n_below = int(np.sum(limits <= mu_grid[0]))
    lim_q = {k: float(np.quantile(limits, p)) for k, p in QUANT.items()}
    floor = lim_q["minus1sigma"]
    pc = np.maximum(limits, floor)
    pc_q = {k: float(np.quantile(pc, p)) for k, p in QUANT.items()}
    # median precision: bootstrap
    brng = np.random.default_rng(1)
    boot = np.array([np.median(brng.choice(limits, limits.size)) for _ in range(2000)])
    # asymptotic reference: two-sided bounded q has 90% point between chi2_1 (2.706) values
    report = {
        "scope": "PLR toy-MC Neyman construction at one mass, WS2024 3.3 t-y templates; scratch output, not a declared result.",
        "mass_key": mass_key,
        "mass_gev": mass, "signal_events_at_reference_sigma": s_ref,
        "analysis_bins": int(lk.support.sum()),
        "settings": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}, "seed_entropy": entropy,
        "mu_grid_events": mu_grid.tolist(),
        "q_crit": crit.tolist(), "q_crit_order_stat_halfwidth": crit_err,
        "q_crit_asymptotic_chi2_1": float(chi2.ppf(CL, 1)),
        "signal_toy_fraction_mu_hat_zero": frac_zero,
        "fit_failures": {"signal_toys": int(sig_fails), "bkg_toys": int(bkg_fails)},
        "bkg_toys_mu_hat_zero_fraction": float(np.mean(muhat_b == 0.0)),
        "bkg_limits_beyond_grid": n_beyond, "bkg_limits_at_or_below_grid_start": n_below,
        "limit_quantiles_events": lim_q,
        "limit_quantiles_cm2": {k: v * to_sigma for k, v in lim_q.items()},
        "pc_limit_quantiles_cm2": {k: v * to_sigma for k, v in pc_q.items()},
        "median_bootstrap_sigma_cm2": float(boot.std() * to_sigma),
        "limits_events": limits.tolist(),
        "timing_s": {"signal_toys": t_sig - t0, "bkg_toys": t_bkg - t_sig, "total": time.perf_counter() - t0},
        "workers": args.workers,
        "inputs": {"signal": str(args.signal), "background": str(args.background)},
    }
    print(f"q_crit: " + " ".join(f"{m:.2f}:{c:.2f}" for m, c in zip(mu_grid, crit)))
    print(f"limit quantiles (events): " + ", ".join(f"{k} {v:.3f}" for k, v in lim_q.items()))
    print(f"PC median {pc_q['median'] * to_sigma:.4g} cm2 (+-{boot.std() * to_sigma:.2g} boot); "
          f"-1/+1 sigma {pc_q['minus1sigma'] * to_sigma:.3g}/{pc_q['plus1sigma'] * to_sigma:.3g}; "
          f"beyond grid {n_beyond}, at grid start {n_below}; fails {sig_fails}/{bkg_fails}")
    print(f"timing {report['timing_s']}")
    return report


def add_common_args(ap):
    ap.add_argument("--signal", required=True, type=Path)
    ap.add_argument("--background", required=True, type=Path)
    ap.add_argument("--mu-min", type=float, default=0.5)
    ap.add_argument("--mu-max", type=float, default=20.0)
    ap.add_argument("--mu-points", type=int, default=20)
    ap.add_argument("--toys-per-mu", type=int, default=1000)
    ap.add_argument("--bkg-toys", type=int, default=500)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--seed", type=int, default=20261118)


def check_scratch(path):
    if not path.is_absolute() or str(path.resolve()).startswith(str(HERE.parent / "results")):
        raise ValueError("output must be an absolute scratch path outside results/")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    add_common_args(ap)
    ap.add_argument("--mass-key", default="mass_40gev")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    check_scratch(args.output)
    ctx = mp.get_context("spawn")
    with ctx.Pool(args.workers, initializer=_init, initargs=(args.signal, args.background)) as pool:
        report = run_mass(pool, args, args.mass_key, mass_index=None)
    args.output.write_text(json.dumps(report, indent=1) + "\n")


if __name__ == "__main__":
    main()
