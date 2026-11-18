"""Progressive-voxel FC counting-region scan over the WIMP mass grid.

For each mass, bins of the existing 77x91 (S1c, log10 S2c) physical
histograms are ranked by s_i/b_i (b_i = 0 ranks first, ties by s_i) and added
in that order (bins with s_i = 0 never). At each step k the exact
background-only median limit is FC_ints(n_med(B_k), B_k), with n_med the
Poisson(B_k) median (exact: FC upper limits are non-decreasing in n, asserted
per chosen region). sigma_med(k) = 1e-36 * mu_med / S_k; the chosen region
minimizes sigma_med. Bands are FC limits at Poisson count quantiles;
power-constrained values use floor = exp(Phi^-1(0.1)) x median.

Ranking uses raw s_i/b_i with no MC-sparsity guard (user decision 2026-11-17:
a first-pass counting cross-check). At 40 GeV/c2 the optimum sits where the
sparse heavy-weight ER MC (214Pb etc.) has no events; see AGENTS.md.

Read-only on unchanged archives; no NEST, no regeneration, no results/ writes.
Writes scan.json and fc_sensitivity_scan.csv into an existing scratch dir.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time
import warnings

import h5py
import numpy as np
from scipy.stats import norm, poisson

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from validate_backgrounds import TOTAL_SOURCES, load_band
from fc_counting_pilot import EXPECTED_FC_SHA256, sha256

REFERENCE_SIGMA_CM2 = 1e-36
ALPHA = 0.1
B_CAP = 50.0
NEAR_OPTIMUM = 1.05
QUANTILES = {"minus2sigma": float(norm.cdf(-2)), "minus1sigma": float(norm.cdf(-1)),
             "median": 0.5, "plus1sigma": float(norm.cdf(1)), "plus2sigma": float(norm.cdf(2))}


def poisson_quantile(q, b):
    """Smallest n with CDF(n; b) >= q (b = 0 is the point mass at 0)."""
    return 0 if b <= 0 else int(poisson.ppf(q, b))


class FC:
    def __init__(self, module):
        self.fc, self.cache = module, {}

    def upper(self, n, b):
        key = (int(n), float(b))
        if key not in self.cache:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                self.cache[key] = float(self.fc.FC_ints(int(n), float(b), alpha=ALPHA)[1])
            if caught:
                raise RuntimeError(f"FC warning at n={n}, b={b}: {[str(w.message) for w in caught]}")
        return self.cache[key]


def load_background(background_path, response_path, edges):
    src = {}
    with h5py.File(background_path, "r") as f:
        for name in TOTAL_SOURCES:
            g = f[name]
            if not (np.array_equal(g["hist_s1c_edges"][:], edges[0]) and np.array_equal(g["hist_y_edges"][:], edges[1])):
                raise ValueError(f"{name}: histogram geometry differs from signal")
            src[name] = {"h": g["hist_counts"][:], "m": g["hist_mc_counts"][:].astype(float),
                         "v": g["hist_mc_weight_variance"][:]}
    centers = [0.5 * (e[1:] + e[:-1]) for e in edges]
    s1, y = np.meshgrid(*centers, indexing="ij")
    mu, sig = load_band(response_path)(s1)
    return src, s1, (y - mu) / sig


def load_signal(f, mass_key):
    g = f[f"ws2024/{mass_key}"]
    if g.attrs["reference_sigma_cm2"] != REFERENCE_SIGMA_CM2:
        raise ValueError(f"{mass_key}: unexpected reference sigma")
    return {"s": g["hist_counts"][:], "s_mc": g["hist_mc_counts"][:], "s_var": g["hist_mc_weight_variance"][:],
            "mass": float(g.attrs["mass_gev"]), "edges": (g["hist_s1c_edges"][:], g["hist_y_edges"][:])}


def region_summary(mask, d, fc, b_total_hist):
    """Exact background-only limit quantiles and MC diagnostics for one region."""
    b = float(b_total_hist[mask].sum())
    s_ref = float(d["s"][mask].sum())
    to_sigma = REFERENCE_SIGMA_CM2 / s_ref
    counts = {k: poisson_quantile(q, b) for k, q in QUANTILES.items()}
    n_hi = counts["plus2sigma"] + 2
    ups = [fc.upper(n, b) for n in range(n_hi + 1)]
    if any(y < x for x, y in zip(ups, ups[1:])):
        raise AssertionError(f"FC upper limit not monotone in n at b={b}")
    mu = {k: ups[n] for k, n in counts.items()}
    floor_mu = mu["median"] * float(np.exp(norm.ppf(ALPHA)))
    var_b = float(sum(r["v"][mask].sum() for r in d["src"].values()))
    per_source = {n: {"b": float(r["h"][mask].sum()), "mc_events": int(r["m"][mask].sum())}
                  for n, r in d["src"].items()}
    return {
        "n_bins": int(mask.sum()),
        "background_events": b,
        "signal_events_at_reference_sigma": s_ref,
        "signal_fraction_of_full_roi": s_ref / float(d["s"].sum()),
        "quantile_counts": counts,
        "fc_upper_signal_events": mu,
        "sigma_ul_cm2": {k: v * to_sigma for k, v in mu.items()},
        "sigma_ul_power_constrained_cm2": {k: max(v, floor_mu) * to_sigma for k, v in mu.items()},
        "power_floor_cm2": floor_mu * to_sigma,
        "background_mc_events": int(sum(r["m"][mask].sum() for r in d["src"].values())),
        "background_mc_stat_sigma_events": math.sqrt(var_b),
        "background_effective_sample_size": (b * b / var_b) if var_b > 0 else None,
        "signal_mc_events": int(d["s_mc"][mask].sum()),
        "signal_mc_stat_sigma_events": float(math.sqrt(d["s_var"][mask].sum())),
        "per_source": per_source,
    }


def scan(d, b_hist, fc):
    s = d["s"].ravel()
    b = b_hist.ravel()
    cand = np.flatnonzero(s > 0)
    with np.errstate(divide="ignore"):
        ratio = s[cand] / b[cand]
    order = cand[np.lexsort((-s[cand], -ratio))]  # primary: ratio desc; tie: s desc
    cum_b, cum_s = np.cumsum(b[order]), np.cumsum(s[order])
    steps = []
    for k in range(len(order)):
        if cum_b[k] > B_CAP:
            break
        n_med = poisson_quantile(0.5, cum_b[k])
        mu_med = fc.upper(n_med, cum_b[k])
        steps.append({"k": k + 1, "B": float(cum_b[k]), "S": float(cum_s[k]), "n_med": n_med,
                      "mu_med": mu_med, "sigma_med": REFERENCE_SIGMA_CM2 * mu_med / float(cum_s[k])})
    sig = np.array([st["sigma_med"] for st in steps])
    best = int(np.argmin(sig))
    near = np.flatnonzero(sig <= NEAR_OPTIMUM * sig[best])
    mask = np.zeros(s.size, bool)
    mask[order[: best + 1]] = True
    return steps, best, near, order, mask.reshape(d["s"].shape)


def median_jumps(b_max):
    """Background values where the Poisson median steps from n-1 to n."""
    from scipy.special import gammainccinv
    out, n = [], 1
    while True:
        b = float(gammainccinv(n, 0.5))
        if b > b_max:
            return out
        out.append(b)
        n += 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for a in ("--signal", "--background", "--response", "--fc-source", "--outdir"):
        p.add_argument(a, required=True, type=Path)
    args = p.parse_args()
    if not (args.outdir.is_absolute() and args.outdir.is_dir()):
        raise ValueError("--outdir must be an existing absolute scratch directory")
    t0 = time.perf_counter()
    src_path = args.fc_source.resolve()
    if sha256(src_path) != EXPECTED_FC_SHA256:
        raise ValueError("FC source is not the reviewed, source-pinned file")
    sys.path.insert(0, str(src_path.parent))
    import FeldmanCousins
    fc = FC(FeldmanCousins)

    with h5py.File(args.signal, "r") as sf:
        if not sf.attrs.get("completed", False):
            raise ValueError("WIMP archive is not marked complete")
        keys = sorted(sf["ws2024"].keys(), key=lambda k: float(sf[f"ws2024/{k}"].attrs["mass_gev"]))
        signals = [load_signal(sf, k) for k in keys]
    src, s1, u = load_background(args.background, args.response, signals[0]["edges"])
    b_hist = sum(r["h"] for r in src.values())

    per_mass, rows = {}, []
    for sig in signals:
        if not all(np.array_equal(a, b) for a, b in zip(sig["edges"], signals[0]["edges"])):
            raise ValueError("Signal histogram geometry differs between masses")
        d = {**sig, "src": src, "s1": s1, "u": u}
        steps, best, near, order, mask = scan(d, b_hist, fc)
        summ = region_summary(mask, d, fc, b_hist)
        if not math.isclose(summ["sigma_ul_cm2"]["median"], steps[best]["sigma_med"], rel_tol=1e-12):
            raise AssertionError("Scan step and region summary disagree")
        per_mass[f"{sig['mass']:g}"] = {
            "mass_gev": sig["mass"], "chosen_k": steps[best]["k"], "steps_scanned": len(steps),
            "near_optimum_k": [int(steps[i]["k"]) for i in near],
            "near_optimum_B_range": [steps[near.min()]["B"], steps[near.max()]["B"]],
            "region": summ, "steps": steps,
            "selected_flat_indices": np.flatnonzero(mask.ravel()).tolist(),
        }
        pc = summ["sigma_ul_power_constrained_cm2"]
        rows.append([sig["mass"], pc["median"], pc["minus2sigma"], pc["minus1sigma"], pc["plus1sigma"],
                     pc["plus2sigma"], summ["sigma_ul_cm2"]["median"], summ["background_events"],
                     summ["signal_fraction_of_full_roi"], summ["quantile_counts"]["median"], summ["n_bins"]])
        print(f"{sig['mass']:8g} GeV  k={steps[best]['k']:5d} B={summ['background_events']:.3f} "
              f"S_frac={summ['signal_fraction_of_full_roi']:.3f} n_med={summ['quantile_counts']['median']} "
              f"sigma_med={summ['sigma_ul_cm2']['median']:.4g}")

    header = ("mass_gev,median_pc_cm2,minus2sigma_pc_cm2,minus1sigma_pc_cm2,plus1sigma_pc_cm2,plus2sigma_pc_cm2,"
              "median_unconstrained_cm2,background_events,signal_fraction_of_roi,n_med,n_bins")
    csv = args.outdir / "fc_sensitivity_scan.csv"
    csv.write_text(header + "\n" + "\n".join(",".join(f"{x:.10g}" for x in r) for r in rows) + "\n")

    report = {
        "scope": "FC counting-region scan over the 26-mass grid, WS2024 3.3 t-y; scratch output, not a declared result.",
        "alpha": ALPHA, "b_cap_events": B_CAP, "near_optimum_factor": NEAR_OPTIMUM,
        "ranking": "raw s_i/b_i desc (b_i = 0 first), ties by s_i desc; bins with s_i = 0 never added; no MC-sparsity guard (user decision).",
        "median_definition": "exact: FC_ints at Poisson(B) median count (monotonicity asserted)",
        "power_constraint": "max(limit, exp(norm.ppf(0.1)) x median), in signal-event units",
        "poisson_median_jumps_B": median_jumps(B_CAP),
        "masses": per_mass,
        "inputs": {k: {"path": str(p), "sha256": sha256(p)} for k, p in
                   (("signal", args.signal), ("background", args.background), ("response", args.response))},
        "fc_sha256": EXPECTED_FC_SHA256,
        "fc_calls": len(fc.cache),
        "elapsed_seconds": time.perf_counter() - t0,
    }
    out = args.outdir / "scan.json"
    out.write_text(json.dumps(report, indent=1, allow_nan=False) + "\n")
    print(f"Wrote {out} ({len(fc.cache)} FC calls, {report['elapsed_seconds']:.1f} s)")


if __name__ == "__main__":
    main()
