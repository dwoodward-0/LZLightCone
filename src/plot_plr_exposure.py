"""Light-format figure: 2D-PLR vs FC counting median sensitivity at two exposures.

Top: power-constrained medians of the binned 2D PLR (solid, with the UNCONSTRAINED
+-1 sigma toy band) and the exact FC counting median (dashed) at each exposure.
Bottom: PLR / FC median ratio per exposure, with the PLR bootstrap median
uncertainty. No LZ reference (the extended exposure has none). Scratch output only.

Each --run is LABEL,PLR_CSV,FC_CSV; the per-mass PLR JSON sit beside PLR_CSV.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COLORS = ["#1b1b1b", "#2a5db0", "#c0504d"]


def read_csv(path):
    lines = [ln for ln in Path(path).read_text().splitlines() if ln and not ln.startswith("#")]
    return np.genfromtxt(lines, delimiter=",", names=True)


def load_run(spec):
    label, plr_csv, fc_csv = spec.split(",")
    plr, fc = read_csv(plr_csv), read_csv(fc_csv)
    m = plr["mass_gev"]
    if not np.array_equal(m, fc["mass_gev"]):
        raise ValueError(f"{label}: PLR and FC mass grids differ")
    lo, hi = [], []
    for mm in m:
        j = json.loads((Path(plr_csv).parent / f"mass_{mm:g}gev.json").read_text())
        if j["mass_gev"] != mm:
            raise ValueError(f"{label}: JSON mass mismatch at {mm}")
        lo.append(j["limit_quantiles_cm2"]["minus1sigma"])
        hi.append(j["limit_quantiles_cm2"]["plus1sigma"])
    return label, m, plr, fc, np.array(lo), np.array(hi)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", action="append", required=True, help="LABEL,PLR_CSV,FC_CSV (repeatable)")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    runs = [load_run(s) for s in args.run]

    fig, (ax, rx) = plt.subplots(2, 1, figsize=(7.2, 6.8), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1.4], "hspace": 0.06})
    for (label, m, plr, fc, lo, hi), col in zip(runs, COLORS):
        ax.fill_between(m, lo, hi, color=col, alpha=0.12, lw=0)
        ax.plot(m, plr["median_pc_cm2"], color=col, lw=2, label=f"{label}: 2D PLR median (PC)")
        ax.plot(m, fc["median_pc_cm2"], color=col, lw=1.4, ls="--", label=f"{label}: FC counting median")
        ratio = plr["median_pc_cm2"] / fc["median_pc_cm2"]
        err = plr["median_bootstrap_sigma_cm2"] / fc["median_pc_cm2"]
        rx.errorbar(m, ratio, yerr=err, fmt="o-", ms=3.5, lw=1, color=col, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"SI WIMP-nucleon $\sigma$ [cm$^2$], 90% CL")
    ax.grid(which="major", color="#e4e4e4", lw=0.6)
    ax.legend(fontsize=8, loc="upper left", frameon=False)
    ax.set_title(r"Binned 2D PLR vs FC counting: median sensitivity vs exposure"
                 "\n(shaded: PLR $\\pm1\\sigma$ of unconstrained toy limits)", fontsize=10.5)
    rx.axhline(1, color="#8a8a8a", lw=0.8)
    rx.set_ylabel("PLR / FC median")
    rx.set_xlabel(r"WIMP mass [GeV/c$^2$]")
    rx.grid(which="major", color="#e4e4e4", lw=0.6)
    rx.legend(fontsize=8, loc="lower right", frameon=False, ncol=len(runs))
    fig.savefig(args.output, dpi=160, bbox_inches="tight")
    print(f"wrote {args.output}")
    print("mass  " + "  ".join(f"{r[0]:>18s}" for r in runs))
    for i, mm in enumerate(runs[0][1]):
        print(f"{mm:6g} " + "  ".join(
            f"{r[2]['median_pc_cm2'][i] / r[3]['median_pc_cm2'][i]:8.3f} +- "
            f"{r[2]['median_bootstrap_sigma_cm2'][i] / r[3]['median_pc_cm2'][i]:.3f}" for r in runs))


if __name__ == "__main__":
    main()
