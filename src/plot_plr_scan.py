"""Light-format figure: 2D-PLR median sensitivity vs LZ WS2024-only (slide 24).

Top: our power-constrained median with +-1/2 sigma bands of the UNCONSTRAINED toy
limits (LZ's convention; the PC'd -2 sigma edge would just repeat -1 sigma), LZ's WS2024-only
median and band edges, and (optionally) the FC counting median. Bottom: ratio of
medians to LZ's, with our bootstrap median uncertainty. Scratch output only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    lines = [ln for ln in Path(path).read_text().splitlines() if ln and not ln.startswith("#")]
    return np.genfromtxt(lines, delimiter=",", names=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scan", required=True, type=Path, help="plr_sensitivity_scan.csv")
    ap.add_argument("--reference", required=True, type=Path, help="lz_ws2024_only_sensitivity_slide24.csv")
    ap.add_argument("--fc", type=Path, help="optional fc_sensitivity_scan.csv")
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--scan-dir", type=Path, help="per-mass JSON dir (default: --scan's directory)")
    args = ap.parse_args()
    ours, ref = read_csv(args.scan), read_csv(args.reference)
    m = ours["mass_gev"]
    if not np.array_equal(m, ref["mass_gev"]):
        raise ValueError("mass grids differ")

    scan_dir = args.scan_dir or args.scan.parent
    band = {q: [] for q in ("minus2sigma", "minus1sigma", "plus1sigma", "plus2sigma")}
    for mm in m:
        j = json.loads((scan_dir / f"mass_{mm:g}gev.json").read_text())
        if j["mass_gev"] != mm:
            raise ValueError(f"JSON mass mismatch at {mm}")
        for q in band:
            band[q].append(j["limit_quantiles_cm2"][q])
    band = {q: np.array(v) for q, v in band.items()}

    fig, (ax, rx) = plt.subplots(2, 1, figsize=(7.2, 6.6), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.06})
    ax.fill_between(m, band["minus2sigma"], band["plus2sigma"], color="#f6e27a", lw=0,
                    label=r"ours $\pm2\sigma$")
    ax.fill_between(m, band["minus1sigma"], band["plus1sigma"], color="#8fd18a", lw=0,
                    label=r"ours $\pm1\sigma$")
    ax.plot(m, ours["median_pc_cm2"], color="#1b1b1b", lw=2, label="ours: 2D PLR median, PC (3.3 t-y)")
    ax.plot(m, ref["median"], color="#2a5db0", lw=2, ls="--", label="LZ WS2024-only median (slide 24)")
    for col, lab in (("minus1sigma", r"LZ $\pm1\sigma$"), ("plus1sigma", None)):
        ax.plot(m, ref[col], color="#2a5db0", lw=1, ls=":", label=lab)
    fc = read_csv(args.fc) if args.fc else None
    if fc is not None:
        ax.plot(fc["mass_gev"], fc["median_pc_cm2"], color="#c0504d", lw=1.4, ls="-.",
                label="ours: FC counting median")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"SI WIMP-nucleon $\sigma$ [cm$^2$], 90% CL")
    ax.set_ylim(5e-49, 3e-45)
    ax.grid(which="major", color="#e4e4e4", lw=0.6)
    ax.legend(fontsize=8, loc="upper left", frameon=False)
    ax.set_title("Median sensitivity: binned 2D PLR vs LZ WS2024-only", fontsize=11)

    ratio = ours["median_pc_cm2"] / ref["median"]
    err = ours["median_bootstrap_sigma_cm2"] / ref["median"]
    rx.axhspan(0.9, 1.1, color="#efefef", lw=0)
    rx.axhline(1, color="#8a8a8a", lw=0.8)
    rx.errorbar(m, ratio, yerr=err, fmt="o", ms=3.5, color="#1b1b1b", lw=1, label="2D PLR / LZ")
    if fc is not None:
        rx.plot(fc["mass_gev"], fc["median_pc_cm2"] / ref["median"], "s", ms=3, mfc="none",
                color="#c0504d", label="FC counting / LZ")
    rx.set_ylabel("median / LZ")
    rx.set_xlabel(r"WIMP mass [GeV/c$^2$]")
    rx.set_ylim(0.4, 1.8)
    rx.grid(which="major", color="#e4e4e4", lw=0.6)
    rx.legend(fontsize=7.5, loc="upper right", frameon=False, ncol=2)
    fig.savefig(args.output, dpi=160, bbox_inches="tight")
    print(f"wrote {args.output}")
    for mm, r, e in zip(m, ratio, err):
        print(f"{mm:8g}  ratio {r:.3f} +- {e:.3f}")


if __name__ == "__main__":
    main()
