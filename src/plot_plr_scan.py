"""Light-format figure: 2D-PLR median sensitivity vs LZ WS2024-only (slide 24).

Top: our power-constrained median with +-1/2 sigma bands of the UNCONSTRAINED toy
limits (LZ's convention; the PC'd -2 sigma edge would just repeat -1 sigma), LZ's
WS2024-only median and +-1 sigma edges, and (optionally) the FC counting median.
Bottom: ratio of medians to LZ's, with our bootstrap median uncertainty.

Acceptance (sensitivity_vs_lz_validation, like-for-like WS2024-only target):
  A1  at 40 GeV/c2 our PC median is within 25% of LZ's WS2024-only median,
      |ours/LZ - 1| <= 0.25 (user decision 2026-11-19; replaces the original
      |ours - LZ| <= 2 x bootstrap sigma gate, which ignored the uncertainty of
      LZ's curve, its digitization and our binned single-sample likelihood).
      The bootstrap pull is still reported, informationally;
  A2  at every mass where LZ's +-1 sigma band is published (finite), our PC median
      lies inside it.
Both are evaluated, printed and written onto the figure; the script does not fail
on a FAIL verdict (the figure is the record).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stat_decisions import add_decision_args

A1_MASS_GEV = 40.0
A1_MAX_REL_DEV = 0.25
BANDS = ("minus2sigma", "minus1sigma", "plus1sigma", "plus2sigma")


def read_csv(path):
    lines = [ln for ln in Path(path).read_text().splitlines() if ln and not ln.startswith("#")]
    return np.genfromtxt(lines, delimiter=",", names=True)


def unconstrained_bands(ours, scan_csv, scan_dir):
    if all(f"{q}_cm2" in ours.dtype.names for q in BANDS):
        return {q: ours[f"{q}_cm2"] for q in BANDS}
    band = {q: [] for q in BANDS}  # older CSVs: per-mass JSON beside the scan
    for mm in ours["mass_gev"]:
        j = json.loads(((scan_dir or scan_csv.parent) / f"mass_{mm:g}gev.json").read_text())
        if j["mass_gev"] != mm:
            raise ValueError(f"JSON mass mismatch at {mm}")
        for q in band:
            band[q].append(j["limit_quantiles_cm2"][q])
    return {q: np.array(v) for q, v in band.items()}


def acceptance(m, ours, ref):
    i = int(np.flatnonzero(m == A1_MASS_GEV)[0])
    med, boot, lz = ours["median_pc_cm2"][i], ours["median_bootstrap_sigma_cm2"][i], ref["median"][i]
    pull = (med - lz) / boot
    a1 = {"mass_gev": A1_MASS_GEV, "ours_cm2": float(med), "bootstrap_sigma_cm2": float(boot),
          "lz_ws2024_only_cm2": float(lz), "ratio": float(med / lz), "pull": float(pull),
          "max_rel_dev": A1_MAX_REL_DEV, "pass": bool(abs(med / lz - 1.0) <= A1_MAX_REL_DEV)}
    have = np.isfinite(ref["minus1sigma"]) & np.isfinite(ref["plus1sigma"])
    inside = (ours["median_pc_cm2"] >= ref["minus1sigma"]) & (ours["median_pc_cm2"] <= ref["plus1sigma"])
    a2 = {"masses_with_band": int(have.sum()), "inside": int((inside & have).sum()),
          "outside_masses_gev": [float(x) for x in m[have & ~inside]],
          "ratio_range": [float(np.min((ours["median_pc_cm2"] / ref["median"])[have])),
                          float(np.max((ours["median_pc_cm2"] / ref["median"])[have]))],
          "pass": bool((inside | ~have).all())}
    return a1, a2


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scan", required=True, type=Path, help="PLR sensitivity CSV")
    ap.add_argument("--reference", required=True, type=Path, help="lz_ws2024_only_sensitivity_slide24.csv")
    ap.add_argument("--fc", type=Path, help="optional FC sensitivity CSV")
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--scan-dir", type=Path, help="per-mass JSON dir for older CSVs (default: --scan's directory)")
    add_decision_args(ap, ["statistics", "median_sensitivity_definition", "power_constraint", "test_statistic_sidedness"], "both")
    args = ap.parse_args()
    ours, ref = read_csv(args.scan), read_csv(args.reference)
    m = ours["mass_gev"]
    if not np.array_equal(m, ref["mass_gev"]):
        raise ValueError("mass grids differ")
    band = unconstrained_bands(ours, args.scan, args.scan_dir)
    a1, a2 = acceptance(m, ours, ref)

    fig, (ax, rx) = plt.subplots(2, 1, figsize=(7.2, 6.8), sharex=True,
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
    verdict = (f"A1 (40 GeV, within {A1_MAX_REL_DEV:.0%} of LZ): {'PASS' if a1['pass'] else 'FAIL'}, "
               f"ratio {a1['ratio']:.3f}, pull {a1['pull']:+.1f}\n"
               f"A2 (median inside LZ $\\pm1\\sigma$): {'PASS' if a2['pass'] else 'FAIL'}, "
               f"{a2['inside']}/{a2['masses_with_band']} masses")
    ax.text(0.98, 0.04, verdict, transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color="#3a3a3a")

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
    fig.savefig(args.output, dpi=160, bbox_inches="tight", format="png")
    print(f"wrote {args.output}")
    print("acceptance:", json.dumps({"A1": a1, "A2": a2}, indent=1))
    for mm, r, e in zip(m, ratio, err):
        print(f"{mm:8g}  ratio {r:.3f} +- {e:.3f}")


if __name__ == "__main__":
    main()
