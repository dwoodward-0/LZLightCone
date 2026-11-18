"""Light-format figure: FC counting vs binned 2D-PLR median sensitivity (fc_plr_crosscheck).

Top: PLR power-constrained median with its UNCONSTRAINED +-1 sigma toy band, and the
exact FC counting median with its +-1 sigma edges. Bottom: PLR/FC median ratio with
the acceptance envelope.

Acceptance: at every mass, |ln(PLR/FC)| <= sqrt(w_PLR^2 + w_FC^2), with w the upper
(+1 sigma) log half-width of each method's band, ln(plus1sigma / median). Upper edges
are used because the power constraint and the n_med = 0 counting regime pin the lower
edges to the median. Evaluated, printed and written on the figure; the script does
not fail on a FAIL verdict. The FC regions' retained signal fraction is reported so
the "cut-and-count retains efficiency" qualifier can be read off.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stat_decisions import add_decision_args


def read_csv(path):
    lines = [ln for ln in Path(path).read_text().splitlines() if ln and not ln.startswith("#")]
    return np.genfromtxt(lines, delimiter=",", names=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--plr", required=True, type=Path, help="PLR sensitivity CSV")
    ap.add_argument("--fc", required=True, type=Path, help="FC sensitivity CSV")
    ap.add_argument("--output", required=True, type=Path)
    add_decision_args(ap, ["statistics", "median_sensitivity_definition", "power_constraint", "test_statistic_sidedness"], "both")
    args = ap.parse_args()
    plr, fc = read_csv(args.plr), read_csv(args.fc)
    m = plr["mass_gev"]
    if not np.array_equal(m, fc["mass_gev"]):
        raise ValueError("PLR and FC mass grids differ")
    if not np.allclose(plr["exposure_scale"], 1.0):
        raise ValueError("cross-check is defined at the WS2024 exposure")

    p_med, f_med = plr["median_pc_cm2"], fc["median_pc_cm2"]
    w_p = np.log(plr["plus1sigma_cm2"] / plr["median_unconstrained_cm2"])
    w_f = np.log(fc["plus1sigma_pc_cm2"] / f_med)
    env = np.sqrt(w_p ** 2 + w_f ** 2)
    lr = np.log(p_med / f_med)
    ok = np.abs(lr) <= env
    acc = {"masses": int(m.size), "inside": int(ok.sum()), "outside_masses_gev": [float(x) for x in m[~ok]],
           "ratio_range": [float(np.exp(lr).min()), float(np.exp(lr).max())],
           "max_abs_ln_ratio_over_envelope": float(np.max(np.abs(lr) / env)),
           "fc_signal_fraction_range": [float(fc["signal_fraction_of_roi"].min()),
                                        float(fc["signal_fraction_of_roi"].max())],
           "pass": bool(ok.all())}

    fig, (ax, rx) = plt.subplots(2, 1, figsize=(7.2, 6.8), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1.4], "hspace": 0.06})
    ax.fill_between(m, plr["minus1sigma_cm2"], plr["plus1sigma_cm2"], color="#1b1b1b", alpha=0.12, lw=0,
                    label=r"2D PLR $\pm1\sigma$ (unconstrained toys)")
    ax.plot(m, p_med, color="#1b1b1b", lw=2, label="2D PLR median (PC)")
    ax.plot(m, f_med, color="#c0504d", lw=1.6, ls="--", label="FC counting median (PC)")
    ax.plot(m, fc["plus1sigma_pc_cm2"], color="#c0504d", lw=0.9, ls=":", label=r"FC counting $+1\sigma$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"SI WIMP-nucleon $\sigma$ [cm$^2$], 90% CL")
    ax.grid(which="major", color="#e4e4e4", lw=0.6)
    ax.legend(fontsize=8, loc="upper left", frameon=False)
    ax.set_title("Statistical cross-check: FC counting vs binned 2D PLR (WS2024, 3.3 t-y)", fontsize=10.5)
    ax.text(0.98, 0.04, f"{'PASS' if acc['pass'] else 'FAIL'}: {acc['inside']}/{acc['masses']} masses inside "
            f"the combined +1$\\sigma$ envelope;\nPLR/FC {acc['ratio_range'][0]:.2f}-{acc['ratio_range'][1]:.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color="#3a3a3a")

    rx.fill_between(m, np.exp(-env), np.exp(env), color="#efefef", lw=0, label="acceptance envelope")
    rx.axhline(1, color="#8a8a8a", lw=0.8)
    rx.errorbar(m, p_med / f_med, yerr=plr["median_bootstrap_sigma_cm2"] / f_med, fmt="o", ms=3.5, lw=1,
                color="#1b1b1b", label="PLR / FC median")
    rx.set_ylim(0.3, 3.0)
    rx.set_yscale("log")
    rx.set_yticks([0.5, 1, 2])
    rx.set_yticklabels(["0.5", "1", "2"])
    rx.yaxis.set_minor_formatter(NullFormatter())
    rx.set_ylabel("PLR / FC")
    rx.set_xlabel(r"WIMP mass [GeV/c$^2$]")
    rx.grid(which="major", color="#e4e4e4", lw=0.6)
    rx.legend(fontsize=7.5, loc="upper right", frameon=False, ncol=2)
    fig.savefig(args.output, dpi=160, bbox_inches="tight", format="png")
    print(f"wrote {args.output}")
    print("acceptance:", json.dumps(acc, indent=1))


if __name__ == "__main__":
    main()
