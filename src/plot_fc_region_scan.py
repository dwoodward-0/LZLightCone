"""Validation figures for src/fc_region_scan.py (light format, PNG).

fig1_scan_<m>gev.png        sigma_med and signal fraction vs cumulative B_k
fig2_region_<m>gev.png      selected bins in (S1c, log10 S2c) + band overlays
fig3_mc_support_<m>gev.png  background MC support maps + per-source B in region
fig4_sensitivity.png        FC median curve and bands vs LZ WS2024 median/bands
fig5_summary.png            per-mass summary table
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from validate_backgrounds import TOTAL_SOURCES

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
OURS, LZ = "#2a78d6", "#52514e"
plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                     "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.edgecolor": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "font.size": 9, "axes.titlesize": 10, "legend.frameon": False})


def load(args):
    rep = json.loads((args.outdir / "scan.json").read_text())
    key = f"{args.mass:g}"
    with h5py.File(args.signal) as f:
        g = f[f"ws2024/mass_{key}gev"]
        s = g["hist_counts"][:]
        edges = (g["hist_s1c_edges"][:], g["hist_y_edges"][:])
    with h5py.File(args.background) as f:
        h = {n: f[n]["hist_counts"][:] for n in TOTAL_SOURCES}
        m = {n: f[n]["hist_mc_counts"][:] for n in TOTAL_SOURCES}
    bands = {}
    with h5py.File(args.response) as f:
        for k in ("er_ch3t", "nr"):
            bands[k] = (f[f"{k}/band_s1c_centers"][:], f[f"{k}/band_log10s2c_quantiles"][:])
    mk = np.zeros(s.size, bool)
    mk[rep["masses"][key]["selected_flat_indices"]] = True
    return rep, rep["masses"][key], s, h, m, edges, bands, mk.reshape(s.shape)


def overlay_bands(ax, bands, legend=False):
    styles = {"er_ch3t": ("#52514e", "ER band"), "nr": ("#e34948", "NR band")}
    for k, (x, q) in bands.items():
        col, lab = styles[k]
        ax.plot(x, q[2], color=col, lw=1.2, label=f"{lab} median" if legend else None)
        for i in (1, 3):
            ax.plot(x, q[i], color=col, lw=0.8, ls="--", label=f"{lab} ±1σ" if legend and i == 1 else None)
        for i in (0, 4):
            ax.plot(x, q[i], color=col, lw=0.6, ls=":", label=f"{lab} ±2σ" if legend and i == 0 else None)


def fig_scan(rep, v, mass, out):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, height_ratios=[3, 2])
    xmax = 12.0
    for a in (a1, a2):
        for j in rep["poisson_median_jumps_B"]:
            if j < xmax:
                a.axvline(j, color=INK2, lw=0.5, alpha=0.5)
    st = v["steps"]
    B = np.array([x["B"] for x in st])
    sig = np.array([x["sigma_med"] for x in st])
    s_full = v["region"]["signal_events_at_reference_sigma"] / v["region"]["signal_fraction_of_full_roi"]
    frac = np.array([x["S"] for x in st]) / s_full
    k = v["chosen_k"] - 1
    a1.plot(B, sig, color=OURS, lw=1.6)
    a2.plot(B, frac, color=OURS, lw=1.6)
    for a, yv in ((a1, sig[k]), (a2, frac[k])):
        a.plot(B[k], yv, "o", ms=8, color=OURS, mec=SURFACE, mew=2, zorder=5)
    a1.axhline(rep["near_optimum_factor"] * sig[k], color=OURS, lw=0.7, ls=":")
    a1.annotate(f"chosen: B = {B[k]:.3f}, n_med = {st[k]['n_med']}, σ_med = {sig[k]:.3g} cm²",
                (B[k], sig[k]), xytext=(15, -25), textcoords="offset points", fontsize=8, color=INK)
    a1.set_yscale("log")
    a1.set_ylabel("median FC 90% UL  σ_med (cm²)")
    a1.set_title(f"{mass:g} GeV/c² progressive-voxel scan (ranking s/b)\nvertical lines: Poisson-median jumps "
                 "(B = 0.69, 1.68, 2.67, …);  dotted: 1.05 × minimum", loc="left")
    a2.set_ylabel("signal fraction of full ROI")
    a2.set_xlabel("cumulative background B_k (events, 3.3 t·y)")
    a2.set_xlim(0, xmax)
    a2.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(out / f"fig1_scan_{mass:g}gev.png", dpi=150)
    plt.close(fig)


def fig_region(v, mass, s, edges, bands, mask, out):
    e1, e2 = edges
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.pcolormesh(e1, e2, np.where(mask, 1.0, np.nan).T, cmap=ListedColormap([OURS]), alpha=0.7, shading="flat")
    ax.plot([], [], "s", color=OURS, alpha=0.7, label="selected bins")
    centers = [0.5 * (e[1:] + e[:-1]) for e in edges]
    lev = np.quantile(s[s > 0], [0.5, 0.9, 0.99])
    ax.contour(centers[0], centers[1], s.T, levels=lev, colors=INK, linewidths=0.6, alpha=0.5)
    ax.plot([], [], color=INK, lw=0.6, alpha=0.5, label=f"{mass:g} GeV signal density (50/90/99% bin quantiles)")
    overlay_bands(ax, bands, legend=True)
    r = v["region"]
    ax.set_title(f"{mass:g} GeV/c² FC counting region: {r['n_bins']} bins, B = {r['background_events']:.3f} events, "
                 f"signal fraction = {r['signal_fraction_of_full_roi']:.3f}", loc="left")
    ax.set_xlim(3, 80)
    ax.set_ylim(e2[0], 4.0)
    ax.grid(False)
    ax.set_xlabel("S1c (phd)")
    ax.set_ylabel("log10(S2c / phd)")
    ax.legend(loc="upper left", fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out / f"fig2_region_{mass:g}gev.png", dpi=150)
    plt.close(fig)


def fig_mc_support(v, mass, s, h, m, edges, mask, out):
    e1, e2 = edges
    btot, mtot = sum(h.values()), sum(m.values())
    zero_b = (s > 0) & (btot == 0)
    centers = [0.5 * (e[1:] + e[:-1]) for e in edges]
    X, Y = np.meshgrid(*centers, indexing="ij")
    fig = plt.figure(figsize=(13, 9.5))
    gs = fig.add_gridspec(2, 3)
    w_pb = h["214Pb_betas"].sum() / m["214Pb_betas"].sum()
    panels = [(fig.add_subplot(gs[0, 0]), mtot, "Greys", "total background MC events / bin (all sources)"),
              (fig.add_subplot(gs[0, 1]), m["214Pb_betas"], "Greys", f"214Pb MC events / bin ({w_pb:.3g} events each)"),
              (fig.add_subplot(gs[0, 2]), s, "Blues", f"{mass:g} GeV signal (events/bin at 1e-36 cm²)")]
    for ax, arr, cmap, title in panels:
        pc = ax.pcolormesh(e1, e2, np.where(arr > 0, arr, np.nan).T, cmap=cmap, norm=LogNorm(), shading="flat")
        fig.colorbar(pc, ax=ax, fraction=0.046)
        ax.contour(centers[0], centers[1], mask.T.astype(float), levels=[0.5], colors=OURS, linewidths=1.2)
        ax.plot(X[zero_b], Y[zero_b], "x", color=INK, ms=4, mew=0.8)
        ax.set_title(title, loc="left", fontsize=9)
        ax.set_xlim(3, 80)
        ax.set_ylim(e2[0], 4.0)
        ax.set_xlabel("S1c (phd)")
        ax.grid(False)
    a0 = panels[0][0]
    a0.set_ylabel("log10(S2c / phd)")
    a0.plot([], [], "x", color=INK, label="signal bin with zero background MC")
    a0.plot([], [], color=OURS, label="selected region outline")
    a0.legend(loc="upper right", fontsize=7)

    ax = fig.add_subplot(gs[1, :])
    ps = v["region"]["per_source"]
    srcs = list(TOTAL_SOURCES)
    vals = np.array([ps[n]["b"] for n in srcs])
    y = np.arange(len(srcs))
    ax.barh(y, np.maximum(vals, 1e-5), height=0.6, color=OURS)
    for j, n in enumerate(srcs):
        ax.text(max(vals[j], 1e-5) * 1.15, y[j], f"{ps[n]['b']:.3g} events from {ps[n]['mc_events']} MC events",
                va="center", fontsize=7, color=INK2)
    ax.set_yticks(y, srcs)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(1e-5, 10)
    ax.set_xlabel("background events in selected region (log scale; bars at 1e-5 = none)")
    ax.set_title(f"per-source background in the {mass:g} GeV/c² region, with MC support "
                 f"(total B = {v['region']['background_events']:.3f})", loc="left")
    fig.tight_layout()
    fig.savefig(out / f"fig3_mc_support_{mass:g}gev.png", dpi=150)
    plt.close(fig)


def fig_sensitivity(rep, lz_path, out):
    ms = sorted(rep["masses"].values(), key=lambda v: v["mass_gev"])
    mass = np.array([v["mass_gev"] for v in ms])
    pc = {k: np.array([v["region"]["sigma_ul_power_constrained_cm2"][k] for v in ms])
          for k in ("minus2sigma", "minus1sigma", "median", "plus1sigma", "plus2sigma")}
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 8), sharex=True, height_ratios=[3, 1.3])
    # columns: mass, limit, limit_unconstr, -2sigma, -1sigma, median, 1sigma, 2sigma, disco
    raw = np.genfromtxt(lz_path, delimiter=",", skip_header=1)
    lzm, lz_m2, lz_m1, lz_med, lz_p1, lz_p2 = raw[:, 0], raw[:, 3], raw[:, 4], raw[:, 5], raw[:, 6], raw[:, 7]
    a1.fill_between(lzm, lz_m2, lz_p2, color=LZ, alpha=0.10, lw=0, label="LZ WS2024 ±2σ")
    a1.fill_between(lzm, lz_m1, lz_p1, color=LZ, alpha=0.20, lw=0, label="LZ WS2024 ±1σ")
    a1.plot(lzm, lz_med, color=LZ, lw=1.6, ls="--", label="LZ WS2024 median (4.2 t·y, unbinned PLR)")
    a1.fill_between(mass, pc["minus2sigma"], pc["plus2sigma"], color=OURS, alpha=0.12, lw=0, label="this work ±2σ")
    a1.fill_between(mass, pc["minus1sigma"], pc["plus1sigma"], color=OURS, alpha=0.25, lw=0, label="this work ±1σ")
    a1.plot(mass, pc["median"], color=OURS, lw=2, marker="o", ms=3, label="this work: FC counting median (3.3 t·y)")
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_ylabel("SI WIMP-nucleon σ, 90% CL (cm²)")
    i40 = int(np.argmin(np.abs(mass - 40)))
    a1.set_title(f"FC cut-and-count median sensitivity vs LZ WS2024 (power-constrained)\n"
                 f"40 GeV/c²: this work {pc['median'][i40]:.3g} cm², LZ {np.interp(40, lzm, lz_med):.3g} cm²",
                 loc="left")
    a1.legend(loc="upper center", fontsize=7, ncol=2)
    if not np.allclose(lzm, mass):
        raise ValueError("LZ reference masses differ from the scan grid")
    ratio = pc["median"] / lz_med
    a2.axhline(1, color=LZ, lw=0.8)
    a2.fill_between(lzm, lz_m1 / lz_med, lz_p1 / lz_med, color=LZ, alpha=0.2, lw=0, label="LZ ±1σ / median")
    a2.plot(mass, ratio, color=OURS, lw=1.8, marker="o", ms=3, label="this work / LZ median")
    a2.set_yscale("log")
    a2.set_ylabel("ratio to LZ median")
    a2.set_xlabel("WIMP mass (GeV/c²)")
    a2.legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "fig4_sensitivity.png", dpi=150)
    plt.close(fig)
    return mass, pc, lz_med


def fig_summary(rep, lz_med, out):
    ms = sorted(rep["masses"].values(), key=lambda v: v["mass_gev"])
    cols = ["mass (GeV)", "bins", "B (events)", "signal frac", "n_med (−2,−1,+1,+2σ)", "μ_med (events)",
            "σ_med PC (cm²)", "−1σ / +1σ (cm²)", "÷ LZ median", "bkg MC events"]
    cell = []
    for v, lzv in zip(ms, lz_med):
        r = v["region"]
        q, pc = r["quantile_counts"], r["sigma_ul_power_constrained_cm2"]
        cell.append([f"{v['mass_gev']:g}", f"{r['n_bins']}", f"{r['background_events']:.3f}",
                     f"{r['signal_fraction_of_full_roi']:.3f}",
                     "{median} ({minus2sigma},{minus1sigma},{plus1sigma},{plus2sigma})".format(**q),
                     f"{r['fc_upper_signal_events']['median']:.3f}", f"{pc['median']:.3e}",
                     f"{pc['minus1sigma']:.2e} / {pc['plus1sigma']:.2e}", f"{pc['median'] / lzv:.2f}",
                     f"{r['background_mc_events']}"])
    fig, ax = plt.subplots(figsize=(13, 8.5))
    ax.axis("off")
    t = ax.table(cellText=cell, colLabels=cols, loc="center", cellLoc="center")
    t.auto_set_font_size(False)
    t.set_fontsize(7.5)
    t.scale(1, 1.25)
    for (i, j), c in t.get_celld().items():
        c.set_edgecolor(GRID)
        if i == 0:
            c.set_text_props(weight="bold")
    ax.set_title("FC counting regions per mass (90% CL, 3.3 t·y, exact background-only quantiles, ranking s/b)",
                 loc="left")
    fig.tight_layout()
    fig.savefig(out / "fig5_summary.png", dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for a in ("--signal", "--background", "--response", "--lz-reference", "--outdir"):
        p.add_argument(a, required=True, type=Path)
    p.add_argument("--mass", type=float, default=40.0, help="mass for the detail figures 1-3")
    args = p.parse_args()
    rep, v, s, h, m, edges, bands, mask = load(args)
    fig_scan(rep, v, args.mass, args.outdir)
    fig_region(v, args.mass, s, edges, bands, mask, args.outdir)
    fig_mc_support(v, args.mass, s, h, m, edges, mask, args.outdir)
    _, _, lz_med = fig_sensitivity(rep, args.lz_reference, args.outdir)
    fig_summary(rep, lz_med, args.outdir)
    print("wrote figures to", args.outdir)


if __name__ == "__main__":
    main()
