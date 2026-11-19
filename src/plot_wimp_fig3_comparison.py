"""40 GeV/c2 WIMP region comparison vs WS2024 Fig. 3 (All Data, purple).

LZ's Fig. 3 top panel shades the regions containing 68% and 95% of a
40 GeV/c2 WIMP in (S1c, log10 S2c). We extract both purple polygons exactly
from the pinned paper's vector fill paths (tick-label axis calibration from
src/digitize_accidentals_ws2024.py) and compare them with our declared
40 GeV/c2 signal model (wimp_signal_models, ws2024/mass_40gev):

  - containment: the fraction of our ROI-selected MC events falling inside
    each LZ polygon (unbinned, on the archive's uniformly thinned display
    sample: 50,000 of the 59,852 full-ROI events, drawn without replacement,
    equal weights, so an unbiased estimate of the full-sample fraction);
  - figure: left, our full 59,852-event MC histogram (archive grid, 1 phd x
    0.02 dex) under LZ's regions; right, our own 68%/95% highest-density
    contours of a display-only Gaussian-smoothed copy of that histogram
    (1.1 phd, 0.040 dex) against LZ's.

No refit, renormalization, new threshold or new NEST production. The ROI is
the model's corrected-space ROI (3 < S1c < 80 phd, S2c >= 493 phd,
log10 S2c < 4.5).

Acceptance gate (user decision 2026-11-19): the containment inside LZ's 68%
region must lie within 0.68 +/- 0.05, and inside the 95% region within
0.95 +/- 0.03 (absolute fractions). It is a tolerance, not a significance
test: published regions are not a 2D PDF, and LZ's contour/MC/construction
uncertainty is unavailable, so a test against our MC sigma alone would be
meaningless. The pull against our MC sigma is still reported.

One file per invocation, chosen by the --output suffix: .png draws the
figure (wimp_fig3_validation), .json writes the metrics
(wimp_fig3_validation_metrics), including the extracted polygon vertices.

Rebuilt 2026-11-19 from the 2026-11-13 scratch probe
/tmp/lz_wimp_fig3_comparison.py, whose event-level /tmp archive no longer
exists; it read the same 59,852-event realization event by event (68%:
66.00%, 95%: 94.37% inside).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.path import Path as MplPath  # noqa: E402
import numpy as np  # noqa: E402
import pymupdf  # noqa: E402

from digitize_accidentals_ws2024 import SHA256, calibrate  # noqa: E402

MASS_GROUP = "ws2024/mass_40gev"
N_FULL_ROI = 59852
PURPLE = (0.517647, 0.0, 0.803922)
SURFACE, FG, MUTED, GRID = "#fcfcfb", "#20201e", "#62625c", "#e4e3e0"
OUR, LZ = "#D55E00", "#7B61A8"
SMOOTH_PHD, SMOOTH_DEX = 1.1, 0.040
# Acceptance gate: region -> (nominal fraction, absolute tolerance).
GATE = {"q68": (0.68, 0.05), "q95": (0.95, 0.03)}
CAPTION = ("Dark and light gray and purple shading show the regions containing "
           "68% and 95% of the ER portion of the background model and a "
           "40 GeV/c2 WIMP, respectively.")


def regions_from_pdf(pdf):
    """Exact purple fill polygons of Fig. 3's top ('All Data') panel."""
    digest = hashlib.sha256(Path(pdf).read_bytes()).hexdigest()
    assert digest == SHA256, "pinned WS2024 PDF has changed"
    page = pymupdf.open(pdf)[6]  # printed page 7
    flat = " ".join(page.get_text().replace("-\n", "").split())
    assert "a 40 GeV/c2 WIMP" in flat, "Fig. 3 caption is missing"
    xf, yf = calibrate(page)
    regions = []
    for drawing in page.get_drawings():
        fill = drawing["fill"]
        if fill is None or not all(abs(a - b) < .002
                                   for a, b in zip(fill, PURPLE)):
            continue
        rect = drawing["rect"]
        if not (45 < rect.y0 < 206 and rect.y1 < 206):
            continue  # All Data panel, not Radon Tagged
        assert drawing["type"] == "f"
        assert all(item[0] == "l" for item in drawing["items"])
        pts = []
        for item in drawing["items"]:
            pts.extend([(item[1].x, item[1].y), (item[2].x, item[2].y)])
        pts = np.asarray(pts)
        pts = np.column_stack([np.polyval(xf, pts[:, 0]),
                               np.polyval(yf, pts[:, 1])])
        keep = np.ones(len(pts), dtype=bool)
        keep[1:] = np.any(np.diff(pts, axis=0) != 0, axis=1)
        pts = pts[keep]
        assert np.array_equal(pts[0], pts[-1]), "region not closed"
        regions.append(pts[:-1])
    assert len(regions) == 2
    regions.sort(key=lambda p: np.ptp(p[:, 0]))
    regions = {"q68": regions[0], "q95": regions[1]}
    assert MplPath(regions["q95"], closed=True).contains_points(
        regions["q68"]).all()
    return regions


def load_model(path):
    with h5py.File(path) as f:
        g = f[MASS_GROUP]
        hist = g["hist_mc_counts"][()]
        e1, ey = g["hist_s1c_edges"][()], g["hist_y_edges"][()]
        s1, y, w = g["s1c"][()], g["log10_s2c"][()], g["weight"][()]
    assert int(hist.sum()) == N_FULL_ROI, hist.sum()
    assert np.ptp(w) <= 1e-9 * w[0], "display weights not uniform"
    roi = (s1 > 3) & (s1 < 80) & (y >= np.log10(493.0) - 1e-12) & (y < 4.5)
    assert roi.all(), "display sample outside the ROI"
    # The display sample is a subset of the full sample, bin by bin.
    hd = np.histogram2d(s1, y, bins=(e1, ey))[0]
    assert np.all(hd <= hist), "display sample is not a subset"
    return s1, y, hist, (e1, ey)


def containment(s1, y, regions):
    n = len(s1)
    sig = lambda f: float(np.sqrt(f * (1 - f) / n))
    pts = np.column_stack([s1, y])
    out = {}
    for name, poly in regions.items():
        inside = MplPath(poly, closed=True).contains_points(pts)
        f = float(inside.mean())
        out[name] = {
            "fraction_inside": f, "n_inside": int(inside.sum()),
            "mc_sigma_fraction": sig(f),
            "region_bounds_s1c": [float(poly[:, 0].min()),
                                  float(poly[:, 0].max())],
            "region_bounds_log10_s2c": [float(poly[:, 1].min()),
                                        float(poly[:, 1].max())]}
    for key, sel in (("S1c_gt_40", s1 > 40), ("log10_s2c_lt_2.75", y < 2.75)):
        f = float(sel.mean())
        out[key] = {"fraction": f, "n": int(sel.sum()),
                    "mc_sigma_fraction": sig(f)}
    out["central_quantiles"] = {
        "s1c_q16_q50_q84": np.quantile(s1, [.15865, .5, .84135]).tolist(),
        "log10_s2c_q16_q50_q84": np.quantile(y, [.15865, .5, .84135]).tolist()}
    return out


def acceptance(stats):
    out = {}
    for name, (nominal, tol) in GATE.items():
        f = stats[name]["fraction_inside"]
        out[name] = {"nominal": nominal, "tolerance_abs": tol,
                     "measured": f, "deviation_abs": f - nominal,
                     "pull_vs_mc_sigma": (f - nominal)
                     / stats[name]["mc_sigma_fraction"],
                     "pass": bool(abs(f - nominal) <= tol)}
    out["pass"] = all(v["pass"] for v in out.values())
    return out


def smooth2d(hist, sigmas):
    """Display-only separable Gaussian blur (sigmas in bins)."""
    out = hist.astype(float)
    for axis, sigma in enumerate(sigmas):
        k = np.arange(-int(np.ceil(4 * sigma)), int(np.ceil(4 * sigma)) + 1)
        kern = np.exp(-.5 * (k / sigma) ** 2)
        kern /= kern.sum()
        out = np.apply_along_axis(
            lambda v: np.convolve(v, kern, mode="same"), axis, out)
    return out


def hdr_levels(hist):
    flat = np.sort(hist.ravel())[::-1]
    cum = np.cumsum(flat) / flat.sum()
    return {f"q{int(100 * p)}": float(flat[np.searchsorted(cum, p)])
            for p in (.68, .95)}


def smoothing_bins(edges):
    dy = float(np.median(np.diff(edges[1])))
    dx = float(np.median(np.diff(edges[0])))
    return (SMOOTH_PHD / dx, SMOOTH_DEX / dy), dx, dy


def make_figure(s1, regions, stats, acc, hist, edges, out):
    sig_bins, dx, dy = smoothing_bins(edges)
    smooth = smooth2d(hist, sig_bins)
    lev = hdr_levels(smooth)
    cx = .5 * (edges[0][1:] + edges[0][:-1])
    cy = .5 * (edges[1][1:] + edges[1][:-1])
    cmap = LinearSegmentedColormap.from_list(
        "model_density", [SURFACE, "#f4d7bf", OUR])
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 7.0), sharex=True,
                             sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, lw=.6, alpha=.7)
        ax.set_axisbelow(True)
        ax.set_xlim(0, 80)
        ax.set_ylim(2.66, 4.5)
        ax.set_xlabel("S1c [phd]", color=FG)
        ax.set_ylabel("log10(S2c [phd])", color=FG)
        ax.tick_params(colors=FG)
        for spine in ax.spines.values():
            spine.set_color(GRID)
    axes[0].set_title("Our full simulation with LZ's 68% / 95% regions",
                      color=FG, fontsize=11)
    axes[1].set_title("Corresponding contours for our simulated sample",
                      color=FG, fontsize=11)
    mesh = axes[0].pcolormesh(edges[0], edges[1],
                              np.ma.masked_less_equal(hist.T, 0), cmap=cmap,
                              rasterized=True, shading="flat")
    fig.colorbar(mesh, ax=axes[0], pad=.02,
                 label=f"MC events per {dx:g} phd × {dy:g} dex")
    for label, ls, alpha in (("q95", "--", .12), ("q68", "-", .22)):
        poly = regions[label]
        for ax in axes:
            ax.fill(poly[:, 0], poly[:, 1], color=LZ, alpha=alpha, lw=0,
                    zorder=2)
            ax.plot(np.r_[poly[:, 0], poly[0, 0]],
                    np.r_[poly[:, 1], poly[0, 1]], color=LZ, lw=1.65, ls=ls,
                    zorder=4)
    axes[1].contour(cx, cy, smooth.T, levels=[lev["q95"], lev["q68"]],
                    colors=[OUR, OUR], linewidths=1.8,
                    linestyles=["--", "-"])
    ix, iy = np.unravel_index(np.argmax(hist), hist.shape)
    axes[0].plot(cx[ix], cy[iy], marker="+", color=FG, ms=10, mew=1.8)
    axes[0].text(cx[ix] + 1.0, cy[iy] + .018, "density peak", color=FG,
                 va="bottom", ha="left", fontsize=8.5)
    thr = np.log10(493.)
    for ax in axes:
        ax.axhline(thr, color=MUTED, lw=.8, ls=(0, (2, 3)))
        ax.text(77, thr + .008, "our S2c = 493 cut", color=MUTED,
                fontsize=8.5, ha="right", va="bottom")
    n = len(s1)
    verdict = lambda k: "PASS" if acc[k]["pass"] else "FAIL"
    text = (f"Gate: {'PASS' if acc['pass'] else 'FAIL'}  (68% region: 0.68 ± 0.05, "
            "95% region: 0.95 ± 0.03)\n"
            "Our ROI sample within LZ's published region\n"
            f"68% region: {100 * stats['q68']['fraction_inside']:.2f}% "
            f"± {100 * stats['q68']['mc_sigma_fraction']:.2f}% "
            f"({stats['q68']['n_inside']:,} / {n:,})  "
            f"[{verdict('q68')}]\n"
            f"95% region: {100 * stats['q95']['fraction_inside']:.2f}% "
            f"± {100 * stats['q95']['mc_sigma_fraction']:.2f}% "
            f"({stats['q95']['n_inside']:,} / {n:,})  "
            f"[{verdict('q95')}]")
    axes[1].text(.975, .955, text, transform=axes[1].transAxes, color=FG,
                 fontsize=9, ha="right", va="top",
                 bbox=dict(facecolor=SURFACE, alpha=.94, edgecolor=GRID,
                           boxstyle="round,pad=.45"))
    handles = [
        Line2D([], [], color=LZ, ls="-", lw=1.8, label="LZ 68% region"),
        Line2D([], [], color=LZ, ls="--", lw=1.8, label="LZ 95% region"),
        Line2D([], [], color=OUR, ls="-", lw=1.8,
               label="This work: 68% contour"),
        Line2D([], [], color=OUR, ls="--", lw=1.8,
               label="This work: 95% contour")]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               labelcolor=FG, fontsize=9, bbox_to_anchor=(.5, .955))
    fig.suptitle("40 GeV/c² WIMP in (S1c, log10 S2c) vs LZ WS2024 Fig. 3 "
                 "(All Data)", color=FG, fontsize=14, y=.993)
    fig.text(.5, .016,
             f"Left: all {int(hist.sum()):,} full-ROI MC events, raw counts, "
             "no reweighting. Right: smoothed density (1.1 phd, 0.040 dex), "
             "contours enclose 68% / 95% of our ROI sample. Containment: "
             f"unbinned, on the uniform {n:,}-event display subsample.\n"
             "LZ shading is extracted from its vector paths. Only the "
             "published enclosed regions are available, not LZ's full 2D "
             "density or its contour statistical uncertainty.",
             color=MUTED, fontsize=8.7, ha="center")
    fig.tight_layout(rect=[0, .075, 1, .93])
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--models", required=True)
    ap.add_argument("--paper", required=True)
    ap.add_argument("--output", required=True,
                    help=".png for the figure, .json for the metrics")
    args = ap.parse_args()
    suffix = Path(args.output).suffix
    assert suffix in (".png", ".json"), suffix

    regions = regions_from_pdf(args.paper)
    s1, y, hist, edges = load_model(args.models)
    stats = containment(s1, y, regions)
    acc = acceptance(stats)
    if suffix == ".png":
        make_figure(s1, regions, stats, acc, hist, edges, args.output)
        print(f"68%: {stats['q68']['fraction_inside']:.4f}  "
              f"95%: {stats['q95']['fraction_inside']:.4f}")
        return
    sig_bins, dx, dy = smoothing_bins(edges)
    report = {
        "scope": "40 GeV/c2 WIMP vs WS2024 Fig. 3 purple regions",
        "acceptance": acc,
        "acceptance_rule": ("|f68 - 0.68| <= 0.05 and |f95 - 0.95| <= 0.03 "
                            "(absolute fractions; tolerance, not a "
                            "significance test; user decision 2026-11-19)"),
        "reference": {
            "paper": "arXiv:2410.17036v3", "figure": "Fig. 3 top, All Data",
            "sha256": SHA256, "quote": CAPTION,
            "extraction": ("exact polygon vertices of the two purple fill "
                           "paths; tick-label calibration"),
            "vertices": {k: v.tolist() for k, v in regions.items()}},
        "our_model": {
            "archive_group": MASS_GROUP,
            "selection": "3 < S1c < 80 phd; S2c >= 493 phd; log10 S2c < 4.5",
            "n_full_roi_mc": int(hist.sum()),
            "n_containment_sample": len(s1),
            "containment_sample": ("uniform display subsample of the full "
                                   "ROI sample, without replacement, equal "
                                   "weights")},
        "containment": stats,
        "our_visual_contours": {
            "binning": f"{dx:g} phd x {dy:g} dex (archive grid)",
            "smoothing_sigma": f"{SMOOTH_PHD} phd x {SMOOTH_DEX} dex",
            "smoothing_sigma_bins": list(sig_bins),
            "levels": hdr_levels(smooth2d(hist, sig_bins)),
            "purpose": "visual guide only; containment is unbinned"},
        "scratch_probe_2026_11_13": {
            "note": ("same 59,852-event realization read event by event "
                     "from a /tmp archive no longer present"),
            "q68_fraction_inside": 0.6600447771168884,
            "q95_fraction_inside": 0.9437111541803115},
        "cautions": [
            "Published contours are regions, not a 2D PDF or its marginals.",
            "Containment is not a goodness-of-fit test; LZ's contour MC and "
            "construction uncertainties are unavailable.",
            "The MC sigma is binomial on the containment sample only.",
            "The low-S2 bound is the model's corrected-space S2c >= 493 phd, "
            "applied once; no raw-S2 alternative."]}
    Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["acceptance"], indent=2))


if __name__ == "__main__":
    main()
