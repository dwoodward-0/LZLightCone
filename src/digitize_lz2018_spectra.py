"""Digitize the per-source background energy spectra from the LZ 2018
sensitivity projection (arXiv:1802.06039): Fig. 4 (ER spectra, right
0-200 keV close-up panel) and Fig. 5 (NR spectra).

These are the SHAPE inputs for the ER and NR background components
(decision background_components: all_table1); each is sampled in recoil
energy and pushed through the response model. Normalization comes from
the WS2024 Table I expected counts, not from the curves' absolute rates
(the 2018 paper's assumptions differ from WS2024's). The plotted S1c/S2c
distributions of the 2018 paper are NOT used (different field
configuration).

Runs LOCALLY (not under lc): provenance generation for a data/ input,
the same pattern as src/digitize_neutron_spectrum.py.

Color -> curve mapping (established by rendering both panels and reading
the directly-placed labels against line styles; the palette is shared by
the two figures):
  Fig. 4 right panel (0-200 keV ER close-up), smooth polylines:
    navy  solid        222Rn chain  (214Pb beta; flat ~2.1e-5 in panel)
    olive  dashed      220Rn chain  (212Pb beta)
    teal   solid       85Kr beta (flat ~1e-6, lowest)
    orange dash-dot    136Xe 2nu-beta-beta (rising 1e-6 -> 1e-4)
    green  solid       Solar nu (free electron) (RRPA)
    pink   step hist   Det. + Sur. + Env. ER (flat ~2e-6, rises >155 keV)
    black  step hist   Total (not extracted)
  Fig. 5 (0-100 keV NR), step histograms:
    navy  solid        8B nu NR (steep, 1e-3 -> 2e-8, ends ~5 keV)
    teal   solid       hep nu NR (1e-5 -> 1e-10, ends ~7 keV)
    orange dash-dot    Atm nu NR (1e-8 at 2 keV -> 1e-11 at 87 keV)
    olive  dashed      DSN nu NR (not a WS2024 Table I row; not extracted)
    pink   step hist   Det. + Sur. + Env. NR (already digitized separately)
    black  step hist   Total (not extracted)

Solar nu subtlety: the RRPA note states that below 30 keV the
pp+7Be+13N contribution is scaled by the RRPA. In the right panel the
green solid curve below ~25 keV is a SEPARATE path (suppressed rate,
with a visible step at the joining energy) from a dashed path at
slightly higher rate; the dashed piece is the free-electron below-30-keV
continuation (it connects smoothly onto the main solid path). The
published RRPA curve is therefore the SOLID path pair (suppressed below
~25 keV), which is what this script extracts as er_solar_nu; the dashed
free-electron piece is extracted as er_solar_nu_freeelow for reference.

Figure facts (established by inspection of the vector drawings):
- Fig. 4 (page index 7): right panel plot box x ~ 340.5-540.4 pt,
  y ~ 260-393 pt. x labels 0/50/100/150/200 keV; y log decade labels
  1e-4 / 1e-5 / 1e-6 (panel top is 1e-3, clipped in the render).
  Smooth curves are single connected paths (~91 segments); the green
  curve is three paths (solid low-E piece, dashed low-E piece, solid
  main); Det.+Sur.+Env. and Total are step histograms (~91 segments).
- Fig. 5 (page index 9): plot area x ~ 349-547.6, y ~ 69-211 pt;
  y log majors 1e-3..1e-11 every 17.7875 pt; x majors 0-100 keV every
  20 keV (same calibration as src/digitize_neutron_spectrum.py).
"""

import csv
import hashlib
import re
from pathlib import Path

import fitz  # pymupdf
import numpy as np

PDF = Path(__file__).resolve().parent.parent / "data" / "lz_sensitivity2018_paper.pdf"
OUT = Path(__file__).resolve().parent.parent / "data" / "lz2018_component_spectra.csv"
DIAGNOSTIC = Path("/tmp/lz2018_component_spectra_digitization.png")

SHA256 = "5913d97d2d45710e4a6c7fddf8aac3a98f20adc3be3faec7c553a98611f9c058"

NAVY = (0.06666667014360428, 0.2666666805744171, 0.46666666865348816)
OLIVE = (0.6, 0.6, 0.2)
TEAL = (0.2666666805744171, 0.6666666865348816, 0.6)
ORANGE = (0.9764705896377563, 0.34117648005485535, 0.21960784494876862)
PINK = (0.8470588326454163, 0.06666667014360428, 0.3490196168422699)
GREEN = (0.06666667014360428, 0.46666666865348816, 0.2)
BLACK_TOL = 0.05
COLOR_TOL = 0.02


def close(c, ref, tol=COLOR_TOL):
    return c is not None and len(c) >= 3 and all(
        abs(a - b) < tol for a, b in zip(c[:3], ref)
    )


def paths_of_color(page, color, box):
    """Stroke paths of `color`; returns a list of (n_segments, [points]).
    Only the line segments whose BOTH endpoints lie inside `box` are kept,
    so closure verticals dipping below the plot floor are dropped without
    discarding the rest of the path."""
    x0, x1, y0, y1 = box
    out = []
    for d in page.get_drawings():
        if d["type"] != "s" or not close(d.get("color"), color):
            continue
        pts = []
        for item in d["items"]:
            if item[0] != "l":
                continue
            p, q = item[1], item[2]
            if all(x0 <= v.x <= x1 and y0 <= v.y <= y1 for v in (p, q)):
                pts.extend([(p.x, p.y), (q.x, q.y)])
        if len(pts) >= 4:
            out.append((len(pts), pts))
    return out


def fit(vals, refs, tol, what):
    f = np.polyfit(vals, refs, 1)
    resid = np.abs(np.polyval(f, vals) - refs).max()
    assert resid < tol, f"{what}: calibration residual {resid:.4g}"
    return f


def decade_labels(page, x_min, x_max, y_min, y_max):
    """(label_y, exponent) for '10^-k' words whose x span sits inside
    [x_min, x_max] (the axis-label column) and y inside the band."""
    anchors = []
    for w in page.get_text("words"):
        m = re.fullmatch(r"10[−-](\d+)", w[4])
        if m and x_min < w[0] and w[2] < x_max and y_min < w[1] < y_max:
            cy = 0.5 * (w[1] + w[3])
            anchors.append((cy, -int(m.group(1))))
    return anchors


def calibrate_fig4_axes(page):
    """Right-panel close-up: x from the 0..200 keV labels, y from the log
    decade labels (anchored directly to the label centers; residuals are
    asserted small)."""
    # x: label words 0/50/100/150/200 just below the panel
    xlab = {}
    for w in page.get_text("words"):
        if w[4] in {"0", "50", "100", "150", "200"} and 396 < w[1] < 408 and w[0] > 330:
            xlab.setdefault(int(w[4]), 0.5 * (w[0] + w[2]))
    assert sorted(xlab) == [0, 50, 100, 150, 200], f"x labels found: {xlab}"
    xf = fit([xlab[k] for k in sorted(xlab)], [0, 50, 100, 150, 200], 0.3, "fig4 x")

    # y: the panel's own decade labels sit at x ~ 324-338 (1e-4, 1e-5, 1e-6
    # only; words further left belong to the left panel / caption)
    anchors = decade_labels(page, 320, 340, 255, 400)
    assert sorted(e for _, e in anchors) == [-6, -5, -4], (
        f"fig4 y labels unexpected: {anchors}"
    )
    yf = fit([a[0] for a in anchors], [a[1] for a in anchors], 0.02, "fig4 y")
    return xf, yf


def smooth_curve(page, box, color, xf, yf, exclude_y=None):
    """Concatenated (E, rate) samples of the smooth polyline(s) of `color`.

    exclude_y: optional (ylo, yhi) band whose paths are dropped (used to
    drop the dashed free-electron low-E piece of the green curve).
    """
    segs = []
    for n, pts in paths_of_color(page, box=box, color=color):
        if exclude_y is not None:
            ys = [p[1] for p in pts]
            if exclude_y[0] <= min(ys) and max(ys) <= exclude_y[1]:
                continue
        segs.append(pts)
    assert segs, f"no paths for color {color}"
    # concatenate by x order (a curve split across paths)
    pts = [p for seg in segs for p in seg]
    pts = sorted(pts)  # by x, then y
    e = np.polyval(xf, [p[0] for p in pts])
    rate = 10.0 ** np.polyval(yf, [p[1] for p in pts])
    # collapse duplicate x (path joins)
    keep = np.concatenate([[True], np.diff(e) > 1e-6])
    return e[keep], rate[keep]


def step_curve(page, box, color, xf, yf, min_pts=8, floor_pts=None):
    """(bin edges, rates) of the step histogram of `color` (neutron-
    digitizer idiom: classify segments H/V). floor_pts is the y coordinate
    of the plot floor (the lowest y-major tick) — needed to recognize the
    floor-level closure bin that terminates a physically-ending spectrum;
    None disables the trailing-bin drop."""
    cands = paths_of_color(page, box=box, color=color)
    cands = [(n, p) for n, p in cands if n >= min_pts]
    assert len(cands) == 1, (
        f"expected one multi-segment step path for {color}, got {len(cands)}"
    )
    pts = np.array(cands[0][1], dtype=float)
    bins = []
    for i in range(len(pts) - 1):
        dx = abs(pts[i + 1, 0] - pts[i, 0])
        dy = abs(pts[i + 1, 1] - pts[i, 1])
        if dx > dy:  # horizontal segment = one bin
            elo = np.polyval(xf, min(pts[i, 0], pts[i + 1, 0]))
            ehi = np.polyval(xf, max(pts[i, 0], pts[i + 1, 0]))
            lgr = np.polyval(yf, 0.5 * (pts[i, 1] + pts[i + 1, 1]))
            bins.append((elo, ehi, lgr, 0.5 * (pts[i, 1] + pts[i + 1, 1])))
    bins = np.array(bins)
    # a spectrum that ENDS is drawn with its final vertical dropping to the
    # plot floor plus one floor-level horizontal closure bin (e.g. 8B, which
    # physically ends at the 8B endpoint ~4.1 keVnr); drop trailing bins at
    # the floor level when the previous bin is a cliff above them (a smooth
    # genuine tail as in Atm never satisfies both conditions)
    while (
        len(bins) >= 2
        and bins[-1, 2] < bins[-2, 2] - 1.0  # >10x below previous bin
        and floor_pts is not None
        and abs(bins[-1, 3] - floor_pts) < 1.5  # drawn at the plot floor
    ):
        bins = bins[:-1]
    gaps = bins[1:, 0] - bins[:-1, 1]
    assert np.abs(gaps).max() < 0.06, f"bins not contiguous: max gap {gaps.max():.3f}"
    return bins[:, :2], 10.0 ** bins[:, 2]


def calibrate_fig5_axes(page):
    """Same calibration as digitize_neutron_spectrum.py (identical panel)."""
    y_majors, x_majors = [], []
    for d in page.get_drawings():
        if not close(d.get("color"), (0.0, 0.0, 0.0), BLACK_TOL):
            continue
        for item in d["items"]:
            if item[0] != "l":
                continue
            p, q = item[1], item[2]
            length = abs(complex(q.x - p.x, q.y - p.y))
            if not (3.0 < length < 4.5):
                continue
            if abs(p.y - q.y) < 0.3 and 340 < min(p.x, q.x) < 350:
                y_majors.append((p.y + q.y) / 2)
            if abs(p.x - q.x) < 0.3 and 200 < min(p.y, q.y) and max(p.y, q.y) < 215:
                x_majors.append((p.x + q.x) / 2)
    y_majors, x_majors = sorted(y_majors), sorted(x_majors)
    assert len(y_majors) == 9, f"fig5: expected 9 y majors, got {len(y_majors)}"
    assert len(x_majors) == 6, f"fig5: expected 6 x majors, got {len(x_majors)}"
    anchors = []
    for w in page.get_text("words"):
        m = re.fullmatch(r"10[−-](\d+)", w[4])
        if m and w[2] < 350 and w[1] > 60:
            cy = 0.5 * (w[1] + w[3])
            nearest = min(y_majors, key=lambda t: abs(t - cy))
            if abs(nearest - cy) < 3.0:
                anchors.append((nearest, -int(m.group(1))))
    assert len(anchors) == 9, f"fig5: expected 9 decade labels, got {len(anchors)}"
    ys = np.array([t for t, _ in anchors])
    vs = np.array([v for _, v in anchors])
    yf = np.polyfit(ys, vs, 1)
    assert np.abs(np.polyval(yf, ys) - vs).max() < 0.02, "fig5 y residual"
    xf = np.polyfit(x_majors, np.arange(0.0, 101.0, 20.0), 1)
    return xf, yf, max(y_majors)  # plot floor = lowest y-major (1e-11)


def main():
    digest = hashlib.sha256(PDF.read_bytes()).hexdigest()
    assert digest == SHA256, f"PDF sha256 {digest} != pinned {SHA256}"
    doc = fitz.open(PDF)

    rows = []  # (figure, component, kind, e_lo, e_hi, log10_rate)

    # ---------- Fig. 4 right panel (ER close-up) ----------
    page4 = doc[7]
    assert "ER background spectra" in page4.get_text(), "Fig. 4 caption not found"
    box4 = (335, 545, 255, 397)
    xf4, yf4 = calibrate_fig4_axes(page4)
    print(f"fig4 x fit: E = {xf4[0]:.6f}*pt + {xf4[1]:.3f}")
    print(f"fig4 y fit: log10(rate) = {yf4[0]:.6f}*pt + {yf4[1]:.4f}")

    smooth = {
        "er_pb214": NAVY,
        "er_pb212": OLIVE,
        "er_kr85": TEAL,
        "er_xe136": ORANGE,
    }
    for comp, color in smooth.items():
        e, r = smooth_curve(page4, box4, color, xf4, yf4)
        assert np.all(np.diff(e) > 0), f"{comp}: E not monotone"
        assert e[0] < 1.0 and e[-1] > 195.0, f"{comp}: E span [{e[0]}, {e[-1]}]"
        for ei, ri in zip(e, r):
            rows.append(("fig4", comp, "sample", ei, ei, np.log10(ri)))
        print(f"{comp:16s}: {len(e)} samples, E [{e[0]:.2f}, {e[-1]:.2f}], "
              f"rate [{r.min():.3e}, {r.max():.3e}]")

    # solar nu: solid RRPA curve = solid low-E path + main path; drop the
    # dashed free-electron low-E piece (dashes '[ .475 0.784 ]')
    segs_green = paths_of_color(page4, box=box4, color=GREEN)
    assert len(segs_green) == 3, f"expected 3 green paths, got {len(segs_green)}"
    solid, dashed = [], []
    low = [pts for n, pts in segs_green if max(p[0] for p in pts) <= 367.5]
    main = [pts for n, pts in segs_green if max(p[0] for p in pts) > 367.5]
    assert len(low) == 2 and len(main) == 1, "green path split unexpected"
    # both low-E paths span E ~ 0-25; the RRPA (solid) one sits at LARGER
    # y (lower rate) than the dashed free-electron one
    lo_hi = min(low, key=lambda pts: np.mean([p[1] for p in pts]))   # higher rate
    lo_rrpa = max(low, key=lambda pts: np.mean([p[1] for p in pts]))  # lower rate
    pts = sorted(lo_rrpa + main[0])
    e = np.polyval(xf4, [p[0] for p in pts])
    r = 10.0 ** np.polyval(yf4, [p[1] for p in pts])
    keep = np.concatenate([[True], np.diff(e) > 1e-6])
    e, r = e[keep], r[keep]
    for ei, ri in zip(e, r):
        rows.append(("fig4", "er_solar_nu_rrpa", "sample", ei, ei, np.log10(ri)))
    print(f"{'er_solar_nu_rrpa':16s}: {len(e)} samples, rate [{r.min():.3e}, {r.max():.3e}]")
    pts = sorted(lo_hi + main[0])
    e = np.polyval(xf4, [p[0] for p in pts])
    r = 10.0 ** np.polyval(yf4, [p[1] for p in pts])
    keep = np.concatenate([[True], np.diff(e) > 1e-6])
    e, r = e[keep], r[keep]
    for ei, ri in zip(e, r):
        rows.append(("fig4", "er_solar_nu_free", "sample", ei, ei, np.log10(ri)))
    print(f"{'er_solar_nu_free':16s}: {len(e)} samples, rate [{r.min():.3e}, {r.max():.3e}]")

    # Det.+Sur.+Env. ER: step histogram
    edges, rates = step_curve(page4, box4, PINK, xf4, yf4)
    for (elo, ehi), ri in zip(edges, rates):
        rows.append(("fig4", "er_det_gammas", "bin", elo, ehi, np.log10(ri)))
    print(f"{'er_det_gammas':16s}: {len(rates)} bins, E [{edges[0,0]:.2f}, "
          f"{edges[-1,1]:.2f}], rate [{rates.min():.3e}, {rates.max():.3e}]")

    # shape sanity vs the rendered panel
    def rate_at(comp, e_query):
        sub = [(lo, r) for fig, c, kind, lo, hi, lr in rows
               if c == comp for r in [10.0 ** lr] for lo in [lo]]
        sub = sorted(set(sub))
        es = np.array([s[0] for s in sub]); rs = np.array([s[1] for s in sub])
        i = np.argmin(np.abs(es - e_query))
        return rs[i]

    assert abs(np.log10(rate_at("er_pb214", 50)) + 4.68) < 0.08, \
        f"222Rn at 50 keV = {rate_at('er_pb214', 50):.2e} (expect ~2.1e-5)"
    assert abs(np.log10(rate_at("er_kr85", 50)) + 6.0) < 0.1, \
        f"85Kr at 50 keV = {rate_at('er_kr85', 50):.2e} (expect ~1e-6)"
    assert rate_at("er_xe136", 180) > 5 * rate_at("er_xe136", 20), \
        "136Xe should rise across the panel"

    # ---------- Fig. 5 (NR) ----------
    page5 = doc[9]
    assert "NR background spectra" in page5.get_text(), "Fig. 5 caption not found"
    box5 = (345, 552, 60, 212)
    xf5, yf5, floor5 = calibrate_fig5_axes(page5)
    print(f"\nfig5 x fit: E = {xf5[0]:.6f}*pt + {xf5[1]:.3f}")
    print(f"fig5 y fit: log10(rate) = {yf5[0]:.6f}*pt + {yf5[1]:.4f}")

    for comp, color in {"nr_b8": NAVY, "nr_hep": TEAL, "nr_atm": ORANGE}.items():
        edges, rates = step_curve(page5, box5, color, xf5, yf5, min_pts=8,
                                  floor_pts=floor5)
        for (elo, ehi), ri in zip(edges, rates):
            rows.append(("fig5", comp, "bin", elo, ehi, np.log10(ri)))
        print(f"{comp:16s}: {len(rates)} bins, E [{edges[0,0]:.2f}, "
              f"{edges[-1,1]:.2f}], rate [{rates.min():.3e}, {rates.max():.3e}]")

    assert rows and all(np.isfinite(r[5]) for r in rows)

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "# Per-source background energy spectra, digitized from the LZ 2018",
            "# sensitivity projection, arXiv:1802.06039 (in-project copy",
            "# data/lz_sensitivity2018_paper.pdf, sha256 " + SHA256 + ").",
            "# Fig. 4 right panel (0-200 keV electronic-recoil close-up, page 8):",
            "#   er_pb214 (navy solid, '222Rn'), er_pb212 (olive dashed, '220Rn'),",
            "#   er_kr85 (teal solid, '85Kr'), er_xe136 (orange dash-dot, '136Xe'),",
            "#   er_solar_nu_rrpa (green solid incl. its suppressed below-~25-keV",
            "#     RRPA piece; see the RRPA note), er_solar_nu_free (the dashed",
            "#     free-electron below-30-keV piece spliced onto the same main",
            "#     path, kept for reference), er_det_gammas (pink steps,",
            "#     'Det. + Sur. + Env.').",
            "# Fig. 5 (0-100 keV nuclear-recoil, page 10):",
            "#   nr_b8 (navy steps, '8B'), nr_hep (teal steps, 'hep'),",
            "#   nr_atm (orange dash-dot steps, 'Atm').",
            "# Note: nr_b8's first bin (0.14-1.14 keV) is drawn clipped at the",
            "# axis top (1e-3); its value is a LOWER BOUND. Harmless for the",
            "# model: that bin lies far below the NR efficiency turn-on.",
            "# Extracted programmatically from the vector graphics (pymupdf",
            "# page.get_drawings); axes calibrated from tick marks anchored to",
            "# the axis labels. Color->label mapping verified against a rendered",
            "# overlay of both panels (line styles: dashed=220Rn/DSN,",
            "# dash-dot=136Xe/Atm; the 8B/hep curves end where their labels sit).",
            "# kind 'sample': one row per polyline vertex, rate at e_lo=e_hi.",
            "# kind 'bin': step histogram, rate constant across [e_lo, e_hi).",
            "# Rates are counts/kg/day/keV as plotted (5.6-t fiducial volume,",
            "# single-scatter, no detector efficiency or ROI cuts applied).",
            "# SHAPE inputs for the WS2024 background model (decision",
            "# background_components: all_table1): sampled in recoil energy and",
            "# pushed through the response model; normalized to the WS2024",
            "# Table I counts, NOT to these curves' absolute rates.",
            "# Provenance: src/digitize_lz2018_spectra.py.",
        ])
        w.writerow(["figure", "component", "kind", "e_lo_kev", "e_hi_kev",
                    "log10_rate", "rate"])
        for fig, comp, kind, elo, ehi, lgr in rows:
            w.writerow([fig, comp, kind, f"{elo:.4f}", f"{ehi:.4f}",
                        f"{lgr:.5f}", f"{10.0 ** lgr:.4e}"])
    print(f"\nwrote {OUT} ({len(rows)} rows)")

    # diagnostic: extracted curves over the rendered panels + data coords
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
        clip = fitz.Rect(330, 250, 550, 415)
        pix = page4.get_pixmap(clip=clip, dpi=200)
        axes[0].imshow(
            np.frombuffer(pix.samples, dtype=np.uint8).reshape(
                pix.height, pix.width, pix.n)[:, :, :3],
            extent=(clip.x0, clip.x1, clip.y1, clip.y0),
        )
        axes[0].set_title("Fig. 4 right panel (render)")
        clip5 = fitz.Rect(340, 60, 555, 220)
        pix = page5.get_pixmap(clip=clip5, dpi=200)
        axes[1].imshow(
            np.frombuffer(pix.samples, dtype=np.uint8).reshape(
                pix.height, pix.width, pix.n)[:, :, :3],
            extent=(clip5.x0, clip5.x1, clip5.y1, clip5.y0),
        )
        axes[1].set_title("Fig. 5 (render)")
        styles = {
            "er_pb214": ("navy", "-"), "er_pb212": ("olive", "--"),
            "er_kr85": ("teal", "-"), "er_xe136": ("orangered", "-."),
            "er_solar_nu_rrpa": ("green", "-"), "er_solar_nu_free": ("limegreen", ":"),
            "er_det_gammas": ("crimson", "-"),
        }
        for comp, (col, ls) in styles.items():
            sub = [(lo, hi, 10.0 ** lr) for f_, c, k, lo, hi, lr in rows if c == comp]
            if not sub:
                continue
            if sub[0][0] != sub[0][1]:
                edges = [s[0] for s in sub] + [sub[-1][1]]
                vals = [s[2] for s in sub]
                axes[2].step(edges, vals + [vals[-1]], where="post", color=col,
                             ls=ls, lw=1, label=comp)
            else:
                es = [s[0] for s in sub]; rs = [s[2] for s in sub]
                axes[2].plot(es, rs, color=col, ls=ls, lw=1, label=comp)
        for comp, col in {"nr_b8": "navy", "nr_hep": "teal", "nr_atm": "orangered"}.items():
            sub = [(lo, hi, 10.0 ** lr) for f_, c, k, lo, hi, lr in rows if c == comp]
            edges = [s[0] for s in sub] + [sub[-1][1]]
            vals = [s[2] for s in sub]
            axes[2].step(edges, vals + [vals[-1]], where="post", color=col,
                         ls="-", lw=1, label=comp)
        axes[2].set_yscale("log")
        axes[2].set_xlabel("recoil energy (keV)")
        axes[2].set_ylabel("rate [counts/kg/day/keV]")
        axes[2].legend(fontsize=6, ncol=2)
        axes[2].set_title("digitized component spectra")
        fig.tight_layout()
        fig.savefig(DIAGNOSTIC, dpi=130)
        print(f"diagnostic plot: {DIAGNOSTIC}")
    except Exception as e:  # noqa: BLE001 - diagnostic is best-effort
        print(f"diagnostic plot failed: {e}")


if __name__ == "__main__":
    main()
