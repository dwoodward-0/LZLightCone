"""Digitize the projected detector+surface+environment neutron NR spectrum
(Fig. 5 of arXiv:1802.06039, the LZ 2018 sensitivity projection).

Extracts the pink "Det. + Sur. + Env." step histogram from Fig. 5 via
page.get_drawings(), calibrates both axes from the tick MARKS (anchored to
the decade tick labels), and writes data/lz_neutron_spectrum_2018.csv: the
NR recoil-energy spectrum in counts/kg/day/keV, one row per histogram bin.

This is the SHAPE input for the detector-neutron background component
(decision background_components: all_table1): it is sampled in recoil
energy and pushed through the response model — the S1c/S2c distributions
plotted in the 2018 paper are NOT used (different field configuration).
Normalization comes from the WS2024 Table I expected neutron counts, not
from this curve's absolute rate.

Runs LOCALLY (not under lc): provenance generation for a data/ input, the
same pattern as src/digitize_ws2022_bands.py.

Figure facts (established by prior inspection of this page):
- Fig. 5 is on page index 9, plot area x ~ 349-547.6, y ~ 69-211 (PDF pts).
- The pink curve is the only stroked path in color (0.847, 0.067, 0.349)
  with more than a handful of segments (the second pink drawing is the
  label leader line). It is a closed step histogram: baseline closure
  segments at both ends, then alternating horizontal (bin) / vertical
  (edge) segments; ~98 bins tiling 0-100 keV.
- Y axis: log scale, major ticks every 17.7875 pt, decade labels 10^-3
  (y ~ 69) through 10^-11 (y ~ 211); all nine anchor their major ticks.
- X axis: linear, major ticks (3.56 pt) at 0/20/40/60/80/100 keV, minors
  every 5 keV.
"""

import csv
import hashlib
import re
from pathlib import Path

import fitz  # pymupdf
import numpy as np

PDF = Path(__file__).resolve().parent.parent / "data" / "lz_sensitivity2018_paper.pdf"
OUT = Path(__file__).resolve().parent.parent / "data" / "lz_neutron_spectrum_2018.csv"
DIAGNOSTIC = Path("/tmp/lz2018_fig5_digitization.png")

SHA256 = "5913d97d2d45710e4a6c7fddf8aac3a98f20adc3be3faec7c553a98611f9c058"
CAPTION_SNIPPET = "NR background spectra in the 5.6-tonne"
PAGE = 9
PINK = (0.8470588326454163, 0.06666667014360428, 0.3490196168422699)
COLOR_TOL = 0.02
BLACK_TOL = 0.05


def close(c, ref, tol=COLOR_TOL):
    return c is not None and len(c) >= 3 and all(
        abs(a - b) < tol for a, b in zip(c[:3], ref)
    )


def collect_ticks(page):
    """Major tick marks on both axes, from the black stroked segments."""
    y_majors, x_majors = [], []
    for d in page.get_drawings():
        if not close(d.get("color"), (0.0, 0.0, 0.0), BLACK_TOL):
            continue
        for item in d["items"]:
            if item[0] != "l":
                continue
            p, q = item[1], item[2]
            length = abs(complex(q.x - p.x, q.y - p.y))
            if not (3.0 < length < 4.5):  # majors; minors are 1.78 pt
                continue
            if abs(p.y - q.y) < 0.3 and 340 < min(p.x, q.x) < 350:
                y_majors.append((p.y + q.y) / 2)
            if abs(p.x - q.x) < 0.3 and 200 < min(p.y, q.y) and max(p.y, q.y) < 215:
                x_majors.append((p.x + q.x) / 2)
    return sorted(y_majors), sorted(x_majors)


def calibrate_y(page, y_majors):
    """Anchor the major ticks to the visible 10^-8..10^-11 decade labels."""
    assert len(y_majors) == 9, f"expected 9 y majors, got {len(y_majors)}"
    gaps = np.diff(y_majors)
    assert np.all(np.abs(gaps - gaps.mean()) < 0.05), f"uneven y majors: {gaps}"
    anchors = []
    for w in page.get_text("words"):
        m = re.fullmatch(r"10[−-](\d+)", w[4])
        if m and w[2] < 350 and w[1] > 60:
            cy = 0.5 * (w[1] + w[3])
            nearest = min(y_majors, key=lambda t: abs(t - cy))
            if abs(nearest - cy) < 3.0:
                anchors.append((nearest, -int(m.group(1))))
    assert len(anchors) == 9, f"expected 9 decade labels, matched {len(anchors)}"
    exps = sorted(e for _, e in anchors)
    assert exps == list(range(-11, -2)), f"unexpected label decades {exps}"
    ys = np.array([t for t, _ in anchors])
    vs = np.array([v for _, v in anchors])
    yf = np.polyfit(ys, vs, 1)
    resid = np.polyval(yf, ys) - vs
    assert np.abs(resid).max() < 0.02, f"y calibration residuals {resid}"
    return yf


def calibrate_x(x_majors):
    """x majors are 0..100 keV in 20-keV steps (labels 0/20/40/60/80/100)."""
    assert len(x_majors) == 6, f"expected 6 x majors, got {len(x_majors)}"
    gaps = np.diff(x_majors)
    assert np.all(np.abs(gaps - gaps.mean()) < 0.05), f"uneven x majors: {gaps}"
    xf = np.polyfit(x_majors, np.arange(0.0, 101.0, 20.0), 1)
    return xf


def extract_pink(page, xf):
    """Ordered vertices of the pink step histogram, baseline closure dropped."""
    candidates = []
    for d in page.get_drawings():
        if d["type"] != "s" or not close(d.get("color"), PINK):
            continue
        n = sum(1 for item in d["items"] if item[0] == "l")
        if n > 10:
            candidates.append((n, d))
    assert len(candidates) == 1, (
        f"expected one multi-segment pink path, got {len(candidates)}"
    )
    d = candidates[0][1]
    seq = []
    for item in d["items"]:
        if item[0] != "l":
            raise ValueError(f"unexpected drawing item {item[0]!r}")
        p = (round(item[1].x, 3), round(item[1].y, 3))
        q = (round(item[2].x, 3), round(item[2].y, 3))
        if not seq or seq[-1] != p:
            seq.append(p)
        seq.append(q)
    # drop the baseline closure: leading vertical from the floor, trailing
    # vertical back to the floor (both at the path's extreme y)
    y_floor = max(seq[0][1], seq[-1][1])
    assert seq[0][0] == seq[1][0], "path does not open with a vertical segment"
    assert seq[-2][0] == seq[-1][0], "path does not close with a vertical segment"
    assert seq[0][1] == y_floor and seq[-1][1] == y_floor, "closure not at floor"
    seq = seq[1:-1]
    pts = np.array(seq, dtype=float)
    e = np.polyval(xf, pts[:, 0])
    return e, pts[:, 1], pts


def bins_from_steps(pts, xf, yf):
    """Horizontal step segments -> contiguous (e_lo, e_hi, log10 rate) bins.

    Segments are classified H/V in PDF point space (H: |dx| > |dy|) — the
    figure occasionally draws two adjacent bins as one slightly-sloped
    segment, and near-equal bins can carry a sub-pixel vertical.
    """
    lg = np.polyval(yf, pts[:, 1])
    bins = []
    for i in range(len(pts) - 1):
        dx = abs(pts[i + 1, 0] - pts[i, 0])
        dy = abs(pts[i + 1, 1] - pts[i, 1])
        if dx > dy:  # horizontal: one bin (rate at segment midpoint)
            e = np.sort([np.polyval(xf, pts[i, 0]), np.polyval(xf, pts[i + 1, 0])])
            bins.append((e[0], e[1], 0.5 * (lg[i] + lg[i + 1])))
    bins = np.array(bins)
    # the path starts at the left frame and ends at the right frame
    gaps = bins[1:, 0] - bins[:-1, 1]
    assert np.abs(gaps).max() < 0.05, f"bins not contiguous: max gap {gaps.max()}"
    return bins


def main():
    digest = hashlib.sha256(PDF.read_bytes()).hexdigest()
    assert digest == SHA256, f"PDF sha256 {digest} != pinned {SHA256}"
    doc = fitz.open(PDF)
    page = doc[PAGE]
    assert CAPTION_SNIPPET in page.get_text(), "Fig. 5 caption not on this page"

    y_majors, x_majors = collect_ticks(page)
    yf = calibrate_y(page, y_majors)
    xf = calibrate_x(x_majors)
    print(f"y fit: log10(rate) = {yf[0]:.6f}*pt + {yf[1]:.4f}  "
          f"({-1 / yf[0]:.3f} pt/decade)")
    print(f"x fit: E = {xf[0]:.6f}*pt + {xf[1]:.4f}")
    assert abs(np.polyval(xf, 349.04)) < 0.2 and abs(np.polyval(xf, 547.57) - 100) < 0.2
    assert abs(np.polyval(yf, 157.87) + 8) < 0.02

    e, y, pts = extract_pink(page, xf)
    bins = bins_from_steps(pts, xf, yf)
    rate = 10.0 ** bins[:, 2]
    print(f"\n{len(bins)} bins, E [{bins[0, 0]:.2f}, {bins[-1, 1]:.2f}] keV")
    print(f"log10(rate) [{bins[:, 2].min():.3f}, {bins[:, 2].max():.3f}]")
    slope = np.polyfit(bins[:, 0], bins[:, 2], 1)[0]
    print(f"overall log10(rate) slope {slope:.4f}/keV (expect negative)")
    assert slope < 0
    assert bins[0, 0] < 0.2 and abs(bins[-1, 1] - 100) < 0.2
    total = np.sum(rate * np.diff(bins[:, :2], axis=1)[0])
    print(f"integrated rate {total:.3e} counts/kg/day "
          f"({total * 365.25 * 1000:.1f} counts/t/y before cuts)")

    with open(OUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "# Projected detector + surface + environment neutron NR energy",
                "# spectrum, digitized from Fig. 5 (pink 'Det. + Sur. + Env.' curve)",
                "# of arXiv:1802.06039 (LZ 2018 sensitivity projection), page 10.",
                "# Single-scatter events in the 5.6-tonne fiducial volume with",
                "# neither a xenon skin nor an OD veto signal; no detector",
                "# efficiency or ROI (S1c) cuts applied (figure caption).",
                "# Extracted programmatically from the vector graphics (pymupdf",
                "# page.get_drawings); axes calibrated from the major tick marks,",
                "# anchored to the 10^-8..10^-11 decade labels (y) and 0-100 keV",
                "# labels (x). Step histogram: each row is one bin, rate constant",
                "# across [e_lo, e_hi). Bins are 1 keV (a few 2 keV at high E).",
                "# Columns: e_lo_kev, e_hi_kev, log10_rate, rate",
                "# [counts/kg/day/keV]. SHAPE input for the detector-neutron",
                "# background component; normalization comes from the WS2024",
                "# Table I expected counts, not from this curve's absolute rate.",
                "# Provenance: src/digitize_neutron_spectrum.py; in-project copy",
                "# data/lz_sensitivity2018_paper.pdf, sha256 " + SHA256 + ".",
            ]
        )
        writer.writerow(["e_lo_kev", "e_hi_kev", "log10_rate", "rate"])
        for lo, hi, l in bins:
            writer.writerow([f"{lo:.3f}", f"{hi:.3f}", f"{l:.4f}",
                             f"{10.0 ** l:.4e}"])
    print(f"\nwrote {OUT}")

    # diagnostic render: extracted curve over the page raster + data coords
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        clip = fitz.Rect(310, 55, 570, 245)
        pix = page.get_pixmap(clip=clip, dpi=200)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        axes[0].imshow(
            np.frombuffer(pix.samples, dtype=np.uint8).reshape(
                pix.height, pix.width, pix.n)[:, :, :3],
            extent=(clip.x0, clip.x1, clip.y1, clip.y0),
        )
        e_step = np.repeat(pts[:, 0], 2)[1:-1]  # PDF-space x
        y_step = np.repeat(pts[:, 1], 2)[1:-1]
        axes[0].plot(e_step, y_step, "-", color="lime", lw=0.7)
        axes[0].set_title("extracted pink curve over Fig. 5 raster (PDF coords)")
        axes[1].step(np.append(bins[:, 0], bins[-1, 1]),
                     np.append(rate, rate[-1]),
                     where="post", color="crimson", lw=0.8)
        axes[1].set_yscale("log")
        axes[1].set_xlabel("NR recoil energy (keV)")
        axes[1].set_ylabel("rate [counts/kg/day/keV]")
        axes[1].set_title("digitized Det.+Sur.+Env. neutron spectrum")
        fig.tight_layout()
        fig.savefig(DIAGNOSTIC, dpi=130)
        print(f"diagnostic plot: {DIAGNOSTIC}")
    except Exception as e:  # noqa: BLE001 - diagnostic is best-effort
        print(f"diagnostic plot failed: {e}")


if __name__ == "__main__":
    main()
