"""Digitize the SR1 accidentals-PDF marginal functions (Fig. 12 of
arXiv:2211.17120, "Background Determination for the LUX-ZEPLIN (LZ) Dark
Matter Experiment").

Fig. 12 shows the SR1 accidental-coincidence PDF built by LZ as the outer
product of the smoothed isolated-S1 and isolated-S2 marginal distributions
(text: "The two-dimensional PDF was constructed by taking the outer
product of the smoothed out versions of the S1 and S2 projections"). The
top and right panels plot those two marginal functions (orange curves;
caption: "The orange lines represent the functions that were used to build
the two-dimensional PDF") in events/tonne/yr against S1c [phd] (top) and
log10(S2c) (right). This script extracts both orange curves from the
embedded raster by color segmentation, calibrates all axes from the major
tick marks, and writes data/lz_accidentals_marginals_sr1.csv.

These marginals are the CROSS-CHECK input for the accidentals background
component: the WS2024 model is reconstructed as a factorized
P(S1c)*P(log10 S2c) PDF constrained by the WS2024 Fig. 3 quantile contours
(src/digitize_accidentals_ws2024.py), and its morphology is compared
against these SR1 marginals (LZ states the WS2024 method is identical to
the SR1 one). SR1 configuration: SR1 cuts/ROI (the SR1 raw-S2 600-phd
bound causes the right-panel fall-off at small log10(S2c)).

Runs LOCALLY (not under lc): provenance generation for a data/ input.

Figure facts (established by raster inspection):
- Single embedded 2998x2173 RGB PNG (xref 1218 of page index 20).
- Three panels sharing the S1c x axis (main heatmap rows 599-1922, top
  marginal strip rows 71-554, cols ~305-2376) plus a right marginal panel
  (cols 2441-2926, sharing the main panel's y axis).
- Ticks point INWARD: x majors span rows ~1880-1918 (minors only
  ~1895-1918), anchored under the y majors point right from the left
  spine at cols ~309-345.
- Main/top x: S1c 0-80 phd, majors every 10. Main/right y:
  log10(S2c) 2.75-4.50, majors every 0.25. Top y: log10 rate, majors
  10^-1, 10^-2, 10^-3. Right x: log10 rate, majors 10^-3, 10^-2, 10^-1.
- The orange marginal curves are RGB (255,165,0); orange appears only in
  the two marginal panels (the heatmap itself uses a viridis-like map).
- A gray band covers S1c < 3 (below the ROI bound) inside both left
  panels; the orange S1c marginal therefore starts at the ROI edge.
"""

import csv
import hashlib
import io
from pathlib import Path

import numpy as np
from PIL import Image

PDF = Path(__file__).resolve().parent.parent / "data" / "lz_backgrounds_paper.pdf"
OUT = Path(__file__).resolve().parent.parent / "data" / "lz_accidentals_marginals_sr1.csv"
DIAGNOSTIC = Path("/tmp/lz_sr1_fig12_digitization.png")

SHA256 = "52223054fcdcb65c14787ec4c82c5a2685d222cfd189aa29bdbc230829c7c2e7"
PAGE = 20  # 0-based; printed page 21


def load_raster(doc):
    """The embedded 2998x2173 Fig. 12 PNG as an RGB int array."""
    page = doc[PAGE]
    (xref,) = [im[0] for im in page.get_images(full=True)]
    raw = doc.extract_image(xref)
    return np.asarray(Image.open(io.BytesIO(raw["image"]))).astype(int)


def runs(mask):
    """Contiguous True runs of a 1-D mask as (start, end) pairs."""
    out, start = [], None
    for i, v in enumerate(mask):
        if v and start is None:
            start = i
        elif not v and start is not None:
            out.append((start, i - 1))
            start = None
    if start is not None:
        out.append((start, len(mask) - 1))
    return out


def col_groups(mask_cols, min_hits):
    """Connected column groups (tick marks) with per-column hits >= min_hits."""
    hits = mask_cols.sum(axis=0)
    return [0.5 * (a + b) for a, b in runs(hits >= min_hits)]


def row_groups(mask_rows, min_hits):
    hits = mask_rows.sum(axis=1)
    return [0.5 * (a + b) for a, b in runs(hits >= min_hits)]


def fit(centers, values, tol, what):
    f = np.polyfit(centers, values, 1)
    resid = np.abs(np.polyval(f, centers) - values).max()
    assert resid < tol, f"{what}: calibration residual {resid:.4g}"
    return f


def main():
    digest = hashlib.sha256(PDF.read_bytes()).hexdigest()
    assert digest == SHA256, f"PDF sha256 {digest} != pinned {SHA256}"
    import pymupdf

    doc = pymupdf.open(PDF)
    flat = " ".join(doc[PAGE].get_text().replace("-\n", "").split())
    assert "The accidentals PDF normalized" in flat, "Fig. 12 caption not on page"

    img = load_raster(doc)
    dark = (img[..., 0] < 100) & (img[..., 1] < 100) & (img[..., 2] < 100)
    orange = (
        (img[..., 0] > 200) & (img[..., 1] > 110) & (img[..., 1] < 200)
        & (img[..., 2] < 110)
    )

    # panel spines: continuous dark verticals / horizontals
    cov = dark[599:1923, :].any(axis=0)
    spines_y = [0.5 * (a + b) for a, b in runs(cov) if b - a >= 2]
    covr = dark[:, 310:2370].any(axis=1)
    spines_x = [0.5 * (a + b) for a, b in runs(covr) if b - a >= 2]
    print(f"vertical spines at cols {spines_y}; horizontal spines at rows {spines_x}")

    # --- x calibration (S1c, shared by main + top panels) ---
    # majors: full-length ticks just above the bottom spine (rows 1880-1918).
    # Dark-purple heatmap cells in the bottom-right corner also match, so
    # snap to the arithmetic grid spanned by the extreme major ticks.
    xt_all = [0.5 * (a + b) for a, b in runs(dark[1880:1918, :].sum(axis=0) >= 30)
              if 200 < 0.5 * (a + b) < 2400]
    x0, x8 = 306.5, 2374.5  # S1c = 0 and 80 (the first/last major ticks)
    pred = np.linspace(x0, x8, 9)
    xt = [min(xt_all, key=lambda c: abs(c - p)) for p in pred]
    assert max(abs(c - p) for c, p in zip(xt, pred)) < 5, \
        f"x majors off grid: {xt} vs {pred}"
    xf = fit(xt, np.arange(0.0, 81.0, 10.0), 0.3, "x")

    # --- main/right y calibration (log10 S2c) ---
    # majors point right from the left spine (~col 309); the axis limits
    # 4.50 / 2.75 sit ON the top/bottom spines, so the spine rows merge
    # with the end ticks and 8 groups map to the 8 values 4.50..2.75
    yt = [r + 599 for r in row_groups(dark[599:1923, 312:345], 20)]
    assert len(yt) == 8, f"expected 8 y majors, got {len(yt)}: {yt}"
    yf = fit(yt, np.arange(4.50, 2.74, -0.25), 0.01, "y")

    # --- top panel y calibration (log10 rate) ---
    # interior majors 10^-1, 10^-2, 10^-3; the strip's own spines (at rows
    # 71-74 / 551-554) show up as end groups and are dropped
    tyt_raw = [r + 71 for r in row_groups(dark[71:555, 312:345], 20)]
    assert len(tyt_raw) == 5, f"expected 5 groups in top strip, got {tyt_raw}"
    tyt = sorted(tyt_raw[1:-1])
    tyf = fit(tyt, [-1.0, -2.0, -3.0], 0.02, "top y")

    # --- right panel x calibration (log10 rate) ---
    # majors: full-length ticks in the rows band, excluding the panel's own
    # spines (cols 2441-2444 / 2923-2926)
    rxt = [c + 2440 for c in col_groups(dark[1880:1918, 2440:2930], 30)
           if 2455 < c + 2440 < 2910]
    assert len(rxt) == 3, f"expected 3 log majors in right panel, got {rxt}"
    rxf = fit(sorted(rxt), [-3.0, -2.0, -1.0], 0.02, "right x")

    print(f"x: {1 / xf[0]:.2f} px/phd;  y: {abs(1 / yf[0]):.1f} px per dex(S2c)")
    print(f"top y: {abs(1 / tyf[0]):.1f} px/decade;  right x: {abs(1 / rxf[0]):.1f} px/decade")

    # --- extract the orange curves ---
    # top panel (S1c marginal): per pixel column, median row of orange pixels
    top = np.zeros_like(orange)
    top[76:549, 312:2372] = True
    m = orange & top
    cols_px = np.nonzero(m.any(axis=0))[0]
    s1 = np.polyval(xf, cols_px)
    rate = 10.0 ** np.polyval(
        tyf, np.array([np.median(np.nonzero(m[:, c])[0]) for c in cols_px])
    )
    print(f"\ntop (S1c marginal): {len(cols_px)} columns, S1c [{s1.min():.2f}, "
          f"{s1.max():.2f}], rate [{rate.min():.2e}, {rate.max():.2e}]")

    # right panel (log10 S2c marginal): per pixel row, median col of orange
    right = np.zeros_like(orange)
    right[603:1918, 2450:2920] = True
    m = orange & right
    rows_px = np.nonzero(m.any(axis=1))[0]
    l10 = np.polyval(yf, rows_px)
    rrate = 10.0 ** np.polyval(
        rxf, np.array([np.median(np.nonzero(m[r, :])[0]) for r in rows_px])
    )
    print(f"right (log10S2c marginal): {len(rows_px)} rows, log10(S2c) "
          f"[{l10.min():.3f}, {l10.max():.3f}], rate [{rrate.min():.2e}, "
          f"{rrate.max():.2e}]")

    with open(OUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "# SR1 accidentals-PDF marginal functions, digitized from Fig. 12 of",
                "# arXiv:2211.17120 (Background Determination for LZ), page 21.",
                "# LZ builds the accidentals PDF as the outer product of the smoothed",
                "# isolated-S1 and isolated-S2 marginal distributions; the orange",
                "# curves in the top and right panels are those two functions",
                "# (caption: 'The orange lines represent the functions that were",
                "# used to build the two-dimensional PDF'). Extracted from the",
                "# embedded 2998x2173 raster by color segmentation (orange =",
                "# (255,165,0)), axes calibrated from the major tick marks.",
                "# s1c_marginal: events/tonne/yr vs S1c [phd] (top panel, log y).",
                "# s2c_marginal: events/tonne/yr vs log10(S2c) (right panel, log x).",
                "# Values carry the plotted normalization (SR1 exposure, SR1 ROI and",
                "# cuts; note the fall-off at small log10(S2c) from the SR1 raw-S2",
                "# 600-phd bound) — use as SHAPE only; the WS2024 accidentals model",
                "# is normalized to the WS2024 Table I count (2.8 +/- 0.6 events).",
                "# Provenance: src/digitize_accidentals_sr1.py; in-project copy",
                "# data/lz_backgrounds_paper.pdf, sha256 " + SHA256 + ".",
            ]
        )
        writer.writerow(["s1c_phd", "s1c_marginal_evt_t_yr"])
        for x, v in zip(s1, rate):
            writer.writerow([f"{x:.3f}", f"{v:.4e}"])
        writer.writerow(["log10_s2c", "s2c_marginal_evt_t_yr"])
        for x, v in zip(l10, rrate):
            writer.writerow([f"{x:.4f}", f"{v:.4e}"])
    print(f"\nwrote {OUT}")

    # diagnostic: extraction over the raster + the two marginals
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
        axes[0].imshow(img.astype(np.uint8)[::6, ::6])
        axes[0].set_title("Fig. 12 raster (1/6 scale)")
        axes[1].semilogy(s1, rate, "-", color="darkorange", lw=1)
        axes[1].set_xlabel("S1c (phd)")
        axes[1].set_ylabel("events/tonne/yr")
        axes[1].set_title("digitized S1c marginal (top panel)")
        axes[2].semilogx(rrate, l10, "-", color="darkorange", lw=1)
        axes[2].set_xlabel("events/tonne/yr")
        axes[2].set_ylabel("log10(S2c)")
        axes[2].set_title("digitized log10(S2c) marginal (right panel)")
        fig.tight_layout()
        fig.savefig(DIAGNOSTIC, dpi=130)
        print(f"diagnostic plot: {DIAGNOSTIC}")
    except Exception as e:  # noqa: BLE001 - diagnostic is best-effort
        print(f"diagnostic plot failed: {e}")


if __name__ == "__main__":
    main()
