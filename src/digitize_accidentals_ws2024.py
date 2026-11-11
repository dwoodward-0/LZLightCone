"""Digitize the WS2024 accidental-coincidence quantile contours (Fig. 3 of
arXiv:2410.17036v3).

Fig. 3 (page 7, left column, two stacked panels) overlays the 68% and 95%
quantile contours of the best-fit accidental-coincidence distribution in
(S1c, log10 S2c) space (orange contours; caption: "Green and orange
contours show the same quantiles for the best-fit 124Xe distribution and
the accidental coincidence events, respectively"). Both contours are
vector-drawn. This script extracts the TOP panel's pair ("All Data" — the
panel matching our pooled exposure under exposure_structure.combined_single;
the bottom panel is radon-tagged only), calibrates the axes, and writes
data/lz_ws2024_accidentals_contours.csv.

The contours are clipped at the ROI boundaries (S1c = 3, the S2c lower
bound): the accidentals PDF lives inside the ROI, so contour segments
running along those edges are expected and carry information.

Runs LOCALLY (not under lc): provenance generation for a data/ input, the
same pattern as src/digitize_neutron_spectrum.py.

Figure facts (established by prior inspection):
- Fig. 3 is on page index 7 (printed page 7). Top panel plot area is
  x ~ 90-296, y ~ 52-205 (PDF pts); x tick labels (0-80 phd, shared with
  the bottom panel) sit below the bottom panel at y ~ 367-377.
- The two accidentals contours are the only stroked paths in orange
  (1.0, 0.7568, 0.2314) inside the top panel: the larger-extent one is the
  95% quantile, the compact low-S1c one the 68% quantile (68% nested
  inside 95%; visually confirmed against the page raster).
- Y axis: log10(S2c), labels 2.75-4.25 in 0.25 steps beside the top panel.
"""

import csv
import hashlib
import re
from pathlib import Path

import fitz  # pymupdf
import numpy as np

PDF = Path(__file__).resolve().parent.parent / "data" / "lz_ws2024_paper.pdf"
OUT = Path(__file__).resolve().parent.parent / "data" / "lz_ws2024_accidentals_contours.csv"
DIAGNOSTIC = Path("/tmp/lz_ws2024_fig3_digitization.png")

SHA256 = "d358bffbaed7eaddd431f9b540790673d650f00d0964282a089a15c84b2f9544"
PAGE = 6  # 0-based; printed page 7
ORANGE = (1.0, 0.7568627450980392, 0.23137254901960785)
PANEL = fitz.Rect(85, 45, 300, 206)  # top panel + margins, PDF pts
COLOR_TOL = 0.01


def close(c, ref, tol=COLOR_TOL):
    return c is not None and len(c) >= 3 and all(
        abs(a - b) < tol for a, b in zip(c[:3], ref)
    )


def calibrate(page):
    """x from the shared S1c tick labels, y from the top-panel log10(S2c) labels."""
    words = page.get_text("words")
    xticks, yticks = [], []
    for x0, y0, x1, y1, text, *_ in words:
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        if 360 < y0 < 380 and PANEL.x0 - 5 < cx < PANEL.x1 + 5:  # below bottom panel
            try:
                v = float(text)
            except ValueError:
                continue
            if v in (0, 10, 20, 30, 40, 50, 60, 70, 80):
                xticks.append((cx, v))
        if PANEL.y0 - 5 < cy < PANEL.y1 + 8 and x1 < 85:  # beside top panel
            m = re.fullmatch(r"[234]\.\d\d", text)
            if m:
                yticks.append((cy, float(text)))
    assert len(xticks) == 9, f"expected 9 x labels, got {len(xticks)}"
    assert len(yticks) == 7, f"expected 7 y labels, got {len(yticks)}"
    xf = np.polyfit([t[0] for t in xticks], [t[1] for t in xticks], 1)
    yf = np.polyfit([t[0] for t in yticks], [t[1] for t in yticks], 1)
    xr = np.abs(np.polyval(xf, [t[0] for t in xticks]) - [t[1] for t in xticks]).max()
    yr = np.abs(np.polyval(yf, [t[0] for t in yticks]) - [t[1] for t in yticks]).max()
    assert xr < 0.3 and yr < 0.01, f"calibration residuals x {xr:.3f}, y {yr:.4f}"
    return xf, yf


def extract_contours(page, xf, yf):
    """The two orange stroked paths in the top panel -> data-space point lists."""
    paths = []
    for d in page.get_drawings():
        if d["type"] != "s" or not close(d.get("color"), ORANGE):
            continue
        r = d["rect"]
        if not (PANEL.x0 <= r.x0 and r.x1 <= PANEL.x1 and
                PANEL.y0 <= r.y0 and r.y1 <= PANEL.y1):
            continue
        pts = []
        for item in d["items"]:
            if item[0] == "l":
                pts.append((item[1].x, item[1].y))
                pts.append((item[2].x, item[2].y))
            elif item[0] == "c":
                pts.append((item[1].x, item[1].y))
                pts.append((item[4].x, item[4].y))
            else:
                raise ValueError(f"unexpected item {item[0]!r}")
        pts = np.array(pts)
        data = np.column_stack([np.polyval(xf, pts[:, 0]), np.polyval(yf, pts[:, 1])])
        # collapse consecutive duplicates
        keep = np.ones(len(data), dtype=bool)
        keep[1:] = np.any(np.abs(np.diff(data, axis=0)) > 1e-9, axis=1)
        paths.append(data[keep])
    assert len(paths) == 2, f"expected 2 orange contours in top panel, got {len(paths)}"
    # classify: the 95% contour has the larger S1c extent
    paths.sort(key=lambda p: p[:, 0].max() - p[:, 0].min())
    return {"q68": paths[0], "q95": paths[1]}


def main():
    digest = hashlib.sha256(PDF.read_bytes()).hexdigest()
    assert digest == SHA256, f"PDF sha256 {digest} != pinned {SHA256}"
    doc = fitz.open(PDF)
    page = doc[PAGE]
    flat = " ".join(page.get_text().replace("-\n", "").split())
    assert "the accidental coincidence events" in flat, \
        "Fig. 3 caption not on this page"

    xf, yf = calibrate(page)
    print(f"x fit: S1c = {xf[0]:.5f}*pt + {xf[1]:.3f}")
    print(f"y fit: log10(S2c) = {yf[0]:.6f}*pt + {yf[1]:.4f}")
    contours = extract_contours(page, xf, yf)
    for name, c in contours.items():
        print(f"\n{name}: {len(c)} points, S1c [{c[:, 0].min():.2f}, "
              f"{c[:, 0].max():.2f}] phd, log10(S2c) [{c[:, 1].min():.3f}, "
              f"{c[:, 1].max():.3f}]")
        # ROI clipping: contours must live inside the ROI
        assert c[:, 0].min() >= 2.9, f"{name}: S1c below ROI bound"
        assert c[:, 1].min() >= 2.68, f"{name}: log10(S2c) below ROI bound"
    # nesting: the q68 contour's points should lie inside the q95 hull extent
    q95, q68 = contours["q95"], contours["q68"]
    for label, c in (("q68", q68), ("q95", q95)):
        print(f"{label} S1c extent {c[:, 0].max() - c[:, 0].min():.1f} phd")
    assert q68[:, 0].max() <= q95[:, 0].max() + 0.5, "q68 wider than q95 in S1c"

    with open(OUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "# WS2024 accidental-coincidence 68% and 95% quantile contours in",
                "# (S1c, log10 S2c) space, digitized from Fig. 3 (top panel, 'All",
                "# Data') of arXiv:2410.17036v3, page 7. Caption: 'Green and orange",
                "# contours show the same quantiles for the best-fit 124Xe",
                "# distribution and the accidental coincidence events, respectively.'",
                "# Extracted programmatically from the vector graphics (pymupdf",
                "# page.get_drawings); axes calibrated from the tick labels (S1c",
                "# 0-80 phd shared axis; log10(S2c) 2.75-4.25). Contours are clipped",
                "# at the ROI boundaries (S1c = 3 phd, the S2c lower bound): the",
                "# accidentals PDF lives inside the ROI.",
                "# Columns: contour, s1c_phd, log10_s2c (points in path order).",
                "# SHAPE constraint for the accidentals background component of",
                "# background_components.all_table1 (LZ's construction is a",
                "# factorized P(S1c)*P(log10 S2c) PDF, outer product of the isolated-",
                "# pulse marginals, per arXiv:2211.17120 Sec. V F); normalization",
                "# from the WS2024 Table I count (2.8 +/- 0.6 events,",
                "# lz_table1_counts_reference).",
                "# Provenance: src/digitize_accidentals_ws2024.py; in-project copy",
                "# data/lz_ws2024_paper.pdf, sha256 " + SHA256 + ".",
            ]
        )
        writer.writerow(["contour", "s1c_phd", "log10_s2c"])
        for name in ("q68", "q95"):
            for s1, l10 in contours[name]:
                writer.writerow([name, f"{s1:.3f}", f"{l10:.4f}"])
    print(f"\nwrote {OUT}")

    # diagnostic render: contours over the page raster + data coords
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        clip = fitz.Rect(50, 45, 300, 210)
        pix = page.get_pixmap(clip=clip, dpi=220)
        img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(
            pix.height, pix.width, pix.n)[:, :, :3]
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        axes[0].imshow(img, extent=(clip.x0, clip.x1, clip.y1, clip.y0))
        for name, col in (("q68", "lime"), ("q95", "red")):
            c = contours[name]
            px = (c[:, 0] - xf[1]) / xf[0]
            py = (c[:, 1] - yf[1]) / yf[0]
            axes[0].plot(px, py, "-", color=col, lw=0.6, label=name)
        axes[0].legend(fontsize=8)
        axes[0].set_title("extracted contours over Fig. 3 raster (PDF coords)")
        for name, col in (("q68", "lime"), ("q95", "red")):
            c = contours[name]
            axes[1].plot(c[:, 0], c[:, 1], "-", color=col, lw=1.0, label=name)
        axes[1].set_xlabel("S1c (phd)")
        axes[1].set_ylabel("log10(S2c)")
        axes[1].set_xlim(0, 80)
        axes[1].set_ylim(2.6, 4.6)
        axes[1].axvline(3, color="gray", lw=0.5, ls="--")
        axes[1].legend(fontsize=8)
        axes[1].set_title("digitized accidentals quantile contours (data coords)")
        fig.tight_layout()
        fig.savefig(DIAGNOSTIC, dpi=130)
        print(f"diagnostic plot: {DIAGNOSTIC}")
    except Exception as e:  # noqa: BLE001 - diagnostic is best-effort
        print(f"diagnostic plot failed: {e}")


if __name__ == "__main__":
    main()
