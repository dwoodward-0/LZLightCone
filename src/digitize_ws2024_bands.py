"""Digitize the LZ WS2024 ER and NR calibration bands from the TeVPA talk.

The WS2024 paper (arXiv:2410.17036) does not publish the calibration bands,
but the collaboration's TeVPA results talk does (slide 9, "Calibrations"):
the (S1c, log10 S2c) plane with the tritium ER and DD-neutron NR calibration
scatter and their band overlays (solid median + dashed quantile pair), at the
WS2024 detector parameters (97 V/cm, g1 = 0.112, g2 = 34.0).

Unlike the WS2022 reference (vector paths in the paper PDF, extracted via
page.get_drawings), this plot is a single pasted 2048x1481 raster image, so
the curves are extracted by COLOR SEGMENTATION: the ER band curves are light
blue and the NR band curves red/crimson, distinct from the navy/orange
calibration scatter, the gray iso-energy lines, and the black axes. The
solid median is the largest connected component of its color mask; the dashed
quantiles are the smaller components, assigned above/below the median. Axes
are calibrated from the tick marks (x: 0-80 phd, majors every 10; y: 2.50-
4.50, majors every 0.25), cross-checked against the frame corners.

Quantile convention: the slide does not label the dashed curves. LZ's
WS2022 Fig. 1 used the 10%/90% quantiles; the same convention is assumed
here and confirmed empirically by the band validations (the simulated 10/90
quantiles either match the dashed curves or they do not).

Runs locally (uv run python). Output: data/lz_ws2024_band_reference.csv,
same schema as data/lz_ws2022_band_reference.csv (1-phd S1c grid, median and
10%/90% quantiles of log10(S2c) per band, NaN where the curve is not drawn).
"""

import csv
from collections import deque
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pymupdf
from PIL import Image

PDF = Path(__file__).resolve().parent.parent / "data" / "lz_ws2024_tevpa_talk.pdf"
OUT = Path(__file__).resolve().parent.parent / "data" / "lz_ws2024_band_reference.csv"
DIAGNOSTIC = Path("/tmp/lz_tevpa_talk/ws2024_band_digitization.png")

PAGE = 8  # slide 9
SHA256 = "617bddb94385af704021fb495c676ea81fb6103bcddf51b392843b80e27301f6"

# color masks (JPEG-tolerant thresholds; probed on the embedded image)
def er_curve_mask(rgb):
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    return (b > 180) & (r < 120) & (g > 100) & (g < 210)


def nr_curve_mask(rgb):
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    return (r > 180) & (g < 90) & (b < 120)


def dilate(mask, it=1):
    out = mask.copy()
    for _ in range(it):
        out = (out | np.roll(out, 1, 0) | np.roll(out, -1, 0)
               | np.roll(out, 1, 1) | np.roll(out, -1, 1))
    return out


def erode(mask, it=1):
    return ~dilate(~mask, it)


def connected_components(mask, min_size=40):
    """4-connected components as list of (ys, xs) index arrays."""
    seen = np.zeros_like(mask, dtype=bool)
    comps = []
    ys, xs = np.nonzero(mask)
    for y0, x0 in zip(ys, xs):
        if seen[y0, x0]:
            continue
        comp_y, comp_x = [], []
        q = deque([(y0, x0)])
        seen[y0, x0] = True
        while q:
            y, x = q.popleft()
            comp_y.append(y)
            comp_x.append(x)
            for ny, nx in ((y - 1, x), (y + 1, x), (y, x - 1), (y, x + 1)):
                if (0 <= ny < mask.shape[0] and 0 <= nx < mask.shape[1]
                        and mask[ny, nx] and not seen[ny, nx]):
                    seen[ny, nx] = True
                    q.append((ny, nx))
        if len(comp_y) >= min_size:
            comps.append((np.array(comp_y), np.array(comp_x)))
    return comps


def find_frame(dark):
    """Locate the axes frame from long dark lines; return (x0, x1, y0, y1)."""
    h, w = dark.shape
    col_frac = dark.sum(axis=0) / h
    row_frac = dark.sum(axis=1) / w
    # frame lines span most of the plot; require >60% coverage
    vlines = np.nonzero(col_frac > 0.6)[0]
    hlines = np.nonzero(row_frac > 0.6)[0]
    if len(vlines) < 2 or len(hlines) < 2:
        raise ValueError(f"frame not found: {len(vlines)} v, {len(hlines)} h")
    return vlines[0], vlines[-1], hlines[0], hlines[-1]


def _runs(hits):
    """Centers of contiguous True runs in a 1-D boolean array."""
    centers = []
    i = 0
    while i < len(hits):
        if hits[i]:
            j = i
            while j < len(hits) and hits[j]:
                j += 1
            centers.append(0.5 * (i + j - 1))
            i = j
        else:
            i += 1
    return np.array(centers)


def main():
    doc = pymupdf.open(PDF)
    page = doc[PAGE]
    imgs = page.get_images(full=True)
    if len(imgs) != 1:
        raise ValueError(f"expected 1 image on page {PAGE}, got {len(imgs)}")
    raw = doc.extract_image(imgs[0][0])
    img = np.asarray(Image.open(__import__("io").BytesIO(raw["image"]))).astype(int)
    h, w, _ = img.shape
    print(f"embedded image: {w}x{h} {raw['ext']}")

    dark = (img[..., 0] < 100) & (img[..., 1] < 100) & (img[..., 2] < 100)
    fx0, fx1, fy0, fy1 = find_frame(dark)
    print(f"frame: x [{fx0}, {fx1}], y [{fy0}, {fy1}]")

    # ticks point INWARD; major and minor strokes have identical length, but
    # there are exactly 41 evenly spaced ticks per axis (x: 0-80 phd step 2;
    # y: 2.50-4.50 step 0.05), so calibrate on all of them
    in_band = 14

    # --- x ticks: vertical strokes just above the bottom edge ---
    band = dark[fy1 - in_band: fy1 - 3, fx0 + 2: fx1 - 1]
    xt = _runs(band.any(axis=0)) + fx0 + 2
    if len(xt) != 41:
        raise ValueError(f"expected 41 x ticks, got {len(xt)} at {xt}")
    xvals = np.arange(0.0, 81.0, 2.0)
    xf = np.polyfit(xt, xvals, 1)
    resid = np.polyval(xf, xt) - xvals
    print(f"x ticks: 41 at [{xt[0]:.1f}, {xt[-1]:.1f}], "
          f"fit residuals max {np.abs(resid).max():.3f} phd")

    # --- y ticks: horizontal strokes just right of the left edge ---
    band = dark[fy0 + 2: fy1 - 1, fx0 + 3: fx0 + 3 + in_band]
    yt = _runs(band.any(axis=1)) + fy0 + 2
    if len(yt) != 41:
        raise ValueError(f"expected 41 y ticks, got {len(yt)} at {yt}")
    yvals = np.arange(4.50, 2.49, -0.05)  # top to bottom
    yf = np.polyfit(yt, yvals, 1)
    resid = np.polyval(yf, yt) - yvals
    print(f"y ticks: 41 at [{yt[0]:.1f}, {yt[-1]:.1f}], "
          f"fit residuals max {np.abs(resid).max():.4f} dex")

    # cross-check: frame corners
    print(f"  frame corners -> x: {np.polyval(xf, [fx0, fx1])}, "
          f"y: {np.polyval(yf, [fy0, fy1])}")

    rgb = img
    inside = np.zeros((h, w), dtype=bool)
    inside[fy0 + 2: fy1 - 1, fx0 + 2: fx1 - 1] = True

    curves = {}
    for name, mask_fn in (("er", er_curve_mask), ("nr", nr_curve_mask)):
        mask = mask_fn(rgb) & inside
        # bridge gaps where scatter points overdraw the curves
        closed = erode(dilate(mask, 3), 3)
        comps = connected_components(closed, min_size=40)
        if not comps:
            raise ValueError(f"{name}: no components found")
        # median = the component with the widest x-span
        spans = [c[1].max() - c[1].min() for c in comps]
        med_i = int(np.argmax(spans))
        med_y, med_x = comps[med_i]
        print(f"{name}: {len(comps)} components; median span "
              f"{spans[med_i]} px ({spans[med_i] / (fx1 - fx0) * 80:.0f} phd), "
              f"{sum(len(c[0]) for c in comps)} px total")
        med_yv = np.polyval(yf, med_y)
        med_xv = np.polyval(xf, med_x)
        order = np.argsort(med_xv)
        med = np.column_stack([med_xv[order], med_yv[order]])

        roles = {"q10": [], "q90": []}
        for ci, (cy, cx) in enumerate(comps):
            if ci == med_i:
                continue
            cxv = np.polyval(xf, cx)
            cyv = np.polyval(yf, cy)
            ref = np.interp(cxv, med[:, 0], med[:, 1])
            off = cyv - ref
            key = "q90" if np.mean(off) > 0 else "q10"
            roles[key].append((cxv, cyv))
        tracks = {"median": med}
        for key in ("q10", "q90"):
            if not roles[key]:
                raise ValueError(f"{name}: no {key} components")
            pts = np.concatenate([np.column_stack(p) for p in roles[key]])
            pts = pts[np.argsort(pts[:, 0])]
            tracks[key] = pts
            print(f"  {key}: {len(roles[key])} dash segments, "
                  f"x [{pts[0, 0]:.1f}, {pts[-1, 0]:.1f}] phd, "
                  f"mean offset {np.mean(pts[:, 1] - np.interp(pts[:, 0], med[:, 0], med[:, 1])):+.3f} dex")
        curves[name] = tracks

    # resample onto the shared 1-phd grid
    grid = np.arange(3.0, 81.0, 1.0)
    cols = {}
    for name in ("er", "nr"):
        xmin = max(t[:, 0].min() for t in curves[name].values())
        xmax = min(t[:, 0].max() for t in curves[name].values())
        span = grid[(grid >= xmin) & (grid <= xmax)]
        for role in ("median", "q10", "q90"):
            vals = np.full_like(grid, np.nan)
            t = curves[name][role]
            ok = np.isin(grid, span)
            vals[ok] = np.interp(grid[ok], t[:, 0], t[:, 1])
            cols[f"{name}_{role}"] = vals

    width_ok = True
    print("\nhalf-widths and curve noise:")
    for name in ("er", "nr"):
        good = np.isfinite(cols[f"{name}_q10"]) & np.isfinite(cols[f"{name}_q90"])
        hw = 0.5 * (cols[f"{name}_q90"][good] - cols[f"{name}_q10"][good])
        # pointwise noise: second differences of the median (smooth curve)
        med = cols[f"{name}_median"][good]
        noise = np.abs(np.diff(med, 2)) / np.sqrt(6)
        print(f"  {name}: defined on {good.sum()} pts, half-width median "
              f"{np.median(hw):.3f} min {hw.min():.3f} max {hw.max():.3f}; "
              f"median 2nd-diff noise {np.median(noise):.4f}, "
              f"p90 {np.quantile(noise, 0.9):.4f} dex")
        width_ok &= np.all(hw > 0)
    if not width_ok:
        raise ValueError("a q10 > q90 inversion occurred")

    with open(OUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "# LZ WS2024 ER and NR calibration bands, digitized from slide 9",
            "# ('Calibrations') of the collaboration's TeVPA results talk:",
            f"#   {PDF.name}, sha256 {SHA256}",
            "#   https://indico.uchicago.edu/event/427/contributions/1325/attachments/359/548/lz_results_tevpa.pdf",
            "# The WS2024 paper (arXiv:2410.17036) does not publish the bands; this",
            "# slide shows the tritium ER and DD-neutron NR calibration data with",
            "# band overlays (solid median + dashed quantile pair) at the WS2024",
            "# detector parameters (97 V/cm, g1 = 0.112, g2 = 34.0). The dashed",
            "# quantiles are assumed to be the 10%/90% quantiles (LZ's WS2022 Fig. 1",
            "# convention); the band validations confirm this empirically.",
            "# Extracted programmatically from the embedded 2048x1481 raster by color",
            "# segmentation (ER curves light blue, NR curves red; distinct from the",
            "# navy/orange scatter) -- not eyeball-digitized; axes calibrated from",
            "# the tick marks.",
            "# Provenance: src/digitize_ws2024_bands.py. Used as the reference for",
            "# the er_band_validation_ws2024 and nr_band_validation_ws2024 outputs.",
            "# Columns: s1c_phd, then median and 10%/90% quantiles of log10(S2c) for",
            "# the tritium (er_*) and DD-neutron (nr_*) bands; NaN where the curve",
            "# is not drawn.",
        ])
        writer.writerow([
            "s1c_phd", "er_median", "er_q10", "er_q90",
            "nr_median", "nr_q10", "nr_q90",
        ])
        for i, s1 in enumerate(grid):
            writer.writerow(
                [f"{s1:.0f}"]
                + ["" if not np.isfinite(cols[c][i]) else f"{cols[c][i]:.4f}"
                   for c in ("er_median", "er_q10", "er_q90",
                             "nr_median", "nr_q10", "nr_q90")]
            )
    print(f"\nwrote {OUT}")

    # diagnostic: extracted curves over the raster + data-coords plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].imshow(np.asarray(Image.open(
        __import__("io").BytesIO(raw["image"]))))
    styles = {"er": ("tab:blue", "-"), "nr": ("tab:red", "-")}
    for name, (color, ls) in styles.items():
        for role in ("median", "q10", "q90"):
            t = curves[name][role]
            px = (t[:, 0] - xf[1]) / xf[0]
            py = (t[:, 1] - yf[1]) / yf[0]
            axes[0].plot(px, py, ls, color=color, lw=1.0,
                         alpha=1.0 if role == "median" else 0.6)
    axes[0].set_title("extracted curves over slide 9 raster (px coords)")
    for name in ("er", "nr"):
        for role in ("median", "q10", "q90"):
            axes[1].plot(grid, cols[f"{name}_{role}"],
                         label=f"{name} {role}")
    axes[1].set_xlabel("S1c (phd)")
    axes[1].set_ylabel("log10(S2c)")
    axes[1].legend(fontsize=8)
    axes[1].set_title("digitized WS2024 bands (data coords)")
    fig.tight_layout()
    fig.savefig(DIAGNOSTIC, dpi=130)
    print(f"diagnostic plot: {DIAGNOSTIC}")


if __name__ == "__main__":
    main()
