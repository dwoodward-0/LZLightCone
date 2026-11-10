"""Digitize the LZ WS2022 ER and NR calibration bands (arXiv:2207.03764v4 Fig. 1).

Extracts the vector-graphic band curves (tritium ER median + 10%/90%
quantiles, DD-neutron NR median + 10%/90% quantiles) from Fig. 1 of the
WS2022 paper via page.get_drawings(), calibrates the axes from the tick
labels, and writes data/lz_ws2022_band_reference.csv.

Runs LOCALLY (not under lc): it reads the paper cache, which the lc
sandbox cannot see. The output CSV is a declared input of the band
validations (er_band_validation, nr_band_validation).

Figure facts (established by prior inspection of this page):
- Fig. 1 is on page index 3, in the clip Rect(40, 370, 380, 612).
- Each band is three stroked paths in the band's color: the median
  (solid, dashes '[] 0') and the two quantiles (dotted,
  '[ 3.47652 1.50336 ] 0'). ~164 points per path.
- ER paths are blue (0.1176, 0.5647, 1.0); NR paths are red
  (0.8627, 0.0784, 0.2353).
"""

import csv
from pathlib import Path

import fitz  # pymupdf
import numpy as np

PDF = Path(
    "/Users/dwoodward1/.cache/astra/papers/"
    "10.48550_arXiv.2207.03764_v4/paper.pdf"
)
OUT = Path(__file__).resolve().parent.parent / "data" / "lz_ws2022_band_reference.csv"
DIAGNOSTIC = Path("/tmp/ws2022_fig1_digitization.png")

PAGE = 3
CLIP = fitz.Rect(40, 370, 380, 612)
ER_COLOR = (0.11764705926179886, 0.5647059082984924, 1.0)
NR_COLOR = (0.8627451062202454, 0.0784313753247261, 0.23529411852359772)
MEDIAN_DASHES = "[] 0"
DOTTED_DASHES = "[ 3.47652 1.50336 ] 0"
COLOR_TOL = 0.02


def close(c, ref):
    return c is not None and len(c) >= 3 and all(
        abs(a - b) < COLOR_TOL for a, b in zip(c[:3], ref)
    )


def extract_paths(page, color):
    """Return {dashes: [array of PDF-space points, one per path]} for a color."""
    groups = {}
    for d in page.get_drawings():
        if d["type"] not in ("s", "fs"):
            continue
        if not close(d.get("color"), color):
            continue
        pts = []
        for item in d["items"]:
            if item[0] != "l":
                raise ValueError(f"unexpected drawing item {item[0]!r}")
            pts.append((item[1].x, item[1].y))
            pts.append((item[2].x, item[2].y))
        pts = np.unique(np.round(np.array(pts), 4), axis=0)
        inside = (
            (pts[:, 0] > CLIP.x0)
            & (pts[:, 0] < CLIP.x1)
            & (pts[:, 1] > CLIP.y0)
            & (pts[:, 1] < CLIP.y1)
        )
        pts = pts[inside]
        pts = pts[np.argsort(pts[:, 0])]
        groups.setdefault(d.get("dashes"), []).append(pts)
    return groups


def axis_calibration(page):
    """Fit the axes from the tick-label word bboxes around the axes box."""
    words = page.get_text("words")
    AX_X0, AX_X1, AX_Y0, AX_Y1 = 82.0, 291.9, 428.9, 577.1
    xticks, yticks = [], []
    for w in words:
        x0, y0, x1, y1, text = w[:5]
        if not (CLIP.x0 < x0 and x1 < CLIP.x1 and CLIP.y0 < y0 and y1 < CLIP.y1):
            continue
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        try:
            v = float(text)
        except ValueError:
            continue
        if y1 > AX_Y1 and AX_X0 <= cx < AX_X1 + 15:  # below the box
            xticks.append((cx, v))
        if x1 < AX_X0 and AX_Y0 - 15 < cy < AX_Y1 + 15:  # left of the box
            yticks.append((cy, v))
    if len(xticks) < 3 or len(yticks) < 3:
        raise ValueError(f"too few ticks: {len(xticks)} x, {len(yticks)} y")
    xticks.sort()
    yticks.sort()
    xf = np.polyfit([t[0] for t in xticks], [t[1] for t in xticks], 1)
    yf = np.polyfit([t[0] for t in yticks], [t[1] for t in yticks], 1)
    return xf, yf, xticks, yticks


def to_data(pts, xf, yf):
    return np.column_stack(
        [np.polyval(xf, pts[:, 0]), np.polyval(yf, pts[:, 1])]
    )


def main():
    doc = fitz.open(PDF)
    page = doc[PAGE]
    xf, yf, xticks, yticks = axis_calibration(page)
    print("x tick labels:", [(round(c, 1), v) for c, v in xticks])
    print("y tick labels:", [(round(c, 1), v) for c, v in yticks])
    print(f"x fit: value = {xf[0]:.6f}*px + {xf[1]:.4f}")
    print(f"y fit: value = {yf[0]:.6f}*px + {yf[1]:.4f}")
    # sanity vs the axes-box mapping established by prior inspection:
    # x=82 -> ~0, x=291.9 -> ~80 ; y=428.9 -> ~4.50, y=577.1 -> ~2.75
    assert abs(np.polyval(xf, 82.0)) < 0.5
    assert abs(np.polyval(xf, 291.9) - 80) < 1.0
    assert abs(np.polyval(yf, 428.9) - 4.50) < 0.1
    assert abs(np.polyval(yf, 577.1) - 2.75) < 0.1

    bands = {}
    for name, color in (("er", ER_COLOR), ("nr", NR_COLOR)):
        groups = extract_paths(page, color)
        print(
            f"\n{name}: "
            + str({k: [len(p) for p in v] for k, v in groups.items()})
        )
        if set(groups) != {MEDIAN_DASHES, DOTTED_DASHES}:
            raise ValueError(f"{name}: unexpected dash patterns {set(groups)}")
        if len(groups[MEDIAN_DASHES]) != 1 or len(groups[DOTTED_DASHES]) != 2:
            raise ValueError(f"{name}: expected 1 median + 2 dotted paths")
        med = to_data(groups[MEDIAN_DASHES][0], xf, yf)
        dotted = [to_data(p, xf, yf) for p in groups[DOTTED_DASHES]]
        # classify the two quantiles by offset from the median
        roles = {}
        for d in dotted:
            offset = d[:, 1] - np.interp(d[:, 0], med[:, 0], med[:, 1])
            roles["q90" if np.mean(offset) > 0 else "q10"] = d
        if set(roles) != {"q10", "q90"}:
            raise ValueError(f"{name}: quantile classification failed {set(roles)}")
        for label, d in [("median", med), *roles.items()]:
            print(
                f"  {label:6s}: x [{d[0, 0]:.1f}, {d[-1, 0]:.1f}] phd, "
                f"y [{d[:, 1].min():.3f}, {d[:, 1].max():.3f}]"
            )
        bands[name] = {"median": med, **roles}

    # resample onto a common S1c grid: 1 phd from 3 to the narrowest
    # x-extent over all six curves (per column: NaN where undefined)
    grid = np.arange(3.0, 81.0, 1.0)
    cols = {}
    for name in ("er", "nr"):
        xmin = max(b[:, 0].min() for b in bands[name].values())
        xmax = min(b[:, 0].max() for b in bands[name].values())
        span = grid[(grid >= xmin) & (grid <= xmax)]
        for role in ("median", "q10", "q90"):
            vals = np.full_like(grid, np.nan)
            vals[np.isin(grid, span)] = np.interp(
                span, bands[name][role][:, 0], bands[name][role][:, 1]
            )
            cols[f"{name}_{role}"] = vals

    width_ok = True
    for name in ("er", "nr"):
        good = np.isfinite(cols[f"{name}_q10"]) & np.isfinite(cols[f"{name}_q90"])
        hw = 0.5 * (cols[f"{name}_q90"][good] - cols[f"{name}_q10"][good])
        print(
            f"{name}: defined on {good.sum()} grid points, "
            f"half-width median {np.median(hw):.3f} in log10(S2c)"
        )
        width_ok &= np.all(hw > 0)
    if not width_ok:
        raise ValueError("a q10 > q90 inversion occurred")

    with open(OUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "# LZ WS2022 ER and NR calibration bands, digitized from Fig. 1 of",
                "# arXiv:2207.03764v4 (LZ WS2022 first result), page 4.",
                "# Extracted programmatically from the vector graphics (page.get_drawings,",
                "# pymupdf) -- not eyeball-digitized; axes calibrated from the tick labels.",
                "# Columns: s1c_phd, then median and 10%/90% quantiles of log10(S2c) for",
                "# the tritium (er_*) and DD-neutron (nr_*) calibration bands; NaN where",
                "# the curve is not drawn.",
                "# Provenance: src/digitize_ws2022_bands.py; cache entry",
                "# 10.48550/arXiv.2207.03764 v4. Used as the reference for the",
                "# er_band_validation and nr_band_validation outputs.",
            ]
        )
        writer.writerow(
            ["s1c_phd", "er_median", "er_q10", "er_q90",
             "nr_median", "nr_q10", "nr_q90"]
        )
        for i, s1 in enumerate(grid):
            writer.writerow(
                [f"{s1:.0f}"]
                + [
                    "" if not np.isfinite(cols[c][i]) else f"{cols[c][i]:.4f}"
                    for c in (
                        "er_median", "er_q10", "er_q90",
                        "nr_median", "nr_q10", "nr_q90",
                    )
                ]
            )
    print(f"\nwrote {OUT}")

    # diagnostic render: extracted curves over the page raster
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        pix = page.get_pixmap(clip=CLIP, dpi=200)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        axes[0].imshow(np.frombuffer(pix.samples, dtype=np.uint8).reshape(
            pix.height, pix.width, pix.n)[:, :, :3], extent=(
            CLIP.x0, CLIP.x1, CLIP.y1, CLIP.y0))
        for name, color in (("er", ER_COLOR), ("nr", NR_COLOR)):
            for role in ("median", "q10", "q90"):
                b = bands[name][role]
                px = (b[:, 0] - xf[1]) / xf[0]
                py = (b[:, 1] - yf[1]) / yf[0]
                axes[0].plot(px, py, "-", color=color, lw=0.8)
        axes[0].set_title("extracted curves over Fig. 1 raster (PDF coords)")
        for name in ("er", "nr"):
            for role in ("median", "q10", "q90"):
                axes[1].plot(grid, cols[f"{name}_{role}"],
                             label=f"{name} {role}")
        axes[1].set_xlabel("S1c (phd)")
        axes[1].set_ylabel("log10(S2c)")
        axes[1].legend(fontsize=8)
        axes[1].set_title("digitized bands (data coords)")
        fig.tight_layout()
        fig.savefig(DIAGNOSTIC, dpi=130)
        print(f"diagnostic plot: {DIAGNOSTIC}")
    except Exception as e:  # noqa: BLE001 - diagnostic is best-effort
        print(f"diagnostic plot failed: {e}")


if __name__ == "__main__":
    main()
