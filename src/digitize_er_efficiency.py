"""Digitize LZ's public WS2024 ER trigger and cumulative S1-stage curves.

Authored-input utility, not a declared Lightcone output. Read only the pinned
public arXiv:2511.17350v1 PDF. Fig. 1 is on zero-based page 5 (printed page
6); use the lower WS2024 panel. The blue curve is the trigger. The yellow
curve is cumulative trigger + threefold coincidence + S1c > 3 phd.

The user verified visually that the green cumulative SS/data-analysis curve
and the black total-1D curve overlap through 15 keV. That observation permits
use of the separate HEPData t3 table as G(E) through 15 keV. This utility
therefore does NOT digitize the green/black/red curves, shade uncertainty
bands, or treat the black high-energy 1D-ROI decline as the green plateau.

--output and --overlay are writable probe paths. After inspection the CSV is
authored into data/. Do not write declared results manually.
"""
import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
import pymupdf

PAGE = 5
SHA256 = "9c4b10d7c93dcbb4323d49c9a892b6b21e9e4c3d7ba55d14eeba6064330d2738"
BLUE = (0.0, 0.0, 1.0)
YELLOW = (1.0, 0.6470588, 0.0)
X0 = 78.12345886230469


def close(a, b, tol=2e-5):
    return a is not None and np.allclose(a, b, atol=tol, rtol=0)


def calibrate(page):
    """Calibrate the lower panel from all nine x and six y major ticks."""
    xticks, yticks = [], []
    for d in page.get_drawings():
        if not close(d.get("fill"), (0, 0, 0)):
            continue
        if abs((d.get("width") or 0.0) - 0.377736) > 0.001:
            continue
        for item in d["items"]:
            if item[0] != "l":
                continue
            a, b = item[1], item[2]
            if (abs(a.x - b.x) < 0.001 and 70 < a.x < 300
                    and 408 < a.y < 415 and abs(b.y - 408.52795) < 0.01):
                xticks.append(a.x)
            if (abs(a.y - b.y) < 0.001 and abs(a.x - X0) < 0.001
                    and 3.7 < b.x - a.x < 3.9 and 250 < a.y < 408):
                yticks.append(a.y)
    if len(xticks) != 9 or len(yticks) != 6:
        raise RuntimeError(f"Unexpected lower-panel major ticks: {len(xticks)}, {len(yticks)}")
    xp = np.column_stack([np.arange(0., 20.1, 2.5), sorted(xticks)])
    yp = np.column_stack([np.arange(0., 1.01, .2), sorted(yticks, reverse=True)])
    xfit = np.polyfit(xp[:, 1], xp[:, 0], 1)
    yfit = np.polyfit(yp[:, 1], yp[:, 0], 1)
    for name, fit, points in (("energy", xfit, xp), ("efficiency", yfit, yp)):
        residual = np.abs(np.polyval(fit, points[:, 1]) - points[:, 0]).max()
        if residual > 1e-5:
            raise RuntimeError(f"Unexpected {name} tick residual {residual}")
        print(f"{name} major-tick max residual: {residual:.3g}")
    return xfit, yfit


def extract(page, name, color, nlines, xfit, yfit):
    matches = []
    for d in page.get_drawings():
        r = d["rect"]
        # The main polyline begins at E=0 and covers the entire panel. The
        # horizontal legend samples have only two line items and are excluded.
        if (close(d.get("color"), color) and d.get("fill") is None
                and r.width > 150 and r.x0 < 79 and 255 < r.y0 < 265
                and len(d["items"]) == nlines
                and abs(d.get("width") - 0.708255) < 0.002):
            matches.append(d)
    if len(matches) != 1:
        raise RuntimeError(f"{name}: main lower-panel curve not unique ({len(matches)})")
    items = matches[0]["items"]
    points = [items[0][1]]
    for item in items:
        if item[0] != "l":
            raise RuntimeError(f"{name}: unexpected non-line path item")
        if abs(item[1].x - points[-1].x) > 1e-3 or abs(item[1].y - points[-1].y) > 1e-3:
            raise RuntimeError(f"{name}: discontinuous path")
        points.append(item[2])
    xy = np.array([[p.x, p.y] for p in points])
    curve = np.column_stack([np.polyval(xfit, xy[:, 0]),
                             np.polyval(yfit, xy[:, 1])])
    # The path starts within 3e-5 keV of the calibrated zero. Preserve that
    # origin instead of letting a negative floating-point roundoff discard it
    # and cause a held first nonzero efficiency below the first sampled vertex.
    curve[np.abs(curve[:, 0]) < abs(xfit[0]) * 0.708255, 0] = 0.0
    if np.any(np.diff(curve[:, 0]) <= 0):
        raise RuntimeError(f"{name}: energy not strictly increasing")
    # PDF clipping retains the last segment outside the 20-keV plot edge.
    if curve[0, 0] < -0.001 or curve[-1, 0] < 20 or curve[-1, 0] > 22:
        raise RuntimeError(f"{name}: unexpected support {curve[[0, -1], 0]}")
    # Interpolate the last retained segment to the exact axis edge instead of
    # dropping it and leaving an unintended extrapolation gap at 20 keV.
    edge_efficiency = float(np.interp(20.0, curve[:, 0], curve[:, 1]))
    curve = curve[(curve[:, 0] >= 0) & (curve[:, 0] <= 20)]
    if curve[-1, 0] < 20:
        curve = np.vstack([curve, [20.0, edge_efficiency]])
    if np.any(curve[:, 1] < -1e-4) or np.any(curve[:, 1] > 1.0001):
        raise RuntimeError(f"{name}: nonphysical efficiency")
    curve[:, 1] = np.clip(curve[:, 1], 0, 1)
    if curve[-1, 0] < 19.99:
        raise RuntimeError(f"{name}: does not reach the 20-keV axis edge")
    if not np.all(np.diff(curve[:, 1]) >= -0.002):
        raise RuntimeError(f"{name}: blue/yellow cumulative curve unexpectedly declines")
    return curve


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--paper", default="data/lz_er_efficiency_paper.pdf")
    ap.add_argument("--output", required=True)
    ap.add_argument("--overlay")
    args = ap.parse_args()
    paper = Path(args.paper)
    digest = hashlib.sha256(paper.read_bytes()).hexdigest()
    if digest != SHA256:
        raise RuntimeError(f"Unexpected public ER-paper PDF hash: {digest}")
    page = pymupdf.open(paper)[PAGE]
    caption = " ".join(page.get_text().split())
    for phrase in ("WS2024", "The trigger (blue)",
                   "three-fold coincidence and > 3 phd on S1c (yellow)",
                   "single scatter reconstruction and data analysis cuts (green)",
                   "1D ROI (black)", "2D ROI (red)"):
        if phrase not in caption:
            raise RuntimeError(f"Fig. 1 caption lacks expected phrase: {phrase}")
    xfit, yfit = calibrate(page)
    blue = extract(page, "trigger", BLUE, 3, xfit, yfit)
    yellow = extract(page, "trigger_coincidence_s1", YELLOW, 13, xfit, yfit)
    grid = np.linspace(0, 20, 2001)
    ordering = np.interp(grid, yellow[:, 0], yellow[:, 1]) - np.interp(grid, blue[:, 0], blue[:, 1])
    stroke = 0.708255
    tolerance = max(0.002, abs(yfit[0]) * stroke * 2)
    if ordering.max() > tolerance:
        raise RuntimeError(f"Cumulative yellow above blue beyond digitization tolerance: {ordering.max()}")
    print(f"yellow minus blue: minimum {ordering.min():.7f}, maximum {ordering.max():.7f} "
          f"(later stage required <= earlier stage within {tolerance:.7f})")
    print(f"yellow values at cumulative-stage vertices: {yellow[:, 1]}")
    lines = [
        "# LZ public arXiv:2511.17350v1, Fig. 1 lower WS2024 panel, printed page 6.",
        f"# Source PDF sha256={SHA256}",
        "# Main-panel vector paths only: trigger (blue, 3 lines); cumulative trigger + threefold coincidence + S1c > 3 phd (yellow, 13 lines).",
        "# X ticks 0..20 keV in 2.5-keV steps; y ticks 0..1 in 0.2 steps; last segment clipped to 20 keV.",
        "# Green is NOT extracted. The user visually confirmed that green cumulative SS/data-analysis and black total-1D curves overlap through 15 keV.",
        "# HEPData 10.17182/hepdata.181350.v1/t3 is the G(E) central-value input on that overlap range; G(E>15) is held at G(15).",
        "# Blue/yellow are independent cumulative stages, not factors to multiply together. Interpolation is linear in energy and efficiency.",
        "# No numerical central-path uncertainty or covariance is recovered here; shaded bands in Fig. 1 are not digitized.",
        f"# Central-path stroke resolution: {abs(xfit[0])*stroke:.6g} keV, {abs(yfit[0])*stroke:.6g} efficiency.",
        "# Energies outside sampled plotted support are held at endpoint values by the analysis loader; reported per production.",
        "stage,energy_kevee,efficiency",
    ]
    with Path(args.output).open("w") as f:
        f.write("\n".join(lines) + "\n")
        writer = csv.writer(f, lineterminator="\n")
        for name, curve in (("trigger", blue), ("trigger_coincidence_s1", yellow)):
            for energy, efficiency in curve:
                writer.writerow([name, f"{energy:.6f}", f"{efficiency:.7f}"])
    print(f"wrote {args.output}")
    if args.overlay:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        image = page.get_pixmap(matrix=pymupdf.Matrix(3, 3), alpha=False)
        im = np.frombuffer(image.samples, dtype=np.uint8).reshape(image.height, image.width, image.n)
        fig, ax = plt.subplots(figsize=(11, 7))
        ax.imshow(im, extent=[0, page.rect.width, page.rect.height, 0])
        for curve, color, label in ((blue, "magenta", "trigger extracted"),
                                    (yellow, "cyan", "cumulative S1 stage extracted")):
            raw = np.column_stack([(curve[:, 0] - xfit[1]) / xfit[0],
                                   (curve[:, 1] - yfit[1]) / yfit[0]])
            ax.plot(raw[:, 0], raw[:, 1], lw=0.8, color=color, label=label)
        ax.set_xlim(68, 310)
        ax.set_ylim(420, 245)
        ax.legend()
        fig.savefig(args.overlay, dpi=160)
        print(f"wrote {args.overlay}")


if __name__ == "__main__":
    main()
