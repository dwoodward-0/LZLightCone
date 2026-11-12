"""Digitize the public LZ WS2024 Fig. 2 staged NR efficiencies.

This is an authored-input utility, not a declared Lightcone output recipe.
The output CSV contains the main-panel orange pre-ROI curve (trigger + SS
reconstruction/data-analysis cuts) and, independently, the black full curve.
The latter is a digitization cross-check against HEPData t4, not a replacement
for its higher-precision samples. The inset and legend paths are excluded.

PROBE: --output must be a writable probe path. After inspection, the CSV is
authored into data/ as lz_ws2024_fig2_efficiency.csv; do not write declared
results manually. --overlay writes a PDF/raster/vector overlay diagnostic to
a separate probe path. No internal/privately supplied LZ documents are read.
"""
import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
import pymupdf


PAGE = 4  # zero-based; printed paper page 5
SHA256 = "d358bffbaed7eaddd431f9b540790673d650f00d0964282a089a15c84b2f9544"
ORANGE = (1.0, 0.6470588, 0.0)
EFF_STAGE = "pre_roi_ss_data_cuts"
MAIN_X0 = 78.792236328125


def close(a, b, tol=2e-5):
    return a is not None and np.allclose(a, b, atol=tol, rtol=0)


def calibrate(page):
    """Fit main-panel axes from all inward major ticks."""
    xticks, yticks = [], []
    for d in page.get_drawings():
        # The major ticks have this exact stroke family. Minor ticks are
        # deliberately excluded rather than assigned invented energies.
        if not close(d.get("fill"), (0, 0, 0)):
            continue
        if abs((d.get("width") or 0.0) - 0.388208) > 0.002:
            continue
        for item in d["items"]:
            if item[0] != "l":
                continue
            a, b = item[1], item[2]
            if abs(a.x - b.x) < 0.01 and 100 < a.x < 290:
                if 204.0 < min(a.y, b.y) and max(a.y, b.y) < 209.0:
                    xticks.append(a.x)
            if (abs(a.y - b.y) < 0.01
                    and abs(min(a.x, b.x) - MAIN_X0) < 0.01
                    and 55 < a.y < 205):
                yticks.append(a.y)
    if len(xticks) != 7 or len(yticks) != 6:
        raise RuntimeError(f"Unexpected main tick counts: {len(xticks)}, {len(yticks)}")
    xf_points = np.array([(0.0, MAIN_X0)]
                         + list(zip(np.arange(10., 80., 10.), sorted(xticks))))
    yf_points = np.array(list(zip(np.arange(0., 1.01, .2),
                                  sorted(yticks, reverse=True))))
    xfit = np.polyfit(xf_points[:, 1], xf_points[:, 0], 1)
    yfit = np.polyfit(yf_points[:, 1], yf_points[:, 0], 1)
    for name, fit, points in (("energy", xfit, xf_points),
                              ("efficiency", yfit, yf_points)):
        residual = np.abs(np.polyval(fit, points[:, 1]) - points[:, 0])
        if max(residual) > 1e-5:
            raise RuntimeError(f"Bad {name} calibration: {max(residual)}")
        print(f"{name} tick max residual {max(residual):.3g}")
    return xfit, yfit


def extract_curve(page, stage, xfit, yfit):
    matches = []
    required = 45 if stage == EFF_STAGE else 98
    for d in page.get_drawings():
        r = d["rect"]
        color = ORANGE if stage == EFF_STAGE else (0, 0, 0)
        if close(d.get("color"), color) and d.get("fill") is None:
            if (r.x1 - r.x0 > 150 and r.y0 < 75 and r.x0 < 80
                    and len(d["items"]) == required):
                matches.append(d)
    if len(matches) != 1:
        raise RuntimeError(f"{stage}: main panel not unique ({len(matches)})")
    drawing = matches[0]
    if len(drawing["items"]) != required:
        raise RuntimeError(f"{stage}: expected {required} line items")
    points = [drawing["items"][0][1]]
    for item in drawing["items"]:
        if item[0] != "l":
            raise RuntimeError("Unexpected non-line item in main curve")
        if abs(item[1].x - points[-1].x) > 1e-3 or abs(item[1].y - points[-1].y) > 1e-3:
            raise RuntimeError("Discontinuous curve path")
        points.append(item[2])
    xy = np.array([[p.x, p.y] for p in points])
    curve = np.column_stack([np.polyval(xfit, xy[:, 0]),
                             np.polyval(yfit, xy[:, 1])])
    if np.any(np.diff(curve[:, 0]) <= 0):
        raise RuntimeError(f"{stage}: energy is not monotonic along the polyline")
    if np.any(curve[:, 1] < -1e-4) or np.any(curve[:, 1] > 1.0001):
        raise RuntimeError(f"{stage}: nonphysical efficiency")
    # Tiny PDF stroke numerical excursions at the endpoints are clipped.
    curve[:, 1] = np.clip(curve[:, 1], 0.0, 1.0)
    if min(curve[:, 0]) < -0.01 or max(curve[:, 0]) > 76.4:
        raise RuntimeError(f"{stage}: unexpected energy range {curve[[0, -1], 0]}")
    return curve, float(drawing["width"])


def read_two_column(path):
    rows = []
    for row in csv.reader(line for line in Path(path).read_text().splitlines()
                          if not line.startswith("#")):
        try:
            rows.append([float(row[0]), float(row[1])])
        except (ValueError, IndexError):
            continue
    return np.array(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--paper", default="data/lz_ws2024_paper.pdf")
    ap.add_argument("--hepdata", default="data/lz_ws2024_efficiency.csv")
    ap.add_argument("--output", required=True)
    ap.add_argument("--overlay")
    args = ap.parse_args()
    paper = Path(args.paper)
    if hashlib.sha256(paper.read_bytes()).hexdigest() != SHA256:
        raise RuntimeError("Unexpected public WS2024 PDF hash")
    page = pymupdf.open(paper)[PAGE]
    caption = " ".join(page.get_text().split())
    for phrase in ("SS reconstruction and analysis cuts",
                   "S2 ROI", "S1 ROI"):
        if phrase not in caption:
            raise RuntimeError(f"Fig. 2 caption lacks {phrase}")
    xfit, yfit = calibrate(page)
    orange, orange_width = extract_curve(page, EFF_STAGE, xfit, yfit)
    black, black_width = extract_curve(page, "total_fig2", xfit, yfit)
    reference = read_two_column(args.hepdata)
    common = ((reference[:, 0] >= black[0, 0])
              & (reference[:, 0] <= black[-1, 0]))
    residual = np.interp(reference[common, 0], black[:, 0], black[:, 1]) - reference[common, 1]
    # The line-width resolution plus the calibration residual is small. The
    # explicitly declared 0.01 absolute floor also includes HEPData rounding
    # and PDF anti-aliasing; it is fixed before inspecting event comparisons.
    stroke_width = max(orange_width, black_width)
    resolution_kev = abs(xfit[0]) * stroke_width
    resolution_eff = abs(yfit[0]) * stroke_width
    tolerance = 0.01
    max_residual = float(np.abs(residual).max())
    rms_residual = float(np.sqrt(np.mean(residual ** 2)))
    print(f"black/HepData common range {black[0, 0]:.3f}..{black[-1, 0]:.3f} keVnr")
    print(f"black/HepData max residual {max_residual:.5f}; RMS {rms_residual:.5f}; "
          f"stroke resolution {resolution_kev:.5f} keV / "
          f"{resolution_eff:.5f} efficiency")
    if max_residual > tolerance:
        raise RuntimeError("Independent Fig. 2 black curve disagrees with HEPData")
    eval_e = np.linspace(max(orange[0, 0], black[0, 0]),
                         min(orange[-1, 0], black[-1, 0]), 1000)
    ordering = np.interp(eval_e, orange[:, 0], orange[:, 1]) - np.interp(eval_e, black[:, 0], black[:, 1])
    if ordering.min() < -tolerance:
        raise RuntimeError(f"Pre-ROI curve below total beyond tolerance: {ordering.min()}")
    print(f"pre-ROI minus total: minimum {ordering.min():.5f}")
    lines = [
        "# LZ WS2024 Fig. 2, public arXiv:2410.17036v3, printed page 5.",
        f"# Source PDF sha256={SHA256}",
        "# Vector curves: main panel only; 45 orange / 98 black line items; inset excluded.",
        "# pre_roi_ss_data_cuts: S2 trigger + SS reconstruction & analysis cuts, before S2/S1 ROI.",
        "# total_fig2: after all listed efficiency stages; independent digitization cross-check only.",
        "# Production total-efficiency input remains HEPData 10.17182/hepdata.155182.v2/t4.",
        f"# Black/HEPData comparison on common support: max abs={max_residual:.6g}, rms={rms_residual:.6g}, tolerance={tolerance}.",
        f"# Tick-calibrated stroke resolution: {resolution_kev:.6g} keVnr, {resolution_eff:.6g} efficiency; linear interpolation thereafter.",
        "# No published orange uncertainty band/covariance: digitization/interpolation uncertainty remains.",
        "# Energies outside the sampled plotted support are held at endpoint values by the analysis loader.",
        "stage,energy_kevnr,efficiency",
    ]
    with Path(args.output).open("w") as f:
        f.write("\n".join(lines) + "\n")
        writer = csv.writer(f, lineterminator="\n")
        for stage, curve in ((EFF_STAGE, orange), ("total_fig2", black)):
            for energy, eff in curve:
                writer.writerow([stage, f"{energy:.6f}", f"{eff:.7f}"])
    print(f"wrote {args.output}")
    if args.overlay:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        image = page.get_pixmap(matrix=pymupdf.Matrix(2, 2), alpha=False)
        im = np.frombuffer(image.samples, dtype=np.uint8).reshape(
            image.height, image.width, image.n)
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.imshow(im, extent=[0, page.rect.width, page.rect.height, 0])
        for curve, color, label in ((orange, "magenta", "orange extracted"),
                                    (black, "cyan", "black extracted")):
            raw = np.column_stack([(curve[:, 0] - xfit[1]) / xfit[0],
                                   (curve[:, 1] - yfit[1]) / yfit[0]])
            ax.plot(raw[:, 0], raw[:, 1], lw=0.8, color=color, label=label)
        ax.set_xlim(65, 310); ax.set_ylim(215, 45); ax.legend()
        fig.savefig(args.overlay, dpi=150)
        print(f"wrote {args.overlay}")


if __name__ == "__main__":
    main()
