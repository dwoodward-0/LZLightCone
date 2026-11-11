"""Digitize the WS2024 Fig. 4 band-median projection (arXiv:2410.17036v3).

Fig. 4 (page 7, right column, two stacked panels) compares WS2024 data to
the best-fit model "in a projection along the ER band median, normalized
to the 1-sigma band width": within each of two S1c ranges (upper panel
"3 < S1c < 40 phd", lower panel "40 < S1c < 80 phd", annotated in-plot),
events are histogrammed in the band-projected coordinate
u = (log10(S2c) - mu_ER(S1c)) / sigma_ER(S1c). The y axis is events per
bin, log scale 10^-2..10^2. Series (from the legend swatches):

    Data            black markers
    Total Bkg       blue   (0, 0, 1)
    214Pb betas     light blue (0, 0.7608, 0.9765)
    Other ER betas  pink   (1, 0.3529, 0.6863)
    124Xe DEC       teal   (0, 0.6235, 0.5059)
    127Xe+125Xe EC  orange-red (1, 0.4314, 0.2275)
    Accidentals     orange (1, 0.7647, 0.2314)
    8B              light pink (1, 0.698, 0.9922)
    40 GeV/c2 WIMP  purple (0.5176, 0, 0.8039)  [dashed]

All series are vector-drawn step outlines; each horizontal run of a path
is one bin's level. Bins are 0.5556 wide in u (18 bins spanning the
plotted -6..4.5 axis; the exact grid is recovered from the data-marker
centers and confirmed by the step edges). Levels drawn below the panel
floor (y-axis minimum 10^-2) are clipped in the figure and are not
digitized.

Digitization note (established during extraction, verified by exact
ratios): the drawn per-category outlines are mutually consistent with
Total Bkg = 1.628 x 214Pb (= 1210/743, the full Table I ratio) only if
"Other ER betas" as DRAWN contains 85Kr/39Ar/det-gammas + 136Xe alone
(sum ratio 0.293 = 217.6/743); the solar-nu, 212Pb and tritium/14C
components (223 events, 18.4% of the total) are absent from the outlines
while included in the Total Bkg curve. The CSV stores what the figure
shows; the interpretation is documented in the consuming validation.

Runs LOCALLY (not under lc): provenance generation for a data/ input,
the same pattern as src/digitize_accidentals_ws2024.py.
"""

import csv
import hashlib
from pathlib import Path

import fitz  # pymupdf
import numpy as np

PDF = Path(__file__).resolve().parent.parent / "data" / "lz_ws2024_paper.pdf"
OUT = Path(__file__).resolve().parent.parent / "data" / "lz_ws2024_fig4_projection.csv"
DIAGNOSTIC = Path("/tmp/lz_ws2024_fig4_digitization.png")

SHA256 = "d358bffbaed7eaddd431f9b540790673d650f00d0964282a089a15c84b2f9544"
PAGE = 6  # 0-based; printed page 7

SERIES = {
    "total_bkg": ((0.0, 0.0, 1.0), False),
    "pb214_betas": ((0.0, 0.7608, 0.9765), False),
    "other_er_betas": ((1.0, 0.3529, 0.6863), False),
    "xe124_dec": ((0.0, 0.6235, 0.5059), False),
    "xe127_125_ec": ((1.0, 0.4314, 0.2275), False),
    "accidentals": ((1.0, 0.7647, 0.2314), False),
    "b8": ((1.0, 0.698, 0.9922), False),
    "wimp_40gev": ((0.5176, 0.0, 0.8039), False),
}
COLOR_TOL = 0.005
# panel y bounds and the y of the 10^2 label anchor (filled by calibrate)
PANEL_X = (346.1, 549.5)


def close(c, ref, tol=COLOR_TOL):
    return c is not None and len(c) >= 3 and all(
        abs(a - b) < tol for a, b in zip(c[:3], ref)
    )


def log10_labels(page, y_lo, y_hi):
    """The five 10^k y-axis labels of one panel: (union-bbox center, k)."""
    out = []
    words = page.get_text("words")
    for x0, y0, x1, y1, text, *_ in words:
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        if not (y_lo < cy < y_hi) or x1 > 346 or cx > 348:
            continue
        if text != "10":
            continue
        # find the superscript digit word just above-right. The exponent
        # minus (for 10^-1 / 10^-2) is not extracted as a word, but it
        # pushes the digit ~3 pt right of the base's center: positive
        # exponents sit ~5.4 pt right, negative ones ~8.4 pt.
        sup = None
        for x0b, y0b, x1b, y1b, tb, *_ in words:
            cbx, cby = 0.5 * (x0b + x1b), 0.5 * (y0b + y1b)
            if tb in "0123456789" and 330 < cbx < 348 and \
                    -8 < (cby - cy) < -1 and abs(cbx - cx) < 12:
                neg = (cbx - cx) > 7
                sup = ("-" + tb if neg else tb, y0b, y1b)
        if sup is None:
            continue
        top = min(y0, sup[1])
        bot = max(y1, sup[2])
        out.append((0.5 * (top + bot), int(sup[0])))
    out.sort(key=lambda t: t[0])
    assert len(out) == 5 and [k for _, k in out] == [2, 1, 0, -1, -2], \
        f"expected 10^2..10^-2 labels, got {out}"
    return out


def calibrate(page):
    """x (u) from the bottom-panel tick labels; y per panel from 10^k labels.

    Returns per panel (yf) plus the shared x fit (uf), where yf maps PDF y
    to log10(events) and uf maps PDF x to u.
    """
    # x: the u = -6..4 tick marks at the bottom of the bottom panel
    # (thin black marks just above y = 392.0). The tick labels are not
    # usable directly: the minus glyphs are not extracted as words, and
    # the negative labels' digit centers sit ~2.15 pt right of their
    # ticks; the labels are used only as an order cross-check below.
    xticks = []
    for d in page.get_drawings():
        r = d["rect"]
        c = d.get("color") or d.get("fill")
        if c is None or tuple(round(v, 3) for v in c[:3]) != (0.0, 0.0, 0.0):
            continue
        if d["type"] not in ("f", "fs"):
            continue
        if r.width > 0.5 or not (3.0 < r.height < 5.0):
            continue
        if not (388.0 < r.y0 < 391.0 and PANEL_X[0] - 2 < r.x0
                and r.x1 < PANEL_X[1] + 2):
            continue
        xticks.append((0.5 * (r.x0 + r.x1), 0.0))
    xticks = sorted(set(round(x, 2) for x, _ in xticks))
    assert len(xticks) == 6, f"expected 6 x ticks, got {len(xticks)}"
    xticks = [(x, v) for x, v in zip(xticks, (-6, -4, -2, 0, 2, 4))]
    uf = np.polyfit([t[0] for t in xticks], [t[1] for t in xticks], 1)
    xr = np.abs(np.polyval(uf, [t[0] for t in xticks])
                - [t[1] for t in xticks]).max()
    assert xr < 0.15, f"x calibration residual {xr:.3f}"
    # cross-check: the three non-negative digit labels (u = 0, 2, 4; the
    # bare digits of negative labels also pass this filter, so match by
    # position against the fitted line instead of counting)
    labs = sorted(
        0.5 * (w[0] + w[2]) for w in page.get_text("words")
        if 393 < 0.5 * (w[1] + w[3]) < 402
        and PANEL_X[0] - 6 < 0.5 * (w[0] + w[2]) < PANEL_X[1] + 6
        and w[4] in ("0", "2", "4"))
    lab_u = sorted(np.polyval(uf, x) for x in labs)
    # the three non-negative labels must land on integers >= 0
    nonneg = [u for u in lab_u if abs(u - round(u)) < 0.05 and u >= -0.05]
    assert len(nonneg) == 3, f"non-negative labels off integer u: {lab_u}"
    assert len(xticks) == 6, f"expected 6 x tick labels, got {len(xticks)}"
    xticks.sort(key=lambda t: t[0])
    uf = np.polyfit([t[0] for t in xticks], [t[1] for t in xticks], 1)
    xr = np.abs(np.polyval(uf, [t[0] for t in xticks])
                - [t[1] for t in xticks]).max()
    assert xr < 0.15, f"x calibration residual {xr:.3f}"

    yf = {}
    for panel, (y_lo, y_hi) in (("top", (86, 241)), ("bottom", (241, 394))):
        labs = log10_labels(page, y_lo, y_hi)
        ys = [t[0] for t in labs]
        ks = [t[1] for t in labs]
        fit = np.polyfit(ys, ks, 1)
        res = np.abs(np.polyval(fit, ys) - np.array(ks)).max()
        assert res < 0.02, f"{panel} y calibration residual {res:.4f}"
        yf[panel] = fit
    return uf, yf


def panel_of(rect):
    return "top" if rect.y0 < 245 else "bottom"


def extract_runs(page, uf, yf):
    """Per panel and series: horizontal runs (u_lo, u_hi, log10 level) of
    the step outlines, plus the black data-marker centers."""
    runs = {p: {s: [] for s in SERIES} for p in ("top", "bottom")}
    verts = {p: [] for p in ("top", "bottom")}  # u of vertical step edges
    for d in page.get_drawings():
        r = d["rect"]
        c = d.get("color") or d.get("fill")
        if r.x0 < PANEL_X[0] - 2 or r.y1 < 90:
            continue  # left column (Fig. 3) or legend swatches
        match = None
        for name, (ref, _) in SERIES.items():
            if close(c, ref):
                match = name
        panel = panel_of(r)
        ylo, yhi = ((86.8, 239.4) if panel == "top" else (239.4, 392.0))
        if match is not None:
            for it in d["items"]:
                if it[0] != "l":
                    continue
                (xa, ya), (xb, yb) = (it[1].x, it[1].y), (it[2].x, it[2].y)
                if abs(ya - yb) <= 0.5:
                    if not (ylo - 0.5 <= ya <= yhi + 0.5):
                        continue  # clipped below the floor / above the top
                    u1, u2 = np.polyval(uf, xa), np.polyval(uf, xb)
                    if abs(0.5 * (u1 + u2)) > 5.3 or abs(u2 - u1) < 0.02:
                        continue
                    runs[panel][match].append(
                        (min(u1, u2), max(u1, u2),
                         float(np.polyval(yf[panel], ya))))
                elif abs(xa - xb) <= 0.5 and ylo - 0.5 <= min(ya, yb):
                    u = np.polyval(uf, xa)
                    if abs(u) <= 5.3:
                        verts[panel].append(u)
    # data markers
    marks = {"top": [], "bottom": []}
    for d in page.get_drawings():
        r = d["rect"]
        c = d.get("color") or d.get("fill")
        if c is None or tuple(round(v, 3) for v in c[:3]) != (0.0, 0.0, 0.0):
            continue
        if d["type"] not in ("f", "fs"):
            continue
        if not (PANEL_X[0] + 1 < r.x0 and r.x1 < PANEL_X[1] - 1):
            continue
        w, h = r.width, r.height
        if not (1.5 < w < 4.5 and abs(w - h) < 0.2):
            continue
        panel = panel_of(r)
        ylo, yhi = ((86.8, 239.4) if panel == "top" else (239.4, 392.0))
        if not (ylo < r.y0 and r.y1 < yhi):
            continue
        u = np.polyval(uf, 0.5 * (r.x0 + r.x1))
        l10 = np.polyval(yf[panel], 0.5 * (r.y0 + r.y1))
        marks[panel].append((u, l10))
    for panel in marks:
        marks[panel].sort()
    return runs, verts, marks


def extract_data_markers(page, uf, yf):
    """Black circular data markers (w == h ~ 2.9 pt) inside the panels."""
    marks = {"top": [], "bottom": []}
    for d in page.get_drawings():
        r = d["rect"]
        c = d.get("color") or d.get("fill")
        if c is None or tuple(round(v, 3) for v in c[:3]) != (0.0, 0.0, 0.0):
            continue
        if d["type"] not in ("f", "fs"):
            continue
        if not (PANEL_X[0] + 1 < r.x0 and r.x1 < PANEL_X[1] - 1):
            continue
        w, h = r.width, r.height
        if not (1.5 < w < 4.5 and abs(w - h) < 0.2):
            continue
        panel = panel_of(r)
        ylo, yhi = ((86.8, 239.4) if panel == "top" else (239.4, 392.0))
        if not (ylo < r.y0 and r.y1 < yhi):
            continue
        u = np.polyval(uf, 0.5 * (r.x0 + r.x1))
        l10 = np.polyval(yf[panel], 0.5 * (r.y0 + r.y1))
        marks[panel].append((u, l10))
    for panel in marks:
        marks[panel].sort()
    return marks


def main():
    digest = hashlib.sha256(PDF.read_bytes()).hexdigest()
    assert digest == SHA256, f"PDF sha256 {digest} != pinned {SHA256}"
    doc = fitz.open(PDF)
    page = doc[PAGE]
    flat = " ".join(page.get_text().replace("-\n", "").split())
    assert "projection along the ER band median" in flat, \
        "Fig. 4 caption not on this page"
    for ann in ("3 < S1c < 40", "40 < S1c < 80"):
        assert ann.replace(" ", "") in flat.replace(" ", ""), \
            f"panel annotation {ann!r} not found"

    uf, yf = calibrate(page)
    print(f"u fit: u = {uf[0]:.6f}*pt + {uf[1]:.4f}")
    for panel, fit in yf.items():
        print(f"{panel} y fit: log10(events) = {fit[0]:.6f}*pt + {fit[1]:.4f}"
              f"  ({-1/fit[0]:.3f} pt/decade)")

    runs, verts, marks = extract_runs(page, uf, yf)
    for panel in marks:
        assert len(marks[panel]) in (12, 13), \
            f"unexpected marker count {panel}: {len(marks[panel])}"

    # bin grid from the vertical step edges: cluster, take the modal
    # spacing as du, anchor the phase on the densest edge
    all_verts = sorted(round(u, 3) for u in verts["top"] + verts["bottom"])
    clusters = []
    for u in all_verts:
        if clusters and abs(u - clusters[-1][-1]) < 0.05:
            clusters[-1].append(u)
        else:
            clusters.append([u])
    edges = [float(np.mean(c)) for c in clusters]
    dus = np.diff(edges)
    du = float(np.median(dus))
    assert abs(du - 0.5556) < 0.01, f"unexpected bin width {du:.4f}"
    # snap all edges to the global phase (the densest cluster anchors it)
    e0 = edges[len(edges) // 2]
    grid = sorted({round(e0 + k * du, 3)
                   for e in edges for k in [round((e - e0) / du)]
                   if abs(e - (e0 + k * du)) < 0.06} | set(edges))
    print(f"bin grid: {len(grid)} edges, width {du:.4f}, "
          f"u [{grid[0]:.3f}, {grid[-1]:.3f}]")

    # levels per (panel, series, bin) by overlap-length weighting
    final = {}
    for panel in ("top", "bottom"):
        for name in SERIES:
            for u1, u2, l10 in runs[panel][name]:
                for e_lo, e_hi in zip(grid[:-1], grid[1:]):
                    ov = min(u2, e_hi) - max(u1, e_lo)
                    if ov > 0.02:
                        final.setdefault((panel, name, round(e_lo, 3)),
                                         []).append((ov, l10))
    levels = {k: float(np.average([v for _, v in vs], weights=[w for w, _ in vs]))
              for k, vs in final.items()}

    # data markers must sit at bin centers (grid holds the EDGES)
    for panel in ("top", "bottom"):
        for u, _ in marks[panel]:
            k = (u - grid[0]) / du - 0.5
            assert abs(k - round(k)) < 0.1, \
                f"marker at u={u:.3f} off bin centers"

    n_rows = 0
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "# WS2024 Fig. 4 band-median projection, digitized from page 7",
            "# (right column) of arXiv:2410.17036v3. Within each panel's S1c",
            "# range (top: 3 < S1c < 40 phd, bottom: 40 < S1c < 80 phd),",
            "# events are histogrammed in u = (log10(S2c) - mu_ER(S1c)) /",
            "# sigma_ER(S1c); y is events per bin, log scale. Bins are",
            "# 0.5556 wide in u (grid from the step-path edges; the black",
            "# data markers sit at the bin centers). Model series are",
            "# step-outline levels per bin; 'data' rows are the markers.",
            "# Levels the figure clips below its 10^-2 floor are absent.",
            "# Digitization note: the drawn outlines satisfy Total = 1.628",
            "# x 214Pb (= 1210/743, the full Table I ratio) only if",
            "# 'other_er_betas' as drawn holds 85Kr/39Ar/gammas + 136Xe",
            "# alone (ratio 0.293 = 217.6/743); solar-nu, 212Pb and",
            "# tritium/14C (223 events, 18.4%) are in Total Bkg but not in",
            "# the drawn outlines. WIMP shown at its power-constrained 90%",
            "# CL upper limit (caption).",
            "# Columns: panel, series, u_center, u_lo, u_hi, log10_events.",
            "# Provenance: src/digitize_fig4_projection.py; in-project copy",
            "# data/lz_ws2024_paper.pdf, sha256 " + SHA256 + ".",
        ])
        w.writerow(["panel", "series", "u_center", "u_lo", "u_hi",
                    "log10_events"])
        for panel in ("top", "bottom"):
            for name in list(SERIES) + ["data"]:
                if name == "data":
                    for u, l10 in marks[panel]:
                        w.writerow([panel, name, f"{u:.4f}",
                                    f"{u - du/2:.4f}", f"{u + du/2:.4f}",
                                    f"{l10:.4f}"])
                        n_rows += 1
                    continue
                for (p, nm, e_lo), l10 in sorted(levels.items(),
                                                 key=lambda kv: (kv[0][2],)):
                    if p != panel or nm != name:
                        continue
                    uc = e_lo + du / 2
                    w.writerow([panel, name, f"{uc:.4f}", f"{e_lo:.4f}",
                                f"{e_lo + du:.4f}", f"{l10:.4f}"])
                    n_rows += 1
    print(f"wrote {OUT} ({n_rows} rows)")

    # sanity: reconstructed panel sums of the total_bkg curve
    for panel in ("top", "bottom"):
        tot = sum(10 ** v for (p, nm, e), v in levels.items()
                  if p == panel and nm == "total_bkg")
        pb = sum(10 ** v for (p, nm, e), v in levels.items()
                 if p == panel and nm == "pb214_betas")
        print(f"{panel}: total_bkg sum {tot:.1f} events, "
              f"214Pb sum {pb:.1f}, ratio {tot/pb:.3f}")

    # diagnostic render
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
        colors = {**{n: ref for n, (ref, _) in SERIES.items()},
                  "data": "black"}
        for ax, panel in zip(axes, ("top", "bottom")):
            for name in list(SERIES) + ["data"]:
                if name == "data":
                    pts = marks[panel]
                    ax.plot([p[0] for p in pts], [p[1] for p in pts],
                            "o", ms=3, color="black", label=name)
                    continue
                xs = sorted(k[2] for k in levels if k[0] == panel
                            and k[1] == name)
                if not xs:
                    continue
                ys = [levels[(panel, name, e)] for e in xs]
                ax.step([e + du / 2 for e in xs], ys,
                        where="mid", color=colors[name], lw=1.0, label=name)
            ax.set_xlabel("u = (log10 S2c - mu_ER) / sigma_ER")
            ax.set_title(f"{panel} panel")
            ax.set_ylim(-4, 2.7)
        axes[0].set_ylabel("log10(events / bin)")
        axes[0].legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(DIAGNOSTIC, dpi=130)
        print(f"diagnostic plot: {DIAGNOSTIC}")
    except Exception as e:  # noqa: BLE001 - diagnostic is best-effort
        print(f"diagnostic plot failed: {e}")


if __name__ == "__main__":
    main()
