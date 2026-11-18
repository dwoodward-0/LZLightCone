"""Digitize LZ's WS2024-only SI limit (TeVPA results talk, slide 24).

Source: data/lz_ws2024_tevpa_talk.pdf (sha256-pinned input), page index 23,
"WS2024-only Spin Independent Limit": the 3.3 t-y WS2024 exposure alone,
frequentist two-sided PLR, power constrained at -1 sigma. The plot is one
embedded 2048x1479 RGB raster, so extraction is by color segmentation:

  * axes: frame spines and inward major ticks (y majors 1e-45..1e-48,
    x majors 10, 100, 1000 GeV; the frame top/right are 1e-44 / 1e4 GeV),
    least-squares fit of pixel -> log10 value per axis;
  * green (+-1 sigma) and yellow (+-2 sigma) bands: per column, the extent
    of green / green-or-yellow fill. Black curves drawn on top of a band
    edge (the power-constrained observed limit follows the -1 sigma edge)
    hide it, so a dark run touching a band edge is absorbed into the band
    and the edge is taken at the run's center line;
  * median expected limit (thin dashed black): dark runs strictly inside
    the green band, at most MAX_MEDIAN_RUN px tall, cleaned with a running
    median in log space; dash gaps are interpolated;
  * observed power-constrained limit (thick solid black): the lowest
    thick dark run at or below the green band's lower edge.

Output: CSV on the analysis mass grid (log-log interpolation of the
per-column traces) plus a per-column trace CSV and an overlay PNG for
visual verification.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pymupdf

PAGE_INDEX = 23
EXPECTED_SIZE = (2048, 1479)
DARK_SUM = 250          # r+g+b below this is "black" (gray 128 lines excluded)
MAX_MEDIAN_RUN = 6      # dashed median line thickness bound, px
SMOOTH_HALF = 7         # running-median half window, columns
EDGE_PX = 25            # columns this close to the right frame are not trusted (median dash
                        # contaminated by frame ticks up to ~20 px in)
UPPER_EDGE_MIN_GEV = 17 # +1/+2 sigma edges occluded by overlaid curves below this (visual check)
GRID = [9, 11, 13, 16, 17, 19, 21, 23, 26, 29, 32, 36, 40, 43, 46, 65, 91,
        129, 182, 256, 361, 508, 1008, 2000, 5000, 10000]


def load_raster(pdf):
    doc = pymupdf.open(pdf)
    page = doc[PAGE_INDEX]
    imgs = page.get_images(full=True)
    if len(imgs) != 1:
        raise ValueError(f"expected one raster on slide 24, found {len(imgs)}")
    pix = pymupdf.Pixmap(doc, imgs[0][0])
    if (pix.width, pix.height) != EXPECTED_SIZE or pix.n != 3:
        raise ValueError("unexpected raster geometry")
    return np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, 3).astype(int)


def runs(mask):
    """(start, stop_inclusive) of True runs in a 1D mask."""
    d = np.diff(np.concatenate([[0], mask.astype(int), [0]]))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1) - 1))


def calibrate(a):
    dark = a.sum(2) < DARK_SUM
    cs, rs = dark.sum(0), dark.sum(1)
    left, right = sorted(np.argsort(cs)[-6:])[1], sorted(np.argsort(cs)[-6:])[-2]
    top, bottom = sorted(np.argsort(rs)[-6:])[1], sorted(np.argsort(rs)[-6:])[-2]
    # inward tick lengths along the inside of the left/bottom spines
    def tick_len(line):
        k = 0
        while k < 60 and line[k]:
            k += 1
        return k
    lx = left + 2
    ylen = np.array([tick_len(dark[r, lx:]) for r in range(dark.shape[0])])
    by = bottom - 2
    xlen = np.array([tick_len(dark[by::-1, c]) for c in range(dark.shape[1])])
    def majors(lengths, lo, hi):
        idx = np.flatnonzero(lengths[lo:hi] >= 27) + lo
        groups = np.split(idx, np.flatnonzero(np.diff(idx) > 2) + 1)
        return [float(g.mean()) for g in groups if len(g)]
    ym = majors(ylen, top + 20, bottom - 20)
    xm = majors(xlen, left + 20, right - 20)
    if len(ym) != 3 or len(xm) != 3:
        raise ValueError(f"major ticks not found as expected: y {ym}, x {xm}")
    # y majors (rows increasing downward) -> 1e-46, 1e-47, 1e-48; frame top 1e-44
    y_pix = np.array([float(top + 1)] + ym)
    y_val = np.array([-44.0, -46.0, -47.0, -48.0])
    x_pix = np.array(xm + [float(right)])
    x_val = np.array([1.0, 2.0, 3.0, 4.0])
    ycoef = np.polyfit(y_pix, y_val, 1)
    xcoef = np.polyfit(x_pix, x_val, 1)
    yres = np.polyval(ycoef, y_pix) - y_val
    xres = np.polyval(xcoef, x_pix) - x_val
    return {"left": int(left), "right": int(right), "top": int(top), "bottom": int(bottom),
            "y_majors_px": ym, "x_majors_px": xm, "ycoef": ycoef, "xcoef": xcoef,
            "y_fit_max_resid_dex": float(np.abs(yres).max()),
            "x_fit_max_resid_dex": float(np.abs(xres).max())}


def cluster_extent(rows, anchor, direction, max_gap):
    """Extent of a gap-tolerant cluster of rows starting next to anchor.

    direction -1 walks upward (decreasing row) from anchor, +1 downward.
    Returns the far edge row, or None if no row lies within max_gap.
    """
    rows = np.sort(rows[(rows < anchor) if direction < 0 else (rows > anchor)])
    if direction < 0:
        rows = rows[::-1]
    edge, last = None, anchor
    for r in rows:
        if abs(r - last) > max_gap:
            break
        edge = last = r
    return edge


def trace(a, cal):
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    green = (g > 225) & (r > 165) & (r < 232) & (b > 105) & (b < 185)
    yellow = (r > 225) & (g > 225) & (b > 150) & (b < 222)
    dark = a.sum(2) < DARK_SUM
    rows_out = []
    for c in range(cal["left"] + 4, cal["right"] - 3):
        col_d = dark[:, c].copy()
        col_d[: cal["top"] + 3] = False
        col_d[cal["bottom"] - 3:] = False
        gr = np.flatnonzero(green[:, c])
        gr = gr[(gr > cal["top"] + 3) & (gr < cal["bottom"] - 3)]
        if gr.size < 5:
            continue
        # main green cluster: the largest gap-tolerant (<= 12 px) group
        groups = np.split(gr, np.flatnonzero(np.diff(gr) > 12) + 1)
        gmain = max(groups, key=len)
        g_top, g_bot = int(gmain.min()), int(gmain.max())
        drun = runs(col_d)
        top_edge, bot_edge = float(g_top), float(g_bot)
        for s, e in drun:  # black line drawn on a green edge: take its center
            if g_top - 4 <= e < g_top and e - s < 12:
                top_edge = 0.5 * (s + e)
            if g_bot < s <= g_bot + 4 and e - s < 12:
                bot_edge = 0.5 * (s + e)
        yr = np.flatnonzero(yellow[:, c])
        up = cluster_extent(yr, g_top, -1, 14)
        lo = cluster_extent(yr, g_bot, +1, 22)
        if up is None or lo is None:
            continue
        # median: thin dark runs strictly inside green
        inner = [(s, e) for s, e in drun
                 if s > g_top + 2 and e < g_bot - 2 and (e - s + 1) <= MAX_MEDIAN_RUN]
        median_row = 0.5 * (inner[0][0] + inner[0][1]) if len(inner) == 1 else np.nan
        # observed: thick dark run at/just below the lower green edge
        obs = [(s, e) for s, e in drun if g_bot - 2 <= s <= g_bot + 6 and (e - s + 1) >= 4]
        obs_row = 0.5 * (obs[0][0] + obs[0][1]) if obs else np.nan
        rows_out.append([c, float(up), top_edge, median_row, bot_edge, float(lo), obs_row])
    return np.array(rows_out, float)


def running_median(y, half):
    out = np.full_like(y, np.nan)
    for i in range(len(y)):
        w = y[max(0, i - half): i + half + 1]
        w = w[np.isfinite(w)]
        if w.size >= 3:
            out[i] = np.median(w)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--talk", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--trace-output", type=Path)
    ap.add_argument("--overlay", type=Path)
    args = ap.parse_args()
    a = load_raster(args.talk)
    cal = calibrate(a)
    t = trace(a, cal)
    lx = np.polyval(cal["xcoef"], t[:, 0])
    ly = np.polyval(cal["ycoef"], t[:, 1:])
    names = ["plus2sigma", "plus1sigma", "median", "minus1sigma", "minus2sigma", "observed_pc"]
    # clean every traced curve: reject points > 0.03 dex from a wide running median (dash gaps on the
    # median; crossing comparison curves spike the band edges), then interpolate over the rejects.
    # Columns within EDGE_PX of the right frame are dropped (frame/tick contamination) and the curve
    # is extended there by a straight log-log fit over the last decade (all curves are linear there).
    near_edge = t[:, 0] > cal["right"] - EDGE_PX
    coverage = {}
    for j, n in enumerate(names):
        y = ly[:, j].copy()
        y[near_edge] = np.nan
        rm = running_median(y, SMOOTH_HALF)
        y[~(np.abs(y - rm) <= 0.03)] = np.nan
        ok = np.isfinite(y)
        coverage[n] = float(ok.mean())
        filled = np.interp(lx, lx[ok], y[ok])
        tail = ok & (lx > lx[ok].max() - 1.0)
        k = np.polyfit(lx[tail], y[tail], 1)
        beyond = lx > lx[ok].max()
        filled[beyond] = np.polyval(k, lx[beyond])
        ly[:, j] = filled
    grid = np.log10(np.array(GRID, float))
    outside = grid < lx[0] - 1e-9
    if outside.any():
        print(f"note: grid masses {[GRID[i] for i in np.flatnonzero(outside)]} GeV lie below the band "
              f"start {10**lx[0]:.2f} GeV; written as nan (not extrapolated)")
    vals = {n: np.where(outside, np.nan, 10 ** np.interp(grid, lx, ly[:, j])) for j, n in enumerate(names)}
    occluded = grid < np.log10(UPPER_EDGE_MIN_GEV) - 1e-9
    for n in ("plus1sigma", "plus2sigma"):
        vals[n] = np.where(occluded, np.nan, vals[n])
    order = ["mass_gev", "minus2sigma", "minus1sigma", "median", "plus1sigma", "plus2sigma", "observed_pc"]
    lines = ["# LZ WS2024-only SI limit, TeVPA talk slide 24 (data/lz_ws2024_tevpa_talk.pdf, sha256 "
             + hashlib.sha256(args.talk.read_bytes()).hexdigest() + "), raster color segmentation; "
             "cm2 vs GeV/c2; src/digitize_ws2024_only_limit.py",
             ",".join(order)]
    for i, m in enumerate(GRID):
        lines.append(",".join([f"{m:g}"] + [f"{vals[n][i]:.6g}" for n in order[1:]]))
    args.output.write_text("\n".join(lines) + "\n")
    print(f"axes: x {10**lx[0]:.2f}-{10**lx[-1]:.0f} GeV; fit max resid y {cal['y_fit_max_resid_dex']:.4f} "
          f"x {cal['x_fit_max_resid_dex']:.4f} dex; coverage after cleaning "
          + ", ".join(f"{n} {v:.2f}" for n, v in coverage.items()))
    i43 = int(np.argmin(np.abs(lx - np.log10(43))))
    imin = int(np.nanargmin(ly[:, 5]))
    print(f"observed min {10**ly[imin,5]:.3g} cm2 at {10**lx[imin]:.1f} GeV (slide: 2.3e-48 @ 43); "
          f"at 43 GeV {10**ly[i43,5]:.3g}")
    if args.trace_output:
        np.savetxt(args.trace_output, np.column_stack([10 ** lx, 10 ** ly]), delimiter=",",
                   header="mass_gev," + ",".join(names), comments="", fmt="%.6g")
    if args.overlay:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(12, 8.7))
        ax.imshow(a.astype(np.uint8))
        cols = t[:, 0]
        for j, (n, sty) in enumerate(zip(names, ["m-", "r-", "b-", "r-", "m-", "c-"])):
            rowpix = (ly[:, j] - cal["ycoef"][1]) / cal["ycoef"][0]
            ax.plot(cols, rowpix, sty, lw=0.6, label=n)
        gx = (grid - cal["xcoef"][1]) / cal["xcoef"][0]
        ax.plot(gx, (np.log10(vals["median"]) - cal["ycoef"][1]) / cal["ycoef"][0], "bo", ms=3)
        ax.legend(loc="upper right", fontsize=7)
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(args.overlay, dpi=150)


if __name__ == "__main__":
    main()
