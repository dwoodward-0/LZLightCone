"""Total-background review: all our sources versus LZ WS2024 Fig. 4.

Our model is the sum of all 12 stored Table I background sources, with the
124Xe LL charge-yield set to our Fig. 4 RMS-optimal q=0.78 via validation
4's analytic S2 rescaling for S1c > 40 events. This affects only 124Xe; it
is a figure-matching diagnostic rather than a likelihood fit to data. Its
red band is +/-1 sigma from finite-MC statistics only (sum of weights
squared, independently for each source). LZ's solid black step outline is its drawn
best-fit total background. The blue band is digitized directly from the
figure's displayed uncertainty envelope, which includes its systematic and
statistical uncertainties; it is not a Table I normalization rescaling.

The PDF constructs the displayed blue envelope with two fills: a translucent
blue patch supplies its upper edge and a following white mask supplies its
lower edge. This script extracts horizontal runs from both paths. At the low-u
edges the lower edge extends beyond Fig. 4's 10^-2 axis floor; those bins fill
to this review's lower axis limit and are not assigned an invented lower error.

Unlike the older component-reference CSV, this direct extraction retains the
full 18-bin, -6 < u < 4 grid shown in Fig. 4, including its left-edge bins.
The panel definitions and ER projection are identical to validation 4.

Standalone: writes one PNG to --output. The declared Validation 4 figure
(background_projection_validation) is drawn through draw_total() from
src/validate_backgrounds.py, with the six-gate verdict as a footer.
"""
import argparse
import hashlib
from pathlib import Path

import fitz  # pymupdf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import digitize_fig4_projection as dig
from plot_backgrounds_per_source import (
    PANEL_S1C, XE124_Q_LL_BEST_RMS, load_band, load_models, u_hist)

MODEL_COLOR = "#d62728"
LZ_BAND_COLOR = "#0072b2"
PLOT_FLOOR = 2e-3


def horizontal_steps(drawing, panel, uf, yf, centers):
    """Extract finite horizontal levels, keyed by nearest standard bin."""
    levels = {}
    for item in drawing["items"]:
        if item[0] != "l":
            continue
        a, b = item[1], item[2]
        if abs(a.y - b.y) > 0.5 or abs(a.x - b.x) <= 1.0:
            continue
        if a.y > 1000.0:
            # Matplotlib's log-scale clipping sentinel in this figure.
            continue
        u = [float(np.polyval(uf, x)) for x in (a.x, b.x)]
        if min(u) < -6.01 or max(u) > 4.01:
            continue
        k = int(np.argmin(np.abs(centers - np.mean(u))))
        if abs(centers[k] - np.mean(u)) > 0.01:
            raise RuntimeError(f"Horizontal run off bin centers: {u}")
        levels[k] = float(np.polyval(yf[panel], a.y))
    return levels


def select_paths(page, predicate):
    return [(i, d) for i, d in enumerate(page.get_drawings())
            if d["rect"].width > 200.0 and d["rect"].x0 >= 340.0
            and 86.0 <= d["rect"].y0 < 392.0 and predicate(d)]


def extract_lz(paper_path):
    """Return per-panel total level, band lower/upper, and vector bin edges."""
    paper_path = Path(paper_path)
    if hashlib.sha256(paper_path.read_bytes()).hexdigest() != dig.SHA256:
        raise RuntimeError(f"Unexpected WS2024 PDF hash: {paper_path}")
    page = fitz.open(paper_path)[dig.PAGE]
    uf, yf = dig.calibrate(page)
    drawings = page.get_drawings()
    centers = -6.0 + (np.arange(18) + 0.5) * (10.0 / 18.0)

    totals = select_paths(
        page, lambda d: dig.close(d.get("color"), dig.SERIES["total_bkg"][0]))
    bands = select_paths(
        page,
        lambda d: dig.close(d.get("fill"), (0.1176, 0.5647, 1.0), 0.002)
        and abs((d.get("fill_opacity") or 0.0) - 0.3) < 0.01)
    if len(totals) != 2 or len(bands) != 2:
        raise RuntimeError(f"Expected two total and band paths; got "
                           f"{len(totals)}, {len(bands)}")

    out = {}
    for panel in PANEL_S1C:
        total_candidates = [d for _, d in totals if dig.panel_of(d["rect"]) == panel]
        band_candidates = [(i, d) for i, d in bands if dig.panel_of(d["rect"]) == panel]
        if len(total_candidates) != 1 or len(band_candidates) != 1:
            raise RuntimeError(f"Non-unique Fig. 4 paths for {panel}")
        total = horizontal_steps(total_candidates[0], panel, uf, yf, centers)
        i, blue = band_candidates[0]
        upper = horizontal_steps(blue, panel, uf, yf, centers)
        if i + 1 >= len(drawings) or not dig.close(drawings[i + 1].get("fill"),
                                                   (1.0, 1.0, 1.0)):
            raise RuntimeError(f"Missing {panel} white lower-band mask")
        lower = horizontal_steps(drawings[i + 1], panel, uf, yf, centers)
        if set(total) != set(range(18)) or set(upper) != set(range(18)):
            raise RuntimeError(f"Incomplete {panel} total/band digitization")
        if any(lower[k] > upper[k] for k in lower):
            raise RuntimeError(f"Inverted uncertainty envelope in {panel}")
        out[panel] = {
            "total": 10 ** np.array([total[k] for k in range(18)]),
            "lower_log": np.array([lower.get(k, np.nan) for k in range(18)]),
            "upper": 10 ** np.array([upper[k] for k in range(18)]),
            "n_lower_unbounded": sum(k not in lower for k in range(18)),
        }
    return out, centers


def draw_total(models_path, response_path, paper_path, output,
               title_suffix="", footer=None):
    """Draw the all-source total vs Fig. 4 figure; return per-panel sums."""
    models = load_models(models_path)
    sources = sorted(models)
    if len(sources) != 12:
        raise RuntimeError(f"Expected 12 background sources, got {sources}")
    band = load_band(response_path)
    lz, centers = extract_lz(paper_path)
    edges = -6.0 + np.arange(19) * (10.0 / 18.0)

    fig, axes = plt.subplots(2, 1, figsize=(9, 7.6), sharex=True,
                             gridspec_kw={"hspace": 0.12})
    sums = {}
    for ax, panel in zip(axes, PANEL_S1C):
        ref = lz[panel]
        lower = np.where(np.isfinite(ref["lower_log"]),
                         10 ** ref["lower_log"], PLOT_FLOOR)
        ax.fill_between(
            centers, lower, ref["upper"], step="mid", color=LZ_BAND_COLOR,
            alpha=0.20, linewidth=0,
            label="LZ Fig. 4 uncertainty (systematic + stat.)")
        ax.step(centers, ref["total"], where="mid", color="0.15", lw=1.8,
                label="LZ WS2024 Fig. 4 total bkg (best fit)")

        h, var = u_hist(models, sources, band, edges, panel,
                        q_ll=XE124_Q_LL_BEST_RMS)
        nz = h > 0
        sigma_mc = np.sqrt(var[nz])
        ax.fill_between(
            centers[nz], np.maximum(h[nz] - sigma_mc, 0.0),
            h[nz] + sigma_mc, step="mid", color=MODEL_COLOR, alpha=0.12,
            linewidth=0, label="This Work ±1σ MC stat.")
        ax.step(centers[nz], h[nz], where="mid", color=MODEL_COLOR,
                lw=1.8, ls="--", label="This Work")

        ax.set_yscale("log")
        ax.set_ylim(PLOT_FLOOR, 2e3)
        ax.set_ylabel("events / bin")
        ax.grid(alpha=0.25, lw=0.5)
        ax.set_title(
            f"TOTAL BACKGROUND — all 12 sources  —  {panel} panel "
            f"({PANEL_S1C[panel][0]:.0f}<S1c<"
            f"{PANEL_S1C[panel][1]:.0f}){title_suffix}", fontsize=11)
        sums[panel] = {"this_work": float(h.sum()),
                       "lz_drawn": float(ref["total"].sum()),
                       "lz_lower_unbounded_bins": int(ref["n_lower_unbounded"])}
        print(f"{panel}: This Work = {h.sum():.3f}; "
              f"LZ drawn curve = {ref['total'].sum():.3f} events; "
              f"LZ lower band below figure floor in "
              f"{ref['n_lower_unbounded']} bins")
    axes[1].set_xlim(-6.0, 4.0)
    axes[1].set_xlabel("u = (log10 S2c - mu_ER(S1c)) / sigma_ER(S1c)")
    axes[0].legend(fontsize=8.5, loc="upper left")
    if footer:
        fig.text(0.5, -0.01, footer, ha="center", va="top", fontsize=8.5,
                 family="monospace")
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", output)
    return sums


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--models", required=True)
    ap.add_argument("--response", required=True)
    ap.add_argument("--paper", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    draw_total(args.models, args.response, args.paper, args.output)


if __name__ == "__main__":
    main()
