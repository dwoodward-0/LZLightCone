"""Validation 3 (efficiency_passthrough): plot the declared LZ WS2024 NR
efficiency curve and check its 50% point against the published value.

Recipe-ready: reads the declared input CSV and writes the figure to exactly
the path handed in via {output}. Fails (nonzero exit) if the 50% crossing
falls outside the acceptance tolerance (5.4 +/- 0.1 keVnr), so the gate in
the output's acceptance criterion is enforced, not just drawn.
"""

import argparse
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3e0"
SERIES_1 = "#2a78d6"

E50_PUBLISHED = 5.4
E50_TOL = 0.1


def read_curve(path):
    energies, effs = [], []
    with open(path) as f:
        for row in csv.reader(f):
            if not row or row[0].startswith("#") or row[0].lower().startswith("energy"):
                continue
            energies.append(float(row[0]))
            effs.append(float(row[1]))
    return energies, effs


def e50_crossing(energies, effs):
    """Linear interpolation of the first upcrossing of 50%."""
    for i in range(1, len(effs)):
        if effs[i - 1] < 0.5 <= effs[i]:
            e0, e1, f0, f1 = energies[i - 1], energies[i], effs[i - 1], effs[i]
            return e0 + (0.5 - f0) / (f1 - f0) * (e1 - e0)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="declared efficiency CSV")
    ap.add_argument("--output", required=True, help="figure path to write")
    args = ap.parse_args()

    energies, effs = read_curve(args.input)
    e50 = e50_crossing(energies, effs)
    if e50 is None:
        print("FAIL: curve never reaches 50% efficiency", file=sys.stderr)
        return 1
    print(f"E50 (interpolated) = {e50:.3f} keVnr; published = {E50_PUBLISHED} keVnr")
    if abs(e50 - E50_PUBLISHED) > E50_TOL:
        print(
            f"FAIL: E50 = {e50:.3f} keVnr outside {E50_PUBLISHED} +/- {E50_TOL} keVnr",
            file=sys.stderr,
        )
        return 1

    fig, ax = plt.subplots(figsize=(7, 4.4), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.plot(energies, effs, color=SERIES_1, lw=2, solid_capstyle="round")

    # 50% reference marker
    ax.plot([e50], [0.5], "o", ms=9, color=SERIES_1, mec=SURFACE, mew=2, zorder=5)
    ax.axhline(0.5, color=INK_2, lw=0.8, ls=(0, (4, 4)), alpha=0.6)
    ax.annotate(
        f"50% at {e50:.1f} keVnr\n(published: {E50_PUBLISHED} keVnr)",
        xy=(e50, 0.5),
        xytext=(e50 + 3.5, 0.42),
        fontsize=9,
        color=INK,
        arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.8),
    )

    ax.set_xlabel("Nuclear recoil energy [keVnr]", fontsize=10, color=INK)
    ax.set_ylabel("Detection efficiency", fontsize=10, color=INK)
    ax.set_title(
        "LZ WS2024 NR data-selection efficiency (declared input, HEPData t4)",
        fontsize=11,
        color=INK,
        loc="left",
        pad=12,
    )
    ax.set_xlim(0, max(energies))
    ax.set_ylim(0, 1.02)
    ax.tick_params(colors=INK_2, labelsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.grid(True, color=GRID, lw=0.7, alpha=0.8)
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(args.output, facecolor=SURFACE)
    print(f"PASS: wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
