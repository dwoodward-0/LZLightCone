"""Validation 3: verify and compare the staged public WS2024 NR inputs.

The orange digitized Fig. 2 curve is trigger + SS/data-analysis cuts, before
S2 and S1 ROI. Only that pre-ROI curve precedes NEST in the baseline; the
observable ROI is then applied once. HEPData t4 is the black cumulative total,
which already includes those ROI stages. The independently digitized black
curve is displayed only as a vector/HEPData consistency check.

The full-curve E50 gate remains 5.4 +/- 0.1 keVnr. There is deliberately no
such gate on the pre-ROI curve, and no curve division or multiplication.
"""
import argparse
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

THEMES = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e",
              "grid": "#e4e3e0", "blue": "#2a78d6", "orange": "#eb6834"},
    "dark": {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7",
             "grid": "#383835", "blue": "#3987e5", "orange": "#d95926"},
}
E50_PUBLISHED = 5.4
E50_TOL = 0.1
BLACK_HEPDATA_TOL = 0.01


def read_total(path):
    with open(path) as stream:
        rows = [row for row in csv.reader(stream)
                if row and not row[0].startswith("#")
                and not row[0].lower().startswith("energy")]
    curve = np.array([[float(row[0]), float(row[1])] for row in rows])
    if len(curve) < 2 or np.any(np.diff(curve[:, 0]) <= 0):
        raise ValueError("Total-efficiency energies must be strictly increasing")
    return curve


def read_stage(path, requested):
    with open(path) as stream:
        rows = [row for row in csv.reader(stream)
                if row and not row[0].startswith("#") and row[0] != "stage"]
    curve = np.array([[float(row[1]), float(row[2])]
                      for row in rows if len(row) >= 3 and row[0] == requested])
    if len(curve) < 2 or np.any(np.diff(curve[:, 0]) <= 0):
        raise ValueError(f"No valid staged curve {requested!r}")
    return curve


def curve_comparisons(total, pre, black):
    common = (total[:, 0] >= black[0, 0]) & (total[:, 0] <= black[-1, 0])
    if not np.any(common):
        raise ValueError("No black/HEPData common support")
    residual = np.interp(total[common, 0], black[:, 0], black[:, 1]) - total[common, 1]
    lo, hi = max(pre[0, 0], total[0, 0]), min(pre[-1, 0], total[-1, 0])
    energies = np.linspace(lo, hi, 2000)
    ordering = np.interp(energies, pre[:, 0], pre[:, 1]) - np.interp(
        energies, total[:, 0], total[:, 1])
    if np.max(np.abs(residual)) > BLACK_HEPDATA_TOL:
        raise ValueError("Independent Fig. 2 black curve disagrees with HEPData")
    if np.min(ordering) < -BLACK_HEPDATA_TOL:
        raise ValueError("Pre-ROI curve lies substantially below total efficiency")
    return float(np.max(np.abs(residual))), float(np.sqrt(np.mean(residual**2)))


def e50_crossing(curve):
    for (e0, f0), (e1, f1) in zip(curve[:-1], curve[1:]):
        if f0 < 0.5 <= f1:
            if f1 == f0:
                raise ValueError("Flat E50 interval cannot be interpolated")
            return e0 + (0.5 - f0) / (f1 - f0) * (e1 - e0)
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="HEPData total CSV")
    parser.add_argument("--staged-input", required=True, help="staged Fig. 2 CSV")
    parser.add_argument("--efficiency-treatment", required=True,
                        choices=["input_pre_roi_curve", "input_published_curve"])
    parser.add_argument("--theme", choices=THEMES, default="light")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    total = read_total(args.input)
    pre = read_stage(args.staged_input, "pre_roi_ss_data_cuts")
    black = read_stage(args.staged_input, "total_fig2")
    maximum, rms = curve_comparisons(total, pre, black)
    e50 = e50_crossing(total)
    if e50 is None or abs(e50 - E50_PUBLISHED) > E50_TOL:
        print(f"FAIL: HEPData E50 {e50} outside 5.4 +/- 0.1 keVnr", file=sys.stderr)
        return 1
    print(f"HEPData E50: {e50:.3f} keVnr; black max|residual| {maximum:.6f}, "
          f"RMS {rms:.6f} (tolerance {BLACK_HEPDATA_TOL})")

    t = THEMES[args.theme]
    fig, ax = plt.subplots(figsize=(8.8, 5.5), dpi=160)
    fig.patch.set_facecolor(t["surface"])
    ax.set_facecolor(t["surface"])
    ax.plot(pre[:, 0], pre[:, 1], color=t["orange"], lw=2.2,
            label="Pre-ROI: trigger + SS/data-analysis cuts (Fig. 2 orange)")
    ax.plot(total[:, 0], total[:, 1], color=t["blue"], lw=2.2,
            label="Total: pre-ROI + S2 ROI + S1 ROI (HEPData t4)")
    ax.plot(black[:, 0], black[:, 1], color=t["ink2"], lw=1.2, ls=(0, (3, 3)),
            label="Independent Fig. 2 black extraction (check only)")
    ax.plot([e50], [0.5], marker="o", ms=9, color=t["blue"],
            markeredgecolor=t["surface"], markeredgewidth=2, zorder=4)
    ax.annotate(f"Total E50 = {e50:.2f} keVnr\n(published 5.4 +/- 0.1)",
                xy=(e50, 0.5), xytext=(7.0, 0.32), color=t["ink"], fontsize=9,
                arrowprops=dict(arrowstyle="-", color=t["ink2"], lw=0.8))
    ax.set(xlim=(0, 18), ylim=(0, 1.02), xlabel="Nuclear recoil energy [keVnr]",
           ylabel="Cumulative NR signal efficiency")
    mode = "baseline uses orange before NEST; ROI once afterwards" \
        if args.efficiency_treatment == "input_pre_roi_curve" \
        else "legacy diagnostic: total before NEST + explicit ROI"
    ax.set_title("LZ WS2024 NR efficiency: staged inputs", loc="left",
                 fontsize=13, color=t["ink"], pad=30)
    ax.text(0.0, 1.02, mode, transform=ax.transAxes,
            color=t["ink2"], fontsize=9.5)
    ax.grid(color=t["grid"], lw=0.7, alpha=0.8)
    ax.tick_params(colors=t["ink2"])
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_color(t["ink2"])
    for spine in ax.spines.values():
        spine.set_color(t["grid"])
    for key in ("top", "right"):
        ax.spines[key].set_visible(False)
    legend = ax.legend(loc="lower right", fontsize=8.8, facecolor=t["surface"])
    for text in legend.get_texts():
        text.set_color(t["ink"])
    ax.text(0.0, -0.19, f"Black / HEPData common-support max residual = {maximum:.4f} "
            f"(gate 0.01). Inputs are interpolated linearly, with endpoint holding.",
            transform=ax.transAxes, fontsize=8.5, color=t["ink2"])
    fig.tight_layout()
    fig.savefig(args.output, facecolor=t["surface"], bbox_inches="tight")
    print(f"PASS: wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
