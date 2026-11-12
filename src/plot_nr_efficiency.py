"""Plot the declared matched NR-efficiency comparison and its MC controls.

All event shapes are unit area inside the explicit observable ROI, including
zero-physical-weight detector neutrons. They are not Table I rate changes.
Efficiency curves have distinct stages: the orange pre-ROI curve is applied
in recoil-energy space before NEST, whereas the black/HEPData total already
contains the S2 and S1 ROI stages and must not precede another explicit ROI.

The main comparison uses complete, matched-budget 12-source productions. The
B8+hep difference panel additionally shows the independent-seed fivefold
study. Difference bands are per-mode multinomial MC uncertainties added in
quadrature; they do not assume cancellation of covariance between the two
rejection samples. No ratio is formed against an empty bin.

A table-view summary is included as text. Light and dark palettes are
separately validated. Legends and line styles make treatment independent of
color, and the thin marks/recessive grids follow the shared chart design.
"""
import argparse
import json

import h5py
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
SOURCES = ("atm_nu_nr", "b8_hep_nu_nr", "detector_neutrons")
LABELS = ("Atmospheric neutrinos", "8B + hep neutrinos", "Detector neutrons")
MODES = {"legacy_total": "Legacy: total efficiency before NEST",
         "corrected_pre_roi": "Corrected: pre-ROI before NEST"}


def add_summary(fig, metrics, axes, t, active):
    fig.suptitle("Correcting the NR selection chain: pre-ROI efficiency, NEST, then one ROI",
                 x=0.065, y=0.995, ha="left", color=t["ink"], fontsize=14)
    subtitle = ("Complete matched-budget 12-source A/B; unit-area shapes conditional on ROI. "
                "Table I physical central rates are unchanged.")
    fig.text(0.065, 0.970, subtitle, ha="left", color=t["ink2"], fontsize=9.5)
    for ax in np.ravel(axes):
        ax.set_facecolor(t["surface"])
        ax.grid(color=t["grid"], lw=0.65, alpha=0.8, zorder=0)
        ax.tick_params(colors=t["ink2"], labelsize=8)
        ax.xaxis.label.set_color(t["ink2"])
        ax.yaxis.label.set_color(t["ink2"])
        ax.title.set_color(t["ink"])
        for spine in ax.spines.values():
            spine.set_color(t["grid"])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax = axes[0, 0]
    total = np.asarray(archive["efficiency_inputs/total_hepdata"])
    pre = np.asarray(archive["efficiency_inputs/pre_roi_fig2_orange"])
    black = np.asarray(archive["efficiency_inputs/independent_total_fig2_black"])
    ax.plot(pre[:, 0], pre[:, 1], color=t["orange"], lw=2,
            label="Pre-ROI trigger + SS/data cuts")
    ax.plot(total[:, 0], total[:, 1], color=t["blue"], lw=2,
            label="Total including S2 + S1 ROI")
    ax.plot(black[:, 0], black[:, 1], color=t["ink2"], lw=1.0, ls="--",
            label="Independent black check")
    ax.set(xlim=(0, 15), ylim=(0, 1.02), xlabel="Recoil energy [keVnr]",
           ylabel="Efficiency", title="1. Do not apply the total efficiency twice")
    residual = metrics["controls"]["black_hepdata_max_abs_efficiency_residual"]
    ax.text(0.04, 0.72, f"Black / HEPData max residual {residual:.4f}\n"
            "Pre-ROI replaces total; curves are not multiplied.",
            transform=ax.transAxes, fontsize=8.5, color=t["ink2"])
    legend(ax, t, fontsize=7.7, loc="lower right")

    x = np.arange(len(SOURCES))
    for column, (key, title) in enumerate((
            ("pre_efficiency_acceptance", "2. Selection before NEST"),
            ("roi_conditional_acceptance", "3. ROI after NEST, conditional on selection"))):
        ax = axes[0, column + 1]
        for pos, mode, color in ((0, "legacy_total", t["blue"]),
                                 (1, "corrected_pre_roi", t["orange"])):
            vals = [metrics["source_metrics"][mode][source][key] for source in SOURCES]
            bars = ax.bar(x + (pos - 0.5) * 0.3, vals, width=0.27, color=color,
                          label=MODES[mode], zorder=2)
            ax.bar_label(bars, fmt="%.3f", fontsize=7.4, color=t["ink2"],
                         padding=2)
        ax.set(xticks=x, xticklabels=("Atm. nu", "8B + hep", "Det. n"),
               ylim=(0, 1.18), ylabel="Acceptance", title=title)
        legend(ax, t, fontsize=7.1, loc="upper left" if column == 0 else "upper right")
    axes[0, 1].text(0.5, -0.22, "Corrected selection keeps more low-energy recoils;\n"
                    "the subsequent observable ROI is therefore more selective.",
                    transform=axes[0, 1].transAxes, ha="center", fontsize=8.3,
                    color=t["ink2"])
    axes[0, 2].text(0.5, -0.22, "Both stages use one common NEST response and\n"
                    "the identical explicit corrected-observable ROI.",
                    transform=axes[0, 2].transAxes, ha="center", fontsize=8.3,
                    color=t["ink2"])


def legend(ax, t, **kwargs):
    handle = ax.legend(facecolor=t["surface"], edgecolor=t["grid"], **kwargs)
    for text in handle.get_texts():
        text.set_color(t["ink"])
    return handle


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="nr_efficiency_comparison HDF5")
    ap.add_argument("--efficiency-treatment", required=True,
                    choices=["input_pre_roi_curve", "input_published_curve"])
    ap.add_argument("--theme", choices=THEMES, default="light")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    t = THEMES[args.theme]

    global archive
    with h5py.File(args.input) as archive:
        metrics = json.loads(archive["summary"].asstr()[()])
        if archive.attrs["active_efficiency_treatment"] != args.efficiency_treatment:
            raise ValueError("Archive and requested efficiency treatment disagree")
        if not metrics["controls"]["active_model_matches_declared_full_model"]:
            raise ValueError("Active model control did not pass")
        if not metrics["validation4"]["corrected_pre_roi"]["acceptance"]["pass"]:
            raise ValueError("Corrected mode does not pass all six Validation 4 gates")

        fig, axes = plt.subplots(3, 3, figsize=(15.7, 13.0), dpi=150,
                                 gridspec_kw={"height_ratios": [1.14, 1, 0.95]})
        fig.patch.set_facecolor(t["surface"])
        add_summary(fig, metrics, axes, t, args.efficiency_treatment)
        colors = {"legacy_total": t["blue"], "corrected_pre_roi": t["orange"]}
        for column, (source, label) in enumerate(zip(SOURCES, LABELS)):
            ax, difference = axes[1, column], axes[2, column]
            for mode in MODES:
                group = archive[f"shapes/{source}/u/top"]
                edges = group["edges"][:]
                centers = 0.5 * (edges[:-1] + edges[1:])
                value_key = ({"legacy_total": "legacy_unit_mass_per_bin",
                              "corrected_pre_roi": "corrected_unit_mass_per_bin"})[mode]
                values = group[value_key][:]
                n = int(archive[f"{mode}/{source}/hist_mc_counts"][:].sum())
                sigma = np.sqrt(np.maximum(values * (1 - values), 0) / n)
                ax.step(centers, values, where="mid", lw=2.0, color=colors[mode],
                        label=MODES[mode])
                ax.fill_between(centers, np.maximum(values - sigma, 0), values + sigma,
                                step="mid", color=colors[mode], alpha=0.13, linewidth=0)
                count = metrics["source_metrics"][mode][source]
                if mode == "corrected_pre_roi":
                    ax.text(0.97, 0.97,
                            f"full ROI MC: {count['full_roi_mc_count']:,}\n"
                            f"pre-selection: {count['pre_efficiency']:,} / {count['generated']:,}",
                            transform=ax.transAxes, va="top", ha="right", fontsize=8,
                            color=t["ink2"])
                if source == "b8_hep_nu_nr" and mode == "corrected_pre_roi":
                    study = metrics["fivefold_study"]
                    ax.text(0.97, 0.76,
                            f"Independent control: seed {study['seed']}, "
                            f"{study['statistics_factor']}x budget",
                            transform=ax.transAxes, va="top", ha="right", fontsize=8,
                            color=t["ink2"])
            ax.set(xlabel="u = (log10 S2c - ER mean) / ER width",
                   ylabel="Unit-area probability / bin",
                   title=f"{label}: Fig. 4 shape, S1c 3-40 phd")
            ax.set_xlim(-6, 4)
            legend(ax, t, fontsize=7.8, loc="lower right")

            delta = "difference_corrected_minus_legacy"
            unc = "difference_mc_sigma"
            d = archive[f"shapes/{source}/u/top/{delta}"][:]
            ds = archive[f"shapes/{source}/u/top/{unc}"][:]
            difference.axhline(0, color=t["ink2"], lw=0.8, zorder=1)
            difference.fill_between(centers, d - ds, d + ds, step="mid",
                                    color=t["orange"], alpha=0.17, linewidth=0,
                                    label="Matched-budget +/-1 MC sigma")
            difference.step(centers, d, where="mid", color=t["orange"], lw=2,
                            label="Corrected - legacy")
            difference.set(xlabel="u", ylabel="Unit-area difference",
                           title="Shape difference, not a physical-rate change")
            difference.set_xlim(-6, 4)
            if source == "b8_hep_nu_nr":
                fd = archive[f"fivefold/shapes/{source}/u/top/{delta}"][:]
                fds = archive[f"fivefold/shapes/{source}/u/top/{unc}"][:]
                difference.step(centers, fd, where="mid", color=t["ink2"], lw=1.3,
                                ls="--", label="Independent fivefold difference")
                difference.fill_between(centers, fd - fds, fd + fds, step="mid",
                                        color=t["ink2"], alpha=0.07, linewidth=0)
            legend(difference, t, fontsize=7.7, loc="lower right")

        legacy_validation = metrics["validation4"]["legacy_total"]
        corrected_validation = metrics["validation4"]["corrected_pre_roi"]
        footer = (
            "Validation 4 (same exclusive-total validator): "
            f"legacy {'PASS' if legacy_validation['acceptance']['pass'] else 'FAIL'}, "
            f"corrected {'PASS' if corrected_validation['acceptance']['pass'] else 'FAIL'}. "
            "Controls: active reproduction, all non-NR samples, normalization and full histograms PASS.\n"
            f"Corrected-minus-legacy total RMS {corrected_validation['residuals_dex']['total']['rms'] - legacy_validation['residuals_dex']['total']['rms']:+.4f} dex; "
            f"total sub-floor residual {corrected_validation['tail_abs_events']['total'] - legacy_validation['tail_abs_events']['total']:+.4f} events. "
            f"B8 spectral-integral fraction {metrics['b8_integral_fractions']['legacy_total']:.3f} -> {metrics['b8_integral_fractions']['corrected_pre_roi']:.3f}."
        )
        fig.text(0.5, 0.007, footer, fontsize=8.2, color=t["ink2"],
                 ha="center", va="bottom", linespacing=1.65)
        fig.tight_layout(rect=(0.02, 0.058, 0.99, 0.955), h_pad=3.0, w_pad=2.3)
        fig.savefig(args.output, facecolor=t["surface"], bbox_inches="tight")

    print(f"PASS: wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
