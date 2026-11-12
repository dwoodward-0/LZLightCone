"""Render the declared complete ER staged-efficiency comparison.

Modes share blue/orange identity: legacy is solid blue, staged is dashed
orange. Unit-shape bands reflect the conservative per-mode multinomial
errors summed in quadrature. The table separates simulation counts from
fixed Table-I physical means. All three observable shapes are shown for
all eight ER sources; per-S1-panel shape projections remain in the HDF5.
"""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h5py
import numpy as np

from validate_backgrounds import (
    CATEGORY_MAP, EXCLUSIVE_CATEGORY_KEYS, PANEL_S1C, load_band,
    load_reference, project)

LIGHT = {"surface": "#fcfcfb", "ink": "#242422", "muted": "#676b73",
         "grid": "#deded9", "legacy": "#2a78d6", "staged": "#eb6834",
         "G": "#7b3e9d", "reference": "#424542"}
DARK = {"surface": "#1a1a19", "ink": "#f2f2ef", "muted": "#b3b5b8",
        "grid": "#434541", "legacy": "#3987e5", "staged": "#d95926",
        "G": "#9b6ab2", "reference": "#c0c2bf"}
OLD = "legacy_roi_only"
NEW = "staged_conditional"
LABELS = {
    "214Pb_betas": "214Pb", "kr85_ar39_det_gammas": "85Kr+39Ar+det.gamma",
    "solar_nu_er": "solar-nu ER", "pb212_po218_betas": "212Pb+218Po",
    "tritium_c14_betas": "3H+14C", "xe136_2nubb": "136Xe 2nu-bb",
    "xe124_dec": "124Xe DEC", "xe127_xe125_ec": "127Xe+125Xe EC",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--response", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--er-efficiency-treatment", required=True,
                    choices=["staged_conditional", "legacy_no_energy_efficiency"])
    ap.add_argument("--dark", action="store_true")
    args = ap.parse_args()
    c = DARK if args.dark else LIGHT
    with h5py.File(args.input) as f:
        summary = json.loads(f["summary"][()])
        if summary["active_er_efficiency_treatment"] != args.er_efficiency_treatment:
            raise ValueError("Plot mode differs from archived active ER decision")
        curves = {key: f["efficiency_inputs/" + key][:] for key in f["efficiency_inputs"]}
        forward = f["forward_checks/fixed_energy"][:]
        forward_columns = [str(x) for x in f["forward_checks/fixed_energy"].attrs["columns"]]
        shapes = {coordinate: {source: {} for source in summary["er_sources"]}
                  for coordinate in ("s1c", "log10_s2c", "u")}
        for source in summary["er_sources"]:
            for coordinate in shapes:
                # This overview sums both Fig.-4 S1c panels. The archive also
                # retains individual top/bottom projections for every source.
                for panel in ("top", "bottom"):
                    group = f[f"shapes/{source}/{coordinate}/{panel}"]
                    for key in group:
                        if key == "edges":
                            shapes[coordinate][source]["edges"] = group[key][:]
                        else:
                            shapes[coordinate][source][key] = (
                                shapes[coordinate][source].get(key, 0) + group[key][:])
                assert len(shapes[coordinate][source]["edges"]) == len(
                    shapes[coordinate][source]["legacy_unit_mass_per_bin"]) + 1
        old_ess = {s: summary["source_metrics"][OLD][s]["full_histogram_effective_sample_size"]
                   for s in summary["er_sources"]}
        new_ess = {s: summary["source_metrics"][NEW][s]["full_histogram_effective_sample_size"]
                   for s in summary["er_sources"]}
        forward_summary = summary["forward_checks"]
        validation = summary["validation4"]
        metrics = summary["source_metrics"]
        active_nr = f[NEW].attrs["efficiency_treatment"]
        assert f[OLD].attrs["efficiency_treatment"] == active_nr
        models = {mode: {} for mode in (OLD, NEW)}
        for mode in models:
            for source in (summary["er_sources"] + CATEGORY_MAP["b8"]
                           + ["accidental_coincidences"]):
                group = f[f"{mode}/{source}"]
                models[mode][source] = (
                    group["s1c"][:], group["log10_s2c"][:],
                    group["weight"][:])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({
        "font.size": 8.3, "axes.labelcolor": c["ink"],
        "axes.edgecolor": c["muted"], "xtick.color": c["muted"],
        "ytick.color": c["muted"], "text.color": c["ink"],
        "axes.titlecolor": c["ink"], "figure.facecolor": c["surface"],
        "axes.facecolor": c["surface"], "savefig.facecolor": c["surface"],
        "legend.facecolor": c["surface"], "legend.edgecolor": c["surface"],
    })
    band = load_band(args.response)
    reference = load_reference(args.reference)
    total_ref = reference[("top", "total_bkg")]
    u0 = float(total_ref[0, 1])
    du = float(np.median(np.diff(total_ref[:, 1])))
    edges = u0 + du * np.arange(len(total_ref) + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    total_physical = {}
    for mode in (OLD, NEW):
        projections = [project(models[mode], CATEGORY_MAP[key], band, edges)
                       for key in EXCLUSIVE_CATEGORY_KEYS]
        total_physical[mode] = {
            panel: sum((values[panel] for values in projections),
                       start=np.zeros(len(centers)))
            for panel in PANEL_S1C}
    total_display_variance = {mode: {panel: np.zeros(len(centers)) for panel in PANEL_S1C}
                              for mode in (OLD, NEW)}
    for mode in total_display_variance:
        for source, (s1c, y, w) in models[mode].items():
            mu, sigma = band(s1c)
            u = (y - mu) / sigma
            for panel, (lo, hi) in PANEL_S1C.items():
                keep = (s1c > lo) & (s1c <= hi)
                total_display_variance[mode][panel] += np.histogram(
                    u[keep], bins=edges, weights=w[keep] ** 2)[0]

    fig = plt.figure(figsize=(20, 34), constrained_layout=False)
    gs = fig.add_gridspec(11, 4, left=.055, right=.98, bottom=.04, top=.956,
                          hspace=.80, wspace=.30, height_ratios=[1, 1, 1.2, 1.25, 1, 1, 1, 1, 1, 1, .2])
    fig.suptitle("Electron-recoil staged selection: full twelve-source A/B",
                 fontsize=19, fontweight="bold", y=.985)

    ax = fig.add_subplot(gs[0, :2])
    for name, color, label, style in (
            ("trigger_B", c["legacy"], "B: trigger", "-"),
            ("cumulative_s1_Y", c["staged"], "Y: trigger + native coincidence + S1c>3", "--"),
            ("G_overlap_then_G15", c["G"], "G: SS/analysis cumulative, G(15) plateau", "-.") ):
        ax.plot(curves[name][:, 0], curves[name][:, 1], color=color, lw=2, ls=style, label=label)
    ax.plot(curves["independent_total_2d_t4"][:, 0],
            curves["independent_total_2d_t4"][:, 1], color=c["reference"],
            lw=1.2, ls=":", label="t4 total 2D (independent diagnostic)")
    ax.set(title="Public stages — each has its own role; no cumulative-curve products",
           xlabel="Sampled true electron-recoil energy [keV]", ylabel="Efficiency",
           xlim=(0, 20), ylim=(0, 1.05))
    ax.legend(loc="lower left", fontsize=8)

    ax = fig.add_subplot(gs[0, 2:])
    ratio = curves["conditional_G_over_Y"]
    keep = (ratio[:, 0] >= .5) & (ratio[:, 1] > 0)
    ax.plot(ratio[keep, 0], ratio[keep, 1], color=c["staged"], lw=2,
            label="C(E) = G(E) / Y(E), applied once after S1c>3")
    ax.set(title="Conditional SS/analysis acceptance at retained true energy",
           xlabel="Sampled true electron-recoil energy [keV]", ylabel="G / Y", xlim=(.5, 20))
    ax.legend(loc="upper left")
    ax.text(.99, .04, f"max raw C = {summary['input_diagnostics']['max_raw_conditional']:.6f};  "
                      f"clipping = {summary['input_diagnostics']['max_conditional_clipping']:.1g}",
            transform=ax.transAxes, ha="right", color=c["muted"], fontsize=8)

    indices = [forward_columns.index(n) for n in (
        "energy_kevee", "published_cumulative_Y", "predicted_cumulative_Y",
        "published_total_2d", "predicted_total_2d", "nest_full_roi_binomial_sigma")]
    e, yref, ymodel, t4ref, t4model, roi_sigma = [forward[:, i] for i in indices]
    ax = fig.add_subplot(gs[1, :2])
    ax.plot(e, yref, color=c["reference"], lw=2, label="Published cumulative Y")
    ax.plot(e, ymodel, color=c["legacy"], lw=2, ls="--",
            label="B × simulated coincidence + S1c>3")
    ax.set(title="Cumulative S1-stage forward check (not a fit)",
           xlabel="True electron-recoil energy [keV]", ylabel="Cumulative efficiency", xlim=(0, 20))
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, .22))
    ax.text(.99, .05, f"Fixed-energy MC: {forward_summary['statistics_per_energy']} per point\n"
                      f"max |difference| for Y≥0.01 = "
                      f"{forward_summary['s1_max_abs_difference_for_y_at_least_0p01']:.4f}; "
                      f"{forward_summary['n_low_y_grid_points']} low-Y points excluded",
            transform=ax.transAxes, ha="right", color=c["muted"], fontsize=8)

    ax = fig.add_subplot(gs[1, 2:])
    ax.plot(e, t4ref, color=c["reference"], lw=2, label="Independent published t4 total 2D")
    ax.plot(e, t4model, color=c["staged"], lw=2, ls="--",
            label="Forward staged chain: B × (G/Y) × simulated ROI")
    mask = e > 0
    c_value = np.interp(e[mask], ratio[:, 0], ratio[:, 1])
    trigger = np.interp(e[mask], curves["trigger_B"][:, 0], curves["trigger_B"][:, 1])
    ax.fill_between(e[mask], np.maximum(0, t4model[mask] - trigger * c_value * roi_sigma[mask]),
                    np.minimum(1, t4model[mask] + trigger * c_value * roi_sigma[mask]),
                    color=c["staged"], alpha=.16, label="Response MC ±1σ only")
    ax.set(title="Independent total-2D consistency diagnostic — never multiplied in",
           xlabel="True electron-recoil energy [keV]", ylabel="Total efficiency", xlim=(0, 20))
    ax.legend(loc="upper left")
    ax.set_title(
        "Independent total-2D consistency diagnostic — never multiplied in\n"
        f"RMS difference {forward_summary['total_2d_rms_difference']:.4f}; "
        f"max {forward_summary['total_2d_max_abs_difference']:.4f}; "
        f"no published covariance (precision allowance "
        f"{forward_summary['central_precision_allowance']:.2f} + 2 MC sigma)",
        fontsize=10)

    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")
    headings = ["ER source", "Trigger", "S1 stage", "G/Y", "ROI MC: old → staged",
                "Display MC: old → staged", "Full ROI ESS: old → staged",
                "Table-I physical mean (both)"]
    table_rows = []
    for source in summary["er_sources"]:
        old, new = metrics[OLD][source], metrics[NEW][source]
        row = [LABELS[source], str(new["n_after_trigger"]),
               str(new["n_after_native_coincidence_and_s1"]),
               str(new["n_after_conditional_g_over_y"]),
               f"{old['full_roi_mc_count']} → {new['full_roi_mc_count']}",
               f"{old['retained_display_count']} → {new['retained_display_count']}",
               f"{old_ess[source]:.0f} → {new_ess[source]:.0f}",
               f"{new['table1_central_physical_count']:.1f} (unchanged)"]
        table_rows.append(row)
    table = ax.table(cellText=table_rows, colLabels=headings, loc="center",
                     cellLoc="center", colLoc="center", colWidths=[.14, .065, .075, .065, .14, .15, .155, .17])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.35)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor(c["grid"])
        cell.get_text().set_color(c["ink"])
        cell.set_facecolor(c["surface"])
        if row == 0:
            cell.set_facecolor(c["grid"])
            cell.get_text().set_fontweight("bold")
    ax.set_title("Simulation stage counts (staged mode); changing multiplicity ≠ changing physical rate",
                 pad=12, loc="left", fontsize=12, fontweight="bold")
    ax.text(0, -.08, "All twelve physical Table-I means are unchanged. NR sources and observable-space accidentals are bit-identical between modes. "
                     "Independent per-source streams improve A/B isolation; they do not give event-level common random numbers after rejection.",
            color=c["muted"], fontsize=9, va="top")

    for col, panel in enumerate(("top", "bottom")):
        ax = fig.add_subplot(gs[3, 2 * col:2 * col + 2])
        reference_values = load_reference(args.reference)[(panel, "total_bkg")]
        ref_centers = reference_values[:, 0]
        ref_levels = 10 ** reference_values[:, 3]
        ax.step(ref_centers, ref_levels, where="mid", lw=1.3, color=c["reference"],
                label="LZ WS2024 Fig. 4 total")
        for mode, key, style, label in ((OLD, "legacy", "-", "Legacy ER: explicit full ROI only"),
                                         (NEW, "staged", "--", "Staged conditional ER")):
            values = total_physical[mode][panel]
            sigma = np.sqrt(total_display_variance[mode][panel])
            ax.step(centers, np.maximum(values, 1e-4), where="mid", lw=2, ls=style,
                    color=c[key], label=label)
            ax.fill_between(centers, np.maximum(sigma * 0 + 1e-4, values - sigma),
                            values + sigma, step="mid", color=c[key], alpha=.16,
                            label=f"{label}: model MC ±1σ")
        val = validation[mode]
        ax.set(yscale="log", ylim=(3e-3, 400),
               title=f"Exclusive 12-source Fig.-4 total: {panel} S1c panel",
               xlabel="u = [log10 S2c − mu_ER(S1c)] / sigma_ER(S1c)",
               ylabel="Events / bin", xlim=(-6, 4))
        ax.legend(loc="upper right", fontsize=7)
        text = [f"{m}: RMS {validation[m]['residuals_dex']['total']['rms']:.4f} dex; "
                f"tail {validation[m]['tail_abs_events']['total']:.3f} events; "
                f"six gates {'PASS' if validation[m]['acceptance']['pass'] else 'FAIL'}"
                for m in (OLD, NEW)]
        ax.text(.02, .03, "\n".join(text), transform=ax.transAxes, color=c["muted"], fontsize=8)

    coordinate_titles = {"s1c": "S1c [phd]", "log10_s2c": "log10(S2c / phd)",
                         "u": "Fig.-4 coordinate u"}
    for row, coordinate in enumerate(("s1c", "log10_s2c", "u")):
        for index, source in enumerate(summary["er_sources"]):
            axis_col, axis_row = index // 2, 4 + 2 * row + index % 2
            ax = fig.add_subplot(gs[axis_row, axis_col])
            shape = shapes[coordinate][source]
            center = .5 * (shape["edges"][:-1] + shape["edges"][1:])
            difference = shape["difference_corrected_minus_legacy"]
            sigma = shape["difference_mc_sigma"]
            ax.axhline(0, color=c["muted"], lw=.7)
            ax.fill_between(center, -sigma, sigma, step="mid", color=c["staged"], alpha=.18,
                            label="±√(σ_legacy² + σ_staged²)")
            ax.step(center, difference, where="mid", lw=1.6, color=c["staged"],
                    label="Staged − legacy, full unit-area histogram")
            if coordinate == "u":
                ax.step(center, shape["legacy_unit_mass_per_bin"] / 10,
                        where="mid", color=c["legacy"], lw=1, ls=":",
                        label="Legacy unit shape / 10 (reference only)")
            ax.set_title(f"{LABELS[source]}: Δ{coordinate_titles[coordinate]}", fontsize=10)
            ax.set_xlabel(coordinate_titles[coordinate])
            ax.set_ylabel("Unit-shape difference")
            if index == 0:
                ax.legend(fontsize=6.6, loc="best")
    for ax in fig.axes:
        if ax.axison:
            ax.tick_params(length=3)
            ax.grid(axis="y", color=c["grid"], lw=.5, alpha=.65)
            ax.set_axisbelow(True)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
    fig.text(.055, .019,
             f"Common NR treatment: {active_nr}. Each ER source uses its own source-seeded stream. Full MC histograms determine shapes; display samples are capped and separately Table-I normalized.\n"
             "Applied once: trigger → NEST coincidence + S1c>3 → conditional G/Y using retained true energy → S1c<80, S2c≥493, log10(S2c)<4.5. "
             "The energy-conditional approximation does not reproduce unpublished S1c/S2c/position dependence.\n"
             "Overview observable shapes sum both S1c panels; per-panel values, finite-MC diagnostics and all six unchanged Validation 4 gates are in the declared HDF5 archive. "
             "No claim of physical rate change or improved Fig.-4 agreement is made.",
             color=c["muted"], fontsize=9, va="top")
    fig.savefig(args.output, dpi=130)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
