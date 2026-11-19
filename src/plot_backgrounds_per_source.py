"""Per-source backgrounds review plots (one PNG per Table I source).

Review companion to validate_backgrounds.py: instead of overlaying every
category on one axes, this writes one figure per background source, each
showing the source's band-projected u = (log10 S2c - mu_ER(S1c))/sigma_ER
distribution in both S1c panels (top 3-40 phd, bottom 40-80 phd), against
LZ's digitized Fig. 4 outline where the WS2024 paper draws one.

Sources LZ does not draw individually (212Pb/218Po, solar-nu ER,
tritium+C14, detector neutrons, atm-nu NR) get their own figure with ours
alone, so their shapes can be inspected without a reference.

For the kr85 review figure, the kr85/39Ar/det-gammas and 136Xe rows are
summed into a single "This Work" curve because digitization shows that
LZ's drawn "Other ER betas" outline corresponds to this sum. The separate
136Xe diagnostic figure retains its own contribution and the combined
curve for reference.

Colored shaded bands are +/-1 sigma from finite-MC statistics (the same
sum-of-weights-squared errors previously drawn as error bars). Where LZ's
outline maps one-to-one to a Table I row, a gray band shows its pre-fit
normalization uncertainty instead: one shared multiplicative rate scaling
of every bin, not a per-bin statistical or shape uncertainty.

For 124Xe, the main red curve is our Fig. 4 RMS-optimal LL charge yield
q=0.78; LZ's fitted q=0.70 +/- 0.04 is retained as a green comparison. As
in validation 4, both are analytic S2 rescalings of the stored q=0.65
baseline for LL-dominated S1c > 40 events, not new NEST productions:
dlog10(S2c) = log10(q/0.65). Every other modeling choice is held fixed.
Our value is a figure-matching diagnostic, not a likelihood fit to data.

Probe: writes PNGs to the given --outdir; not a declared output.
"""
import argparse
import csv

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PANEL_S1C = {"top": (3.0, 40.0), "bottom": (40.0, 80.0)}
MAIN_FLOOR = 1.0  # LZ drawn level above which the outline is informative
MODEL_COLOR = "#d62728"  # red; paired with blue below for CVD separation
SUM_COLOR = "#0072b2"
XE124_Q_FIT_COLOR = "#009e73"  # green; three-color palette validated
XE124_Q_LL_BASELINE = 0.65
XE124_Q_LL_FIT = 0.70
XE124_Q_LL_BEST_RMS = 0.78

OTHER_ER_SOURCES = ["kr85_ar39_det_gammas", "xe136_2nubb"]

# Only one-to-one mappings get an LZ normalization band. In particular,
# other_er_betas is the sum of two Table I rows with unknown correlation,
# so summing their sigmas in quadrature would invent an uncertainty model.
LZ_REFERENCE_SOURCE = {
    "pb214_betas": "214Pb_betas",
    "xe124_dec": "xe124_dec",
    "xe127_125_ec": "xe127_xe125_ec",
    "accidentals": "accidental_coincidences",
}

# source id -> (plot title, LZ Fig. 4 series or None, extra overlays)
# extra overlays: (label, [sources to sum], color)
SOURCE_PLOTS = {
    "214Pb_betas": ("214Pb betas (743 +/- 88)", "pb214_betas", []),
    "kr85_ar39_det_gammas": (
        "85Kr + 39Ar + det. gammas + 136Xe (217.6) — aka Other ERs",
        "other_er_betas", []),
    "xe136_2nubb": (
        "136Xe 2nu-beta-beta (55.6 +/- 8.3)", "other_er_betas",
        [("kr85 + 136Xe (drawn outline is this sum)",
          ["kr85_ar39_det_gammas", "xe136_2nubb"], SUM_COLOR)]),
    "pb212_po218_betas": ("212Pb + 218Po betas (62.7 +/- 7.5)", None, []),
    "solar_nu_er": ("Solar-nu ER, RRPA (102 +/- 6)", None, []),
    "tritium_c14_betas": (
        "Tritium + 14C calib. sources, 8:1 (58.3 +/- 3.3)", None, []),
    "xe124_dec": ("124Xe EC, per-line q (19.4 +/- 2.5)", "xe124_dec", []),
    "xe127_xe125_ec": ("127Xe/125Xe EC (3.2 +/- 0.6)", "xe127_125_ec", []),
    "accidental_coincidences": ("Accidental coincidences (2.8 +/- 0.6)",
                                "accidentals", []),
    "b8_hep_nu_nr": ("8B + hep nu NR (0.06 +/- 0.01)", "b8", []),
    "atm_nu_nr": ("Atm-nu NR (0.12 +/- 0.02)", None, []),
    "detector_neutrons": ("Detector neutrons (0.0 +/- 0.2; zero-weight)",
                          None, []),
}


def load_reference(path):
    """(panel, series) -> (u_centers, u_lo, u_hi, log10 levels)."""
    out = {}
    with open(path) as f:
        rows = [r for r in csv.reader(f) if r and not r[0].startswith("#")
                and r[0] != "panel"]
    for panel, series, uc, ulo, uhi, l10 in rows:
        out.setdefault((panel, series), []).append(
            (float(uc), float(ulo), float(uhi), float(l10)))
    return {k: np.array(sorted(v)) for k, v in out.items()}


def load_band(response_path):
    with h5py.File(response_path) as f:
        centers = f["er_ch3t/band_s1c_centers"][:]
        q = f["er_ch3t/band_log10s2c_quantiles"][:]
    mu, sigma = q[2], 0.5 * (q[3] - q[1])

    def band(s1c):
        return np.interp(s1c, centers, mu), np.interp(s1c, centers, sigma)
    return band


def load_models(path):
    models = {}
    with h5py.File(path) as f:
        for src in f:
            if src == "accidentals_fit":
                continue
            g = f[src]
            models[src] = (g["s1c"][:], g["log10_s2c"][:], g["weight"][:])
    return models


def u_hist(models, sources, band, edges, panel, q_ll=None):
    """Weighted u-histogram per panel plus its Poisson variance.

    Within a source every stored sample carries the same weight
    (Table I mean / N_kept), so sum(w^2) per bin = k * w^2 — exactly the
    Poisson variance of the underlying count k, scaled by w^2.

    q_ll applies validation 4's analytic 124Xe LL charge-yield variation:
    for xe124_dec with S1c > 40 only, dlog10(S2c) = log10(q_ll / baseline).
    Other sources are explicitly unaffected, allowing summed-model use."""
    lo, hi = PANEL_S1C[panel]
    h = np.zeros(len(edges) - 1)
    var = np.zeros(len(edges) - 1)
    for src in sources:
        s1c, y, w = models[src]
        mu, sig = band(s1c)
        if q_ll is not None and src == "xe124_dec":
            y = y + np.where(
                s1c > 40.0, np.log10(q_ll / XE124_Q_LL_BASELINE), 0.0)
        u = (y - mu) / sig
        m = (s1c > lo) & (s1c <= hi)
        h += np.histogram(u[m], bins=edges, weights=w[m])[0]
        var += np.histogram(u[m], bins=edges, weights=w[m] ** 2)[0]
    return h, var


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--models", required=True)
    ap.add_argument("--response", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--table1", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    ref = load_reference(args.reference)
    band = load_band(args.response)
    models = load_models(args.models)
    with open(args.table1) as f:
        table1 = {
            row["source"]: (float(row["prefit_mean"]),
                            float(row["prefit_sigma"]))
            for row in csv.DictReader(
                line for line in f if not line.startswith("#"))
        }

    # The shared u grid, taken from the total_bkg rows as in validation 4.
    # Column 1 is the LEFT EDGE; using the center in column 0 here would
    # shift both the model histogram and drawn reference by half a bin.
    # The CSV's lo/hi values are independently rounded, so reconstruct one
    # arithmetic grid rather than uniquing them into 0.001-wide slivers.
    tot = ref[("top", "total_bkg")]
    u0 = float(tot[0, 1])
    du = float(np.median(np.diff(tot[:, 1])))
    edges = u0 + du * np.arange(len(tot) + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    # LZ outline levels on our grid (nearest center), per panel/series
    lz = {}
    for (panel, series), rows in ref.items():
        lv = np.full(len(centers), np.nan)
        for uc, ulo, uhi, l10 in rows:
            i = int(round((uc - centers[0]) / du))
            if 0 <= i < len(centers):
                lv[i] = l10
        lz[(panel, series)] = lv

    for src, (title, series, extras) in SOURCE_PLOTS.items():
        fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                 gridspec_kw={"hspace": 0.12})
        for ax, panel in zip(axes, PANEL_S1C):
            # LZ outline (solid): a digitized smooth expectation. Where its
            # Table I mapping is unambiguous, scale every bin by one shared
            # pre-fit normalization uncertainty. No shape variation is
            # implied, and these gray envelopes must not be treated as
            # independent per-bin errors in a significance calculation.
            if series is not None:
                lv = lz.get((panel, series))
                if lv is not None:
                    ok = np.isfinite(lv)
                    lz_counts = 10 ** lv[ok]
                    ref_source = LZ_REFERENCE_SOURCE.get(series)
                    if ref_source is not None:
                        mean, sigma = table1[ref_source]
                        delta = lz_counts * sigma / mean
                        ax.fill_between(
                            centers[ok], np.maximum(lz_counts - delta, 0.0),
                            lz_counts + delta, step="mid", color="0.15",
                            alpha=0.14, linewidth=0,
                            label="LZ Table I rate ±1σ (correlated)")
                    ax.step(centers[ok], lz_counts, where="mid",
                            color="0.15", lw=1.6,
                            label="LZ WS2024 Fig.4 (2410.17036)")
            # Extra overlays first (e.g. the kr85+136Xe sum). Shaded bands
            # replace error bars but encode the same +/-1 sigma MC stats.
            for label, srcs, color in extras:
                h, var = u_hist(models, srcs, band, edges, panel)
                nz = h > 0
                sigma_mc = np.sqrt(var[nz])
                ax.fill_between(
                    centers[nz], np.maximum(h[nz] - sigma_mc, 0.0),
                    h[nz] + sigma_mc, step="mid", color=color, alpha=0.12,
                    linewidth=0)
                ax.step(centers[nz], h[nz], where="mid", color=color,
                        lw=1.6, ls="--", label=label)
            # The source itself (dashed, on top) with a shaded MC-stat band.
            # For the primary Other ER review, show only the two-row sum.
            # Other sources plotting is unaffected.
            model_sources = (OTHER_ER_SOURCES
                             if src == "kr85_ar39_det_gammas" else [src])
            h, var = u_hist(models, model_sources, band, edges, panel,
                            q_ll=XE124_Q_LL_BEST_RMS)
            nz = h > 0
            sigma_mc = np.sqrt(var[nz])
            ax.fill_between(
                centers[nz], np.maximum(h[nz] - sigma_mc, 0.0),
                h[nz] + sigma_mc, step="mid", color=MODEL_COLOR, alpha=0.12,
                linewidth=0, label="This Work ±1σ MC stat.")
            model_label = (
                "This Work, q_LL = 0.78 (best Fig.4 RMS; analytic)"
                if src == "xe124_dec" else "This Work")
            ax.step(centers[nz], h[nz], where="mid", color=MODEL_COLOR,
                    lw=1.8, ls="--", label=model_label)
            if src == "xe124_dec" and panel == "bottom":
                q_variations = [
                    (XE124_Q_LL_FIT, XE124_Q_FIT_COLOR, "-.",
                     "This Work, q_LL = 0.70 (LZ fit +/- 0.04; analytic)"),
                ]
                for q, color, linestyle, label in q_variations:
                    h_var, _ = u_hist(
                        models, [src], band, edges, panel, q_ll=q)
                    nz_var = h_var > 0
                    ax.step(centers[nz_var], h_var[nz_var], where="mid",
                            color=color, lw=1.8, ls=linestyle, label=label)
            ax.set_yscale("log")
            ax.set_ylim(2e-3, 2e3)
            ax.set_ylabel("events / bin")
            ax.grid(alpha=0.25, lw=0.5)
            ax.set_title(f"{title}  —  {panel} panel "
                         f"({PANEL_S1C[panel][0]:.0f}<S1c<"
                         f"{PANEL_S1C[panel][1]:.0f})", fontsize=10)
        # Match the published Fig. 4 horizontal range in every component
        # review, so narrow/tiny components retain their visual context.
        axes[1].set_xlim(-6.0, 4.0)
        axes[1].set_xlabel("u = (log10 S2c - mu_ER(S1c)) / sigma_ER(S1c)")
        axes[0].legend(fontsize=8, loc="upper left")
        if src == "xe124_dec":
            axes[1].legend(fontsize=8, loc="upper right")
        if series == "accidentals":
            axes[1].text(0.5, 0.5,
                         "LZ's bottom-panel accidentals outline lies\n"
                         "entirely below the 10^-2 figure floor\n"
                         "(absent from the digitized reference)",
                         transform=axes[1].transAxes, ha="center",
                         fontsize=9, color="0.35", style="italic")
        slug = src.replace("_", "-")
        out = f"{args.outdir}/bkg_{slug}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("wrote", out)


if __name__ == "__main__":
    main()
