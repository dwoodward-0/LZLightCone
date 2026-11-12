"""Standalone 40 GeV/c2 WS2024 Fig.-4 comparison at an exact input limit.

The model is absolutely normalized from NEST and then multiplied only by the
fixed cross-section ratio (LZ's tabulated observed power-constrained 90% CL
40 GeV/c2 limit divided by NEST's 1e-36 cm2 native reference). This is a
comparison convention at 3.3 WS2024 tonne-years, not a fit, an inference of
LZ's displayed amplitude, or a power-constraint calculation. No background
one-event/bin dex gate, background normalization, or reference-zero filling
is used. Missing LZ sub-floor bins are missing, not zero.

The separate --metrics JSON carries all finite-common-bin and absent-reference
sums. With a .json --output, only that declared metric output is written; with
.png, the figure and optional explicitly supplied scratch metrics are written.
The model error band is finite MC only, never LZ's background band.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_backgrounds import PANEL_S1C, load_band, load_reference
from wimp_native import REFERENCE_SIGMA_CM2

EXPECTED_OBSERVED_LIMIT = 2.1816833824484827e-48
EXPECTED_MEDIAN = 5.0587223537493065e-48
MODEL_COLOR = "#D55E00"
LZ_COLOR = "#7B61A8"


def palette(dark):
    return {
        "surface": "#1a1a19" if dark else "#fcfcfb",
        "fg": "#f2f2ef" if dark else "#20201e",
        "muted": "#b4b4ad" if dark else "#62625c",
        "grid": "#484844" if dark else "#ddddda",
    }


def comparison(model_path, reference_path, sensitivity_path,
               expected_decisions=None, require_full_grid=False):
    with open(sensitivity_path) as f:
        rows = list(csv.DictReader(f))
    row = next(x for x in rows if float(x["mass_gev"]) == 40.0)
    observed = float(row["limit"])
    median = float(row["median"])
    if observed != EXPECTED_OBSERVED_LIMIT:
        raise ValueError(f"The display must use the exact 40 GeV observed limit: {observed}")
    if median != EXPECTED_MEDIAN or observed == median:
        raise ValueError("The observed-limit/median distinction in SI reference is invalid")
    ref = load_reference(reference_path)
    counts = {p: len(ref.get((p, "wimp_40gev"), [])) for p in PANEL_S1C}
    # This checks the actual digitized table: eight upper-panel levels are
    # present and nine would be clipped below the figure floor. A discarded
    # transcription suggested seven; only the file determines the comparison.
    if counts != {"top": 8, "bottom": 3}:
        raise ValueError(f"Unexpected digitized WIMP reference support: {counts}")

    with h5py.File(model_path, "r") as f:
        g = f["ws2024/mass_40gev"]
        for name, expected in (expected_decisions or {}).items():
            actual = g.attrs.get(name)
            if actual != expected:
                raise ValueError(
                    f"WIMP decision {name}: expected {expected!r}, found {actual!r}")
        if require_full_grid:
            groups = [name for name in f["ws2024"] if name.startswith("mass_")]
            actual_masses = np.array([
                f[f"ws2024/{name}"].attrs["mass_gev"] for name in groups])
            published_masses = np.array([float(x["mass_gev"]) for x in rows])
            if (len(groups) != 26 or int(f.attrs.get("n_generated_masses", -1)) != 26
                    or not np.array_equal(np.sort(actual_masses), np.sort(published_masses))):
                raise ValueError(
                    f"Production validation requires exactly the 26 published masses; "
                    f"found {len(groups)} groups")
            if "diagnostics/representative_vs_annual_40gev" not in f:
                raise ValueError("Full signal archive lacks the 40 GeV annual control")
        if not f.attrs.get("completed", False):
            raise ValueError("The WIMP archive was not completed")
        if g.attrs["exposure_tonne_years"] != 3.3 or g.attrs["mass_gev"] != 40.0:
            raise ValueError("This figure requires 40 GeV/c2 and 3.3 tonne-years")
        if g.attrs["efficiency_stage"] != "pre_roi_ss_data_cuts":
            raise ValueError("Incorrect ROI-inclusive efficiency in signal model")
        if f.attrs["comparison_limit_40gev_cm2"] != observed:
            raise ValueError("Input SI limit changed between model and validation")
        ref_sigma = float(g.attrs["reference_sigma_cm2"])
        if ref_sigma != REFERENCE_SIGMA_CM2:
            raise ValueError("Native reference cross-section is not 1e-36 cm2")
        hist_mc = g["hist_mc_counts"][:]
        if int(hist_mc.sum()) != int(g.attrs["n_full_roi_mc"]):
            raise AssertionError("Full integer histogram != retained event count")
        if not np.isclose(hist_mc.sum(), g.attrs["n_full_roi_mc"], atol=0, rtol=0):
            raise AssertionError("Noninteger/thinned histogram used for normalization")
        # One preparation in the standard representative-day model. A single
        # constant w0 makes this an independent exact check of the archived
        # finite-population variance, before display thinning.
        if g.attrs["time_mode"] == "representative_day" and g.attrs["n_seeded_blocks"] == 1:
            w0 = g.attrs["preselection_yield_at_reference_sigma"] / g.attrs["n_native_generated"]
            mc_var = hist_mc * (1 - hist_mc / g.attrs["n_native_generated"]) * w0**2
            if not np.allclose(mc_var, g["hist_mc_weight_variance"][:], rtol=1e-12, atol=0):
                raise AssertionError("Full multinomial MC variance does not reproduce n*p*(1-p)")
        edges = g["u_edges"][:]
        # Band geometry itself was produced with the validated ER response.
        # Reference edges are corrected without independently rounding CSV.
        total_edges = ref[("top", "total_bkg")][:, 1]
        nominal_edges = float(total_edges[0]) + float(np.median(np.diff(total_edges))) * np.arange(len(total_edges) + 1)
        if not np.allclose(edges, nominal_edges, rtol=0, atol=1e-12):
            raise AssertionError("Model/reprojection and LZ u-grid differ")
        centers = (edges[:-1] + edges[1:]) / 2
        out_edges = []
        data = {"model_path": str(model_path),
                "mass_gev": 40.0,
                "campaign": "WS2024 only",
                "exposure_tonne_years": 3.3,
                "native_reference_sigma_cm2": ref_sigma,
                "displayed_sigma_cm2": observed,
                "published_median_not_used_cm2": median,
                "native_rate_events_per_kg_day_at_native_reference_sigma": g.attrs["native_rate_events_kg_day"],
                "absolute_roi_events_at_native_reference_sigma": g.attrs["absolute_roi_yield_at_reference_sigma"],
                "absolute_roi_events_at_displayed_sigma": g.attrs["absolute_roi_yield_at_reference_sigma"] * observed / ref_sigma,
                "full_roi_mc_events": int(hist_mc.sum()),
                "mc_relative_sigma_at_native_reference": float(np.sqrt(g.attrs["roi_yield_mc_variance"]) / g.attrs["absolute_roi_yield_at_reference_sigma"]),
                "normalization_outputs": {
                    "direct_accepted_closure": g.attrs["accepted_yield_direct_closure"],
                    "mean_rate_acceptance_closure": g.attrs["roi_yield_mc_closure"],
                },
                "reference_bin_counts": counts,
                "reference_missing_bins_are_missing_not_zero": True,
                "reference_uncertainty": "No LZ background band is treated as WIMP uncertainty",
                "display_semantics": "Fixed LZ observed-power-constrained limit input at WS2024 exposure; no fit/tuning",
                "comparison_acceptance": "Descriptive go/no-go checkpoint, no agreement-improvement target",
                "panels": {}}
        for panel in PANEL_S1C:
            pg = g[f"fig4_{panel}"]
            model = pg["counts"][:] * observed / ref_sigma
            var = pg["mc_variance"][:] * (observed / ref_sigma) ** 2
            lz = np.full(len(centers), np.nan)
            for row_ in ref.get((panel, "wimp_40gev"), []):
                index = int(np.argmin(np.abs(centers - row_[0])))
                # Match validate_backgrounds: CSV-u offsets, u_lo and u_hi
                # are independently rounded, whereas the model uses the
                # reconstructed arithmetic grid; 0.003 measured here.
                if abs(centers[index] - row_[0]) > 0.01:
                    raise AssertionError("WIMP reference bin is not aligned to model grid")
                lz[index] = 10 ** row_[3]
            common = np.isfinite(lz)
            model_comparable = np.where(common, model, 0.0)
            residual = np.full(len(centers), np.nan)
            residual[common] = model[common] - lz[common]
            sigma = np.sqrt(var)
            sr = {"reference_present_bins": int(common.sum()),
                  "reference_present_centers": centers[common].tolist(),
                  "model_at_reference_bins": model[common].tolist(),
                  "lz_at_reference_bins": lz[common].tolist(),
                  "model_sigma_at_reference_bins": sigma[common].tolist(),
                  "signed_model_minus_lz_at_reference_bins": residual[common].tolist(),
                  "model_sum_at_reference_bins": float(model[common].sum()),
                  "lz_sum_at_reference_bins": float(lz[common].sum()),
                  "signed_sum_model_minus_lz_at_reference_bins": float((model[common] - lz[common]).sum()),
                  "model_all_grid_sum": float(model.sum()),
                  "model_reference_absent_sum": float(model[~common].sum()),
                  "model_reference_absent_mc_sigma": float(np.sqrt(var[~common].sum())),
                  "reference_absent_bin_centers": centers[~common].tolist(),
                  "reference_absent_model_counts": model[~common].tolist(),
                  "model_peak_u": float(centers[np.argmax(model)]),
                  "lZ_peak_u": float(centers[np.nanargmax(lz)]),
                  "peak_model_minus_lz_u": float(centers[np.argmax(model)] - centers[np.nanargmax(lz)]),
                  "max_absolute_residual_at_reference_bins": float(np.max(np.abs(residual[common]))),
                  "reference_level_min_max": [float(np.min(lz[common])), float(np.max(lz[common]))],
                  "model_level_min_max_at_reference_bins": [float(np.min(model[common])), float(np.max(model[common]))],
                  "bins": [{"u": float(centers[i]), "model": float(model[i]),
                            "model_mc_sigma": float(sigma[i]),
                            "reference": float(lz[i]) if common[i] else None,
                            "model_minus_reference": float(residual[i]) if common[i] else None}
                           for i in range(len(centers))]}
            data["panels"][panel] = sr
        if "diagnostics/representative_vs_annual_40gev" in f:
            data["annual_response_check"] = json.loads(f["diagnostics/representative_vs_annual_40gev/summary"][()])
    return data


def render(data, out_path, dark=False):
    c = palette(dark)
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.4),
                             height_ratios=[2.65, 1.0], sharex="col")
    plt.subplots_adjust(hspace=0.08, wspace=0.22)
    fig.patch.set_facecolor(c["surface"])
    for col, panel in enumerate(PANEL_S1C):
        d = data["panels"][panel]
        bins = d["bins"]
        centers = np.array([b["u"] for b in bins])
        model = np.array([b["model"] for b in bins])
        mc = np.array([b["model_mc_sigma"] for b in bins])
        lz = np.array([np.nan if b["reference"] is None else b["reference"] for b in bins])
        residual = np.array([np.nan if b["model_minus_reference"] is None else b["model_minus_reference"] for b in bins])
        ax, lower = axes[0, col], axes[1, col]
        for a in (ax, lower):
            a.set_facecolor(c["surface"])
            a.tick_params(colors=c["fg"], labelsize=9)
            a.xaxis.label.set_color(c["fg"])
            a.yaxis.label.set_color(c["fg"])
            for spine in a.spines.values():
                spine.set_color(c["grid"])
            a.grid(True, color=c["grid"], lw=0.5, alpha=0.55)
            a.set_axisbelow(True)
        step_kw = dict(where="mid", lw=1.8, zorder=4)
        ax.step(centers, model, color=MODEL_COLOR, ls="--", label="This Work, absolute yield", **step_kw)
        show = model > 0
        ax.fill_between(centers[show], np.maximum(model[show] - mc[show], 0),
                        model[show] + mc[show], step="mid", color=MODEL_COLOR,
                        alpha=0.16, linewidth=0, label="This Work ±1σ MC stat.")
        ax.step(centers, lz, color=LZ_COLOR, label="LZ WS2024 Fig. 4: 40 GeV/c² WIMP", **step_kw)
        ax.set_yscale("log")
        ax.set_ylim(3e-3, max(3.0, np.nanmax(model) * 3.5))
        ax.set_ylabel("events per 0.556-wide u bin")
        ax.set_title("3 < S1c < 40 phd" if panel == "top" else "40 < S1c < 80 phd",
                     color=c["fg"], fontsize=11, fontweight="bold")
        ax.text(0.98, 0.94, f"both: comparable bins {d['model_sum_at_reference_bins']:.4g} / {d['lz_sum_at_reference_bins']:.4g}\n"
                f"model, reference absent: {d['model_reference_absent_sum']:.4g}",
                transform=ax.transAxes, ha="right", va="top", fontsize=8.2, color=c["fg"])
        lower.step(centers, residual, color=MODEL_COLOR, **step_kw)
        finite = np.isfinite(residual) & np.isfinite(mc)
        lower.fill_between(centers[finite], -mc[finite], mc[finite], step="mid",
                           color=MODEL_COLOR, alpha=0.12, linewidth=0)
        lower.axhline(0, color=c["muted"], lw=0.8)
        lower.set_ylabel("This Work − LZ")
        lower.set_xlabel("u = [log10(S2c) − median ER] / ER band width")
        finite_r = residual[np.isfinite(residual)]
        extent = max(1e-4, float(np.max(np.abs(finite_r))) * 1.4)
        lower.set_ylim(-extent, extent)
        lower.text(0.02, 0.06, f"±{extent:.3g} events/bin; MC band only",
                   transform=lower.transAxes, fontsize=8, color=c["muted"])
    axes[0, 0].legend(loc="lower left", fontsize=8.5, facecolor=c["surface"],
                      edgecolor=c["grid"], labelcolor=c["fg"], framealpha=0.9)
    fig.suptitle("WIMP model check: LZ WS2024 Fig. 4 projection, 40 GeV/c²",
                 color=c["fg"], fontsize=15, y=0.975)
    fig.text(0.5, 0.007,
             "WS2024 exposure 3.3 t·yr; σ = 2.1816833824484827e−48 cm² = LZ tabulated observed power-constrained limit.\n"             "Native NEST absolute rate; no fit/rediscovery of amplitude/tuning. Blank LZ bins are sub-floor/missing, not zero. No background uncertainty assigned to the WIMP.",
             ha="center", va="bottom", fontsize=8.5, color=c["muted"])
    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor=c["surface"])
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--models", required=True)
    ap.add_argument("--response", required=True,
                    help="Validated ER band calibration used by the stored projection; verified dependency")
    ap.add_argument("--fig4-reference", required=True)
    ap.add_argument("--sensitivity", required=True)
    ap.add_argument("--output", required=True, help="Exactly one declared .png or .json output")
    ap.add_argument("--metrics", help="Explicit scratch/numerical JSON path, separate from a PNG output")
    ap.add_argument("--require-full-grid", action="store_true",
                    help="Refuse anything but the exact 26 published mass groups and the annual control")
    for name in ("time-mode", "support-treatment", "efficiency-treatment",
                 "nest-parameter-set", "detector-effects", "halo-model",
                 "exposure-structure", "response-model"):
        ap.add_argument(f"--{name}", required=True,
                        help="Active ASTRA decision, checked against the stored model metadata")
    ap.add_argument("--dark", action="store_true", help="Render on the validated dark chart surface")
    args = ap.parse_args()
    output = Path(args.output)
    if output.suffix not in {".json", ".png"}:
        raise ValueError("--output must be exactly one .json or .png file")
    if output.suffix == ".json" and args.dark:
        raise ValueError("--dark is meaningful only for PNG rendering")
    if args.metrics and Path(args.metrics) == output:
        raise ValueError("PNG output and metrics must be distinct paths")
    expected = {name.replace("-", "_"): getattr(args, name.replace("-", "_"))
                for name in ("time-mode", "support-treatment", "efficiency-treatment",
                             "nest-parameter-set", "detector-effects", "halo-model",
                             "exposure-structure", "response-model")}
    data = comparison(args.models, args.fig4_reference, args.sensitivity,
                      expected_decisions=expected,
                      require_full_grid=args.require_full_grid)
    data["decision_checks_pass"] = expected
    data["validated_response_calibration_dependency"] = args.response
    if not Path(args.response).is_file():
        raise FileNotFoundError(args.response)
    mu_probe, sigma_probe = load_band(args.response)(
        np.array([3.5, 10.0, 40.0, 79.5]))
    if np.any(~np.isfinite(mu_probe)) or np.any(~np.isfinite(sigma_probe)) \
            or np.any(sigma_probe <= 0):
        raise ValueError("Response calibration does not supply finite ER bands")
    data["response_band_dependency_check"] = {
        "s1c_phd": [3.5, 10.0, 40.0, 79.5],
        "er_median_log10_s2c": mu_probe.tolist(),
        "er_sigma_log10_s2c": sigma_probe.tolist(),
        "meaning": "Finite-response input check; model counts remain full stored projections",
    }
    data["storage_and_closure_checks_pass"] = True
    data["no_LZ_agreement_fit_or_stated_improvement_gate"] = True
    if output.suffix == ".png":
        render(data, output, args.dark)
    if output.suffix == ".json" or args.metrics:
        out_json = output if output.suffix == ".json" else Path(args.metrics)
        out_json.write_text(json.dumps(data, indent=2) + "\n")
    else:
        out_json = None
    for panel, d in data["panels"].items():
        print(f"{panel}: model comparable={d['model_sum_at_reference_bins']:.8g}, "
              f"LZ={d['lz_sum_at_reference_bins']:.8g}, "
              f"model - LZ={d['signed_sum_model_minus_lz_at_reference_bins']:+.8g} events; "
              f"model all grid={d['model_all_grid_sum']:.8g}, "
              f"ref-absent={d['model_reference_absent_sum']:.8g}; "
              f"peak model - LZ={d['peak_model_minus_lz_u']:+.6f}")
    print("full ROI at display limit:", data["absolute_roi_events_at_displayed_sigma"])
    print("wrote", args.output)
    if out_json is not None:
        print("wrote", out_json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
