"""Bounded 40 GeV/c2 FC integration pilot, with unchanged input archives.

This is not a sensitivity validation or a region optimization. The fixed pilot
region was stated before evaluating its counts: stored histogram-bin centers
with 3 < S1c < 20 phd and u < -2, where u is the existing ER-band coordinate.
The background-only median sensitivity is computed exactly, without toys:
the FC upper limit is a non-decreasing function of the observed count n at
fixed b (checked here over the evaluated counts), and for any non-decreasing f
the median of f(N) equals f(median N). So the median limit is FC_ints evaluated
at the median of Poisson(B_region), and the +-1/2 sigma bands are FC_ints at
the corresponding Poisson count quantiles (smallest n with CDF >= q). The same
holds after the power constraint, since max(limit, floor) is also
non-decreasing. Limits are converted to cross sections with the selected
region's absolute WIMP reference yield. Note: this shortcut is specific to
counting; the 2D binned PLR still needs toys.

No NEST simulation, background fitting, parameter scanning, model regeneration,
ROOT operation, or results/ writing occurs here. The optional FC interpolator
is not used. The report is an exploratory JSON, normally under /tmp.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
import warnings

import h5py
import numpy as np
import scipy
from scipy.stats import norm, poisson

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from validate_backgrounds import TOTAL_SOURCES, load_band

EXPECTED_FC_SHA256 = "2e8a566d3940d0e19550afecf335bc94164d26436545b1b2a00c164023f7296d"
REFERENCE_SIGMA_CM2 = 1e-36
MASS_GEV = 40.0
ALPHA = 0.1
FC_GRID_SIGNAL_EVENTS = 0.005
S1C_UPPER_PHD = 20.0
U_UPPER = -2.0


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_region(signal_path, background_path, response_path):
    """Read full physical histograms; never use capped display weights."""
    expected_edges = (
        np.arange(3.0, 80.5, 1.0),
        np.concatenate(([np.log10(493.0)], np.arange(2.70, 4.51, 0.02))),
    )
    band = load_band(response_path)

    with h5py.File(signal_path, "r") as sf:
        sg = sf["ws2024/mass_40gev"]
        if not sf.attrs.get("completed", False):
            raise ValueError("WIMP archive is not marked complete")
        for key, expected in {
            "mass_gev": MASS_GEV,
            "exposure_tonne_years": 3.3,
            "reference_sigma_cm2": REFERENCE_SIGMA_CM2,
            "efficiency_treatment": "input_pre_roi_curve",
            "efficiency_stage": "pre_roi_ss_data_cuts",
            "nest_parameter_set": "ws2024_nest240",
            "detector_effects": "uniform_1d",
            "response_model": "nest_library",
            "halo_model": "shm_standard",
            "time_mode": "representative_day",
        }.items():
            if sg.attrs.get(key) != expected:
                raise ValueError(f"Signal {key}: expected {expected}, found {sg.attrs.get(key)}")
        signal_hist = sg["hist_counts"][:]
        signal_mc_hist = sg["hist_mc_counts"][:]
        signal_mc_weight_variance = sg["hist_mc_weight_variance"][:]
        if not math.isclose(signal_hist.sum(), sg.attrs["absolute_roi_yield_at_reference_sigma"],
                            rel_tol=1e-12, abs_tol=0):
            raise ValueError("Signal full histogram and absolute ROI yield disagree")
        if signal_mc_hist.sum() != sg.attrs["n_full_roi_mc"]:
            raise ValueError("Signal integer histogram is not the unthinned generation")
        signal_edges = (sg["hist_s1c_edges"][:], sg["hist_y_edges"][:])
        signal_edges = (sg["hist_s1c_edges"][:], sg["hist_y_edges"][:])

    background_hist = np.zeros_like(signal_hist)
    background_mc_hist = np.zeros_like(signal_mc_hist)
    background_mc_weight_variance = np.zeros_like(signal_mc_weight_variance)
    source_records = []
    source_histograms = []
    with h5py.File(background_path, "r") as bf:
        for key, expected in {
            "efficiency_treatment": "input_pre_roi_curve",
            "er_efficiency_treatment": "staged_conditional",
        }.items():
            if bf.attrs.get(key) != expected:
                raise ValueError(f"Background {key}: expected {expected}, found {bf.attrs.get(key)}")
        extra_groups = set(bf.keys()) - set(TOTAL_SOURCES) - {"accidentals_fit"}
        if extra_groups or any(source not in bf for source in TOTAL_SOURCES):
            raise ValueError(f"Background source groups changed: missing/extra = {extra_groups}")
        for source in TOTAL_SOURCES:
            bg = bf[source]
            hist = bg["hist_counts"][:]
            if hist.shape != signal_hist.shape:
                raise ValueError(f"{source}: signal/background histogram shape mismatch")
            background_hist += hist
            mc_hist = bg["hist_mc_counts"][:]
            source_histograms.append((hist, mc_hist))
            background_mc_hist += mc_hist
            background_mc_weight_variance += bg["hist_mc_weight_variance"][:]
            source_records.append({
                "source": source,
                "table1_mean_events": float(bg.attrs["table1_mean"]),
                "table1_sigma_events": float(bg.attrs["table1_sigma"]),
                "full_physical_histogram_sum_events": float(hist.sum()),
                "source_seed": int(bg.attrs["source_seed"]),
            })
        background_edges = (bg["hist_s1c_edges"][:], bg["hist_y_edges"][:])

    for actual, expected in zip(signal_edges + background_edges, expected_edges * 2):
        if not np.array_equal(actual, expected):
            raise ValueError("Histogram geometry differs from the declared corrected-observable ROI")
    if np.any(~np.isfinite(signal_hist)) or np.any(signal_hist < 0):
        raise ValueError("Invalid physical signal histogram")
    if np.any(~np.isfinite(background_hist)) or np.any(background_hist < 0):
        raise ValueError("Invalid physical background histogram")

    centers = tuple(0.5 * (edges[1:] + edges[:-1]) for edges in expected_edges)
    s1_grid, y_grid = np.meshgrid(*centers, indexing="ij")
    mu_er, sigma_er = band(s1_grid)
    if np.any(~np.isfinite(mu_er)) or np.any(sigma_er <= 0):
        raise ValueError("Invalid ER-band interpolation at histogram-bin centers")
    u_grid = (y_grid - mu_er) / sigma_er
    region_mask = (s1_grid > 3.0) & (s1_grid < S1C_UPPER_PHD) & (u_grid < U_UPPER)
    if not np.any(region_mask) or np.count_nonzero(region_mask) == region_mask.size:
        raise ValueError("Unexpected pilot region mask")

    # Region yields are sums of the existing absolute physical hist_counts.
    # P(A)/efficiency/exposure factors must not be applied a second time.
    selected_signal = float(signal_hist[region_mask].sum())
    selected_background = float(background_hist[region_mask].sum())
    if selected_signal <= 0 or selected_background <= 0:
        raise ValueError("The fixed pilot requires nonzero signal and background expectations")
    if not math.isclose(background_hist.sum(),
                        sum(row["table1_mean_events"] for row in source_records),
                        rel_tol=1e-12, abs_tol=0):
        raise ValueError("Background histogram sums no longer match the Table-I means")

    for row, (hist, mc_hist) in zip(source_records, source_histograms):
        row["selected_physical_events"] = float(hist[region_mask].sum())
        row["selected_generated_mc_events"] = int(mc_hist[region_mask].sum())

    return {
        "region_mask": region_mask,
        "selected_signal_reference": selected_signal,
        "selected_background": selected_background,
        "selected_signal_full_mc_events": int(signal_mc_hist[region_mask].sum()),
        "selected_background_full_mc_events": int(background_mc_hist[region_mask].sum()),
        "selected_signal_mc_weight_variance": float(signal_mc_weight_variance[region_mask].sum()),
        "selected_background_hist_mc_weight_variance_archive_convention": float(
            background_mc_weight_variance[region_mask].sum()),
        "full_signal_reference": float(signal_hist.sum()),
        "full_background": float(background_hist.sum()),
        "selected_bin_count": int(np.count_nonzero(region_mask)),
        "source_records": source_records,
        "s1c_centers_in_region_phd": sorted(set(s1_grid[region_mask].tolist())),
        "u_range_in_selected_bin_centers": [float(u_grid[region_mask].min()),
                                             float(u_grid[region_mask].max())],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signal", required=True, type=Path)
    parser.add_argument("--background", required=True, type=Path)
    parser.add_argument("--response", required=True, type=Path)
    parser.add_argument("--fc-source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.suffix != ".json" or not args.output.is_absolute():
        raise ValueError("The exploratory report must be an absolute JSON path, normally under /tmp")
    if not args.output.parent.is_dir():
        raise ValueError("Supply an existing scratch output directory; do not create results/ outputs")

    start = time.perf_counter()
    source = args.fc_source.resolve()
    source_hash = sha256(source)
    if source_hash != EXPECTED_FC_SHA256:
        raise ValueError("FC source is not the reviewed, source-pinned file")
    # The collaborator's unchanged source remains the implementation. Import
    # only its one-file directory, not the original Downloads directory.
    sys.path.insert(0, str(source.parent))
    import FeldmanCousins as fc

    region = load_region(args.signal, args.background, args.response)
    b = region["selected_background"]
    s_reference = region["selected_signal_reference"]
    # Explicit unit conversion: at sigma, s(sigma)=s_ref*(sigma/1e-36).
    # Solve mu_UL=s(sigma_UL), so sigma_UL=1e-36*mu_UL/s_ref.
    signal_events_to_cross_section = REFERENCE_SIGMA_CM2 / s_reference
    example_sigma = 1e-48
    example_signal_events = s_reference * example_sigma / REFERENCE_SIGMA_CM2
    cross_section_from_inverse_example = example_sigma / example_signal_events
    if not math.isclose(cross_section_from_inverse_example * s_reference,
                        REFERENCE_SIGMA_CM2, rel_tol=1e-14, abs_tol=0):
        raise AssertionError("Signal-yield/cross-section inverse conversion does not close")

    controls = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for n_control, b_control, expected in [(9, 3.7, (1.490, 11.595)),
                                                (0, 0.0, (0.0, 2.435))]:
            control = fc.FC_ints(int(n_control), float(b_control), alpha=ALPHA)
            if not np.allclose(control, expected, rtol=0, atol=1e-13):
                raise AssertionError(f"Preserved FC source changed at the tested control: {control}")
            controls.append({"n": n_control, "b": b_control,
                             "interval_signal_events": [float(x) for x in control]})
        if caught:
            raise RuntimeError(f"Unexpected FC-control warnings: {[str(x.message) for x in caught]}")

    # Exact background-only quantiles: no toys. Evaluate FC over every count
    # from 0 to the 97.725% Poisson quantile (+2 for the monotonicity check).
    quantile_probabilities = {
        "minus2sigma": float(norm.cdf(-2.0)), "minus1sigma": float(norm.cdf(-1.0)),
        "median": 0.5,
        "plus1sigma": float(norm.cdf(1.0)), "plus2sigma": float(norm.cdf(2.0)),
    }
    quantile_counts = {k: int(poisson.ppf(q, b)) for k, q in quantile_probabilities.items()}
    n_max = quantile_counts["plus2sigma"] + 2
    upper_by_n = []
    for n in range(n_max + 1):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            lower, upper = fc.FC_ints(n, float(b), alpha=ALPHA)
        if caught:
            raise RuntimeError(f"FC warning at n={n}: {[str(x.message) for x in caught]}")
        upper_by_n.append({"n": n, "fc_lower_signal_events": float(lower),
                           "fc_upper_signal_events": float(upper),
                           "poisson_cdf": float(poisson.cdf(n, b))})
    uppers = np.array([row["fc_upper_signal_events"] for row in upper_by_n])
    if np.any(np.diff(uppers) < 0):
        raise AssertionError("FC upper limit is not non-decreasing in n; median-of-limits shortcut invalid")

    ordinary_mu = {k: float(uppers[n]) for k, n in quantile_counts.items()}
    ordinary_sigma = {k: v * signal_events_to_cross_section for k, v in ordinary_mu.items()}
    # LZ-style 90%-CL power floor: exp(Phi^-1(0.10)) times the
    # background-only median cross section. Apply at assembly, not inside FC.
    pc_multiplier = float(np.exp(norm.ppf(ALPHA)))
    pc_median_floor_mu = ordinary_mu["median"] * pc_multiplier
    pc_median_floor_sigma = pc_median_floor_mu * signal_events_to_cross_section
    # Max in signal-event units, convert once (avoids ulp-level round-trip drift).
    pc_mu = {k: max(v, pc_median_floor_mu) for k, v in ordinary_mu.items()}
    pc_sigma = {k: v * signal_events_to_cross_section for k, v in pc_mu.items()}
    # Probability that a background-only experiment's limit falls below the floor.
    n_below_floor = [row["n"] for row in upper_by_n
                     if row["fc_upper_signal_events"] < pc_median_floor_mu]
    p_floor_binds = float(poisson.cdf(max(n_below_floor), b)) if n_below_floor else 0.0

    summary = {
        "quantile_definition": "Poisson count quantile n_q = smallest n with CDF(n; B) >= q; limit quantile = FC upper limit at n_q (exact by monotonicity).",
        "quantile_probabilities": quantile_probabilities,
        "quantile_counts": quantile_counts,
        "ordinary_fc_upper_limit_signal_events": ordinary_mu,
        "ordinary_cross_section_upper_limit_cm2": ordinary_sigma,
        "power_constrained_cross_section_upper_limit_cm2": pc_sigma,
        "discreteness_note": "Because n is an integer, the median sensitivity is a step function of B; this is real counting-experiment behavior, not noise.",
    }

    report = {
        "scope": "Integration/normalization pilot only; not LZ sensitivity validation, FC coverage validation, optimized-region production, or a declared materialized output.",
        "region": {
            "fixed_before_count_evaluation": True,
            "selection": "Stored histogram-bin centers: 3 < S1c < 20 phd and u < -2; u=(log10(S2c)-mu_ER)/sigma_ER from validated response ER-band quantiles.",
            "bin_center_convention": "Whole existing 77x91 physical histogram bins, including the narrow first log10(S2c) bin from log10(493) to 2.70; no new event generation or partial-bin split.",
            "optimization_performed": False,
            "selected_bin_count": region["selected_bin_count"],
            "s1c_centers_phd": region["s1c_centers_in_region_phd"],
            "selected_u_center_range": region["u_range_in_selected_bin_centers"],
        },
        "coordinate_inputs": {
            "response_path": str(args.response),
            "response_sha256": sha256(args.response),
            "band": "er_ch3t/band_log10s2c_quantiles; median q[2], half-width (q[3]-q[1])/2",
        },
        "physical_expectations": {
            "background_mean_selected_region_events": b,
            "signal_mean_selected_region_at_reference_sigma_events": s_reference,
            "signal_reference_sigma_cm2": REFERENCE_SIGMA_CM2,
            "exposure_tonne_years": 3.3,
            "signal_full_roi_at_reference_sigma_events": region["full_signal_reference"],
            "background_full_roi_events": region["full_background"],
            "background_source_count": len(TOTAL_SOURCES),
            "sources": region["source_records"],
            "physical_rate_note": "Background physical means are fixed Table-I expectations. Background-only counts are Poisson(sum) with fixed mean. Detector-neutron zero physical mean remains zero.",
        },
        "generated_mc_diagnostics_not_physical_rate_uncertainties": {
            "signal_selected_full_mc_count": region["selected_signal_full_mc_events"],
            "background_selected_full_mc_count": region["selected_background_full_mc_events"],
            "signal_selected_hist_mc_weight_variance": region["selected_signal_mc_weight_variance"],
            "background_selected_hist_mc_weight_variance_archive_convention": region["selected_background_hist_mc_weight_variance_archive_convention"],
            "caveat": "Signal variance is the archive's per-bin finite-native-count convention. Background hist_mc_weight_variance uses its documented raw-w^2 convention and is not a physical-rate variance. Neither is inserted as a Gaussian prior, an extra toy width, or a normalization correction in this fixed-background FC pilot.",
        },
        "fc_interface": {
            "source": str(source),
            "imported_module": fc.__file__,
            "sha256": source_hash,
            "direct_calls": True,
            "optional_interpolated_lookup_used": False,
            "alpha": ALPHA,
            "confidence_level": 1.0 - ALPHA,
            "native_search_grid_signal_events": FC_GRID_SIGNAL_EVENTS,
            "output_units": "Continuous expected signal events, not cross sections or integer Poisson count quantiles.",
            "two_controls_confirm_unchanged_source": controls,
            "fc_warnings": 0,
            "independent_coverage_or_algorithm_validation_not_performed": True,
        },
        "background_only_quantiles": {
            "method": "exact; no toys",
            "model": "Background-only N~Poisson(B_region), fixed B",
            "argument": "FC upper limit is non-decreasing in n at fixed b (asserted over n=0..n_max), so median_N f(N) = f(median N); likewise for other quantiles and for max(f, floor).",
            "n_max_evaluated": n_max,
            "fc_limits_by_count": upper_by_n,
        },
        "normalization_closure": {
            "formula": "sigma_UL_cm2 = 1e-36 * mu_UL_signal_events / mu_WIMP_selected_at_1e-36",
            "signal_events_to_cross_section_cm2": signal_events_to_cross_section,
            "example_sigma_cm2": example_sigma,
            "example_signal_events_from_absolute_reference_yield": example_signal_events,
            "inverse_example_sigma_per_signal_event_cm2": cross_section_from_inverse_example,
            "additional_efficiency_or_exposure_factor_applied": False,
            "table_i_or_figure_signal_normalization_applied": False,
            "signal_hist_counts_used_not_display_weights": True,
        },
        "power_constraint": {
            "convention": "max(ordinary_sigma_UL, exp(norm.ppf(0.1))*ordinary_median_sigma_UL)",
            "pc_multiplier": pc_multiplier,
            "ordinary_median_based_floor_sigma_cm2": pc_median_floor_sigma,
            "counts_with_limit_below_floor": n_below_floor,
            "background_only_probability_floor_binds": p_floor_binds,
            "applied_at_result_assembly_not_inside_fc": True,
        },
        "results": summary,
        "input_provenance": {
            "signal_path": str(args.signal),
            "signal_sha256": sha256(args.signal),
            "background_path": str(args.background),
            "background_sha256": sha256(args.background),
            "models_opened_read_only": True,
            "numpy_version": np.__version__,
            "scipy_version": scipy.__version__,
            "python_version": platform.python_version(),
        },
        "elapsed_seconds": time.perf_counter() - start,
        "next_step_not_executed": "Review this pilot before any progressive voxel-region selection, 26-mass FC curve, or PLR production.",
    }
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))
    print(f"\nWrote bounded exploratory report: {args.output}")


if __name__ == "__main__":
    main()
